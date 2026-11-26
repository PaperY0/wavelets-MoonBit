// 性质测试：不依赖黄金数据，验证数学不变量。
// 核心是“完美重构”全扫描：106 个离散小波 × 9 种模式 × 9 种长度。

///|
test "perfect reconstruction: all discrete wavelets x all modes x lengths" {
  let modes = [
    @wavelets.Mode::Zero,
    @wavelets.Mode::Constant,
    @wavelets.Mode::Symmetric,
    @wavelets.Mode::Periodic,
    @wavelets.Mode::Smooth,
    @wavelets.Mode::Periodization,
    @wavelets.Mode::Reflect,
    @wavelets.Mode::AntiSymmetric,
    @wavelets.Mode::AntiReflect,
  ]
  let lengths = [1, 2, 3, 5, 8, 9, 16, 31, 64]
  for name in @wavelets.wavelist(kind=Some(@wavelets.Discrete)) {
    // dmey 是 Meyer 的 FIR 近似，非严格正交（pywt 自身重构误差 ~0.03），
    // 单独用近似容差测试，不参与精确完美重构扫描。
    if name == "dmey" {
      continue
    }
    let w = @wavelets.Wavelet::new(name)
    for mode in modes {
      for n in lengths {
        // 上游 dwt 对 reflect/antireflect 要求 N>=2，跳过 n=1 的该组合。
        if n < 2 &&
          (
            mode is @wavelets.Mode::Reflect ||
            mode is @wavelets.Mode::AntiReflect
          ) {
          continue
        }
        let x = Array::makei(n, i => {
          (i * 37 % 23).to_double() - 11.0 + 0.001 * i.to_double()
        })
        let (a, d) = @wavelets.dwt(x, w, mode~)
        let r = @wavelets.idwt(Some(a), Some(d), w, mode~)
        let expected_len = if mode is @wavelets.Mode::Periodization {
          2 * ((n + 1) / 2)
        } else {
          2 * a.length() - w.rec_len() + 2
        }
        if r.length() != expected_len {
          fail("recon length \{name}/\{mode.name()}/n=\{n}")
        }
        for i in 0..<n {
          if tabs(r[i] - x[i]) > 1.0e-8 {
            fail(
              "recon \{name}/\{mode.name()}/n=\{n} at \{i}: \{r[i]} != \{x[i]}",
            )
          }
        }
      }
    }
  }
}

///|
test "multilevel roundtrip preserves original samples" {
  let modes = [
    @wavelets.Mode::Symmetric,
    @wavelets.Mode::Periodization,
    @wavelets.Mode::Smooth,
  ]
  for name in ["haar", "db2", "sym4", "coif2", "bior3.5"] {
    let w = @wavelets.Wavelet::new(name)
    for mode in modes {
      for n in [32, 64] {
        let x = Array::makei(n, i => (i * 53 % 31).to_double() - 15.0)
        let coeffs = @wavelets.wavedec(x, w, mode~)
        let r = @wavelets.waverec(coeffs, w, mode~)
        for i in 0..<n {
          if tabs(r[i] - x[i]) > 1.0e-8 {
            fail("multilevel \{name}/\{mode.name()}/n=\{n} at \{i}")
          }
        }
      }
    }
  }
}

///|
test "2d roundtrip preserves original samples" {
  for name in ["haar", "db2", "bior3.5"] {
    let w = @wavelets.Wavelet::new(name)
    for mode in [@wavelets.Mode::Symmetric, @wavelets.Mode::Periodization] {
      for (rows, cols) in [(8, 8), (7, 5)] {
        let m = @wavelets.Matrix::new(
          rows,
          cols,
          Array::makei(rows * cols, i => (i * 29 % 17).to_double() - 8.0),
        )
        let c = @wavelets.dwt2(m, w, mode~)
        let r = @wavelets.idwt2(c, w, mode~)
        for i in 0..<rows {
          for j in 0..<cols {
            if tabs(r.get(i, j) - m.get(i, j)) > 1.0e-8 {
              fail("2d \{name}/\{mode.name()} \{rows}x\{cols} at \{i},\{j}")
            }
          }
        }
        // 多层 2D 往返：仅用 haar（小滤波器，8x8/7x5 支持 level 2）。
        if name == "haar" {
          let dec = @wavelets.wavedec2(m, w, mode~, level=Some(2))
          let r2 = @wavelets.waverec2(dec, w, mode~)
          for i in 0..<rows {
            for j in 0..<cols {
              if tabs(r2.get(i, j) - m.get(i, j)) > 1.0e-8 {
                fail("2d multilevel \{name}/\{mode.name()} at \{i},\{j}")
              }
            }
          }
        }
      }
    }
  }
}

///|
test "orthogonal filter banks: unit energy and double-shift orthogonality" {
  for name in ["haar", "db2", "db8", "sym4", "coif1"] {
    let w = @wavelets.Wavelet::new(name)
    let h = w.dec_lo
    // 单位能量。
    let mut energy = 0.0
    for v in h {
      energy += v * v
    }
    if tabs(energy - 1.0) > 1.0e-12 {
      fail("\{name}: dec_lo energy \{energy} != 1")
    }
    // 双移位正交：sum_k h[k] h[k+2j] = 0 (j != 0)。
    for j in 1..<(h.length() / 2) {
      let mut acc = 0.0
      for k in 0..<(h.length() - 2 * j) {
        acc += h[k] * h[k + 2 * j]
      }
      if tabs(acc) > 1.0e-12 {
        fail("\{name}: double-shift orthogonality violated at j=\{j}")
      }
    }
    // 正交小波：rec_lo 是 dec_lo 的反转。
    let rl = w.rec_lo
    for k in 0..<h.length() {
      if tabs(rl[k] - h[h.length() - 1 - k]) > 1.0e-15 {
        fail("\{name}: rec_lo != reverse(dec_lo)")
      }
    }
  }
}

///|
test "dmey is an approximate (FIR) orthogonal wavelet" {
  // dmey 滤波器能量略大于 1，重构存在 ~0.03 量级误差——与 pywt 完全一致，
  // 这里验证误差落在已知范围内，而非要求精确完美重构。
  let w = @wavelets.Wavelet::new("dmey")
  let mut energy = 0.0
  for v in w.dec_lo {
    energy += v * v
  }
  if tabs(energy - 1.0022448292411643) > 1.0e-9 {
    fail("dmey energy \{energy}")
  }
  for
    mode in [
      @wavelets.Mode::Zero,
      @wavelets.Mode::Symmetric,
      @wavelets.Mode::Periodization,
    ] {
    for n in [1, 16, 64] {
      let x = Array::makei(n, i => (i * 37 % 23).to_double() - 11.0)
      let (a, d) = @wavelets.dwt(x, w, mode~)
      let r = @wavelets.idwt(Some(a), Some(d), w, mode~)
      for i in 0..<n {
        if tabs(r[i] - x[i]) > 0.05 {
          fail(
            "dmey approx recon \{mode.name()}/n=\{n} at \{i}: \{r[i]} != \{x[i]}",
          )
        }
      }
    }
  }
}

///|
test "cwt is linear in the signal" {
  let w = @wavelets.Wavelet::new("mexh")
  let scales = [1.0, 2.0, 4.0]
  let x = Array::makei(24, i => (i * 13 % 11).to_double() - 5.0)
  let y = Array::makei(24, i => (i * 7 % 5).to_double() - 2.0)
  let combo = Array::makei(24, i => 2.0 * x[i] - 0.5 * y[i])
  let (cx, _) = @wavelets.cwt(x, scales, w)
  let (cy, _) = @wavelets.cwt(y, scales, w)
  let (cc, _) = @wavelets.cwt(combo, scales, w)
  for s in 0..<scales.length() {
    for i in 0..<24 {
      let expect = 2.0 * cx.get(s, i) - 0.5 * cy.get(s, i)
      if tabs(cc.get(s, i) - expect) > 1.0e-8 {
        fail("cwt linearity at scale \{s}, sample \{i}")
      }
    }
  }
}

///|
test "single-branch idwt sums to full reconstruction" {
  for name in ["haar", "db2", "sym4", "bior3.5"] {
    let w = @wavelets.Wavelet::new(name)
    let x = Array::makei(16, i => (i + 1).to_double())
    let (a, d) = @wavelets.dwt(x, w)
    let ra = @wavelets.idwt(Some(a), None, w)
    let rd = @wavelets.idwt(None, Some(d), w)
    let full = @wavelets.idwt(Some(a), Some(d), w)
    for i in 0..<full.length() {
      if tabs(ra[i] + rd[i] - full[i]) > 1.0e-9 {
        fail("\{name}: branch additivity at \{i}")
      }
    }
  }
}
