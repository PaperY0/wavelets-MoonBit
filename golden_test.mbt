// 黄金差分测试：全部用例来自 tools/gen_golden.py 导出的 PyWavelets 1.8.0
// 参考值（golden_data_test.mbt）。每个 API 一组，逐用例对照，
// 失败时报告用例序号与小波/模式名，便于定位。

///|
test "dwt matches pywt (910 cases)" {
  let cases = jarr(jget(golden_root, "dwt"))
  for i, case in cases {
    let wname = jstr(jget(case, "w"))
    let mname = jstr(jget(case, "mode"))
    let w = @wavelets.Wavelet::new(wname)
    let mode = @wavelets.Mode::parse(mname)
    let (a, d) = @wavelets.dwt(jnums(jget(case, "x")), w, mode~)
    if !arr_close(a, jnums(jget(case, "a")), 1.0e-9) {
      fail("dwt[\{i}] \{wname}/\{mname}: cA mismatch")
    }
    if !arr_close(d, jnums(jget(case, "d")), 1.0e-9) {
      fail("dwt[\{i}] \{wname}/\{mname}: cD mismatch")
    }
  }
}

///|
test "idwt matches pywt (405 cases, incl. None branches)" {
  let cases = jarr(jget(golden_root, "idwt"))
  for i, case in cases {
    let wname = jstr(jget(case, "w"))
    let mname = jstr(jget(case, "mode"))
    let w = @wavelets.Wavelet::new(wname)
    let mode = @wavelets.Mode::parse(mname)
    let ja = jget(case, "a")
    let jd = jget(case, "d")
    let ca = if ja is Null { None } else { Some(jnums(ja)) }
    let cd = if jd is Null { None } else { Some(jnums(jd)) }
    let r = @wavelets.idwt(ca, cd, w, mode~)
    if !arr_close(r, jnums(jget(case, "r")), 1.0e-9) {
      fail("idwt[\{i}] \{wname}/\{mname}: reconstruction mismatch")
    }
  }
}

///|
test "wavedec/waverec match pywt (104 cases)" {
  let cases = jarr(jget(golden_root, "wavedec"))
  for i, case in cases {
    let wname = jstr(jget(case, "w"))
    let mname = jstr(jget(case, "mode"))
    let w = @wavelets.Wavelet::new(wname)
    let mode = @wavelets.Mode::parse(mname)
    let level = jint(jget(case, "level"))
    let coeffs = @wavelets.wavedec(
      jnums(jget(case, "x")),
      w,
      mode~,
      level=Some(level),
    )
    let expected = jarr(jget(case, "coeffs"))
    if coeffs.length() != expected.length() {
      fail("wavedec[\{i}] \{wname}/\{mname}: level count mismatch")
    }
    for k in 0..<coeffs.length() {
      if !arr_close(coeffs[k], jnums(expected[k]), 1.0e-9) {
        fail("wavedec[\{i}] \{wname}/\{mname}: coeff[\{k}] mismatch")
      }
    }
    let rec = @wavelets.waverec(coeffs, w, mode~)
    if !arr_close(rec, jnums(jget(case, "rec")), 1.0e-8) {
      fail("waverec[\{i}] \{wname}/\{mname}: reconstruction mismatch")
    }
  }
}

///|
test "upcoef/downcoef match pywt (36 cases)" {
  let cases = jarr(jget(golden_root, "updown"))
  for i, case in cases {
    let wname = jstr(jget(case, "w"))
    let w = @wavelets.Wavelet::new(wname)
    let part = if jstr(jget(case, "part")) == "a" {
      @wavelets.Approx
    } else {
      @wavelets.Detail
    }
    let level = jint(jget(case, "level"))
    let r = if jstr(jget(case, "op")) == "up" {
      @wavelets.upcoef(
        part,
        jnums(jget(case, "c")),
        w,
        level~,
        take=jint(jget(case, "take")),
      )
    } else {
      @wavelets.downcoef(part, jnums(jget(case, "c")), w, level~)
    }
    if !arr_close(r, jnums(jget(case, "r")), 1.0e-9) {
      fail("updown[\{i}] \{wname}: mismatch")
    }
  }
}

///|
test "dwt2/idwt2 match pywt (27 cases)" {
  let cases = jarr(jget(golden_root, "dwt2"))
  for i, case in cases {
    let wname = jstr(jget(case, "w"))
    let mname = jstr(jget(case, "mode"))
    let w = @wavelets.Wavelet::new(wname)
    let mode = @wavelets.Mode::parse(mname)
    let m = @wavelets.Matrix::from_rows(jmat(jget(case, "x")))
    let c = @wavelets.dwt2(m, w, mode~)
    let check = (got : @wavelets.Matrix?, key : String) => {
      if !mat_close(got.unwrap().to_rows(), jmat(jget(case, key)), 1.0e-9) {
        fail("dwt2[\{i}] \{wname}/\{mname}: \{key} mismatch")
      }
    }
    check(c.approx, "a")
    check(c.horizontal, "h")
    check(c.vertical, "v")
    check(c.diagonal, "d")
    let rec = @wavelets.idwt2(c, w, mode~)
    if !mat_close(rec.to_rows(), jmat(jget(case, "rec")), 1.0e-8) {
      fail("idwt2[\{i}] \{wname}/\{mname}: reconstruction mismatch")
    }
  }
}

///|
test "wavedec2/waverec2 match pywt (8 cases)" {
  let cases = jarr(jget(golden_root, "wavedec2"))
  for i, case in cases {
    let wname = jstr(jget(case, "w"))
    let mname = jstr(jget(case, "mode"))
    let w = @wavelets.Wavelet::new(wname)
    let mode = @wavelets.Mode::parse(mname)
    let m = @wavelets.Matrix::from_rows(jmat(jget(case, "x")))
    let dec = @wavelets.wavedec2(
      m,
      w,
      mode~,
      level=Some(jint(jget(case, "level"))),
    )
    if !mat_close(dec.approx.to_rows(), jmat(jget(case, "approx")), 1.0e-9) {
      fail("wavedec2[\{i}] \{wname}: approx mismatch")
    }
    let expected = jarr(jget(case, "details"))
    if dec.details.length() != expected.length() {
      fail("wavedec2[\{i}] \{wname}: detail count mismatch")
    }
    for k in 0..<dec.details.length() {
      let tup = jarr(expected[k])
      let d = dec.details[k]
      if !mat_close(d.horizontal.to_rows(), jmat(tup[0]), 1.0e-9) ||
        !mat_close(d.vertical.to_rows(), jmat(tup[1]), 1.0e-9) ||
        !mat_close(d.diagonal.to_rows(), jmat(tup[2]), 1.0e-9) {
        fail("wavedec2[\{i}] \{wname}: detail[\{k}] mismatch")
      }
    }
    let rec = @wavelets.waverec2(dec, w, mode~)
    if !mat_close(rec.to_rows(), jmat(jget(case, "rec")), 1.0e-8) {
      fail("waverec2[\{i}] \{wname}: reconstruction mismatch")
    }
  }
}

///|
test "cwt matches pywt (20 cases, conv method, precision 12)" {
  let cases = jarr(jget(golden_root, "cwt"))
  for i, case in cases {
    let wname = jstr(jget(case, "w"))
    let w = @wavelets.Wavelet::new(wname)
    let scales = jnums(jget(case, "scales"))
    let (coefs, freqs) = @wavelets.cwt(
      jnums(jget(case, "x")),
      scales,
      w,
      sampling_period=jnum(jget(case, "dt")),
      precision=jint(jget(case, "precision")),
    )
    if !mat_close(coefs.to_rows(), jmat(jget(case, "coefs")), 1.0e-6) {
      fail("cwt[\{i}] \{wname}: coefficients mismatch")
    }
    if !arr_close(freqs, jnums(jget(case, "freqs")), 1.0e-9) {
      fail("cwt[\{i}] \{wname}: frequencies mismatch")
    }
  }
}

///|
test "wavefun matches pywt (32 cases)" {
  let cases = jarr(jget(golden_root, "wavefun"))
  for i, case in cases {
    let wname = jstr(jget(case, "w"))
    let w = @wavelets.Wavelet::new(wname)
    let wf = w.wavefun(level=jint(jget(case, "level")))
    let kind = jstr(jget(case, "kind"))
    if !arr_close(wf.psi, jnums(jget(case, "psi")), 1.0e-9) {
      fail("wavefun[\{i}] \{wname}: psi mismatch")
    }
    if !arr_close(wf.x, jnums(jget(case, "x")), 1.0e-12) {
      fail("wavefun[\{i}] \{wname}: x grid mismatch")
    }
    match kind {
      "orthogonal" =>
        if !arr_close(wf.phi.unwrap(), jnums(jget(case, "phi")), 1.0e-9) {
          fail("wavefun[\{i}] \{wname}: phi mismatch")
        }
      "biorthogonal" =>
        if !arr_close(wf.phi.unwrap(), jnums(jget(case, "phi")), 1.0e-9) ||
          !arr_close(wf.phi_r.unwrap(), jnums(jget(case, "phi_r")), 1.0e-9) ||
          !arr_close(wf.psi_r.unwrap(), jnums(jget(case, "psi_r")), 1.0e-9) {
          fail("wavefun[\{i}] \{wname}: biorthogonal sides mismatch")
        }
      _ => ()
    }
  }
}

///|
test "central_frequency matches pywt (20 cases)" {
  let cases = jarr(jget(golden_root, "centrfrq"))
  for i, case in cases {
    let wname = jstr(jget(case, "w"))
    let w = @wavelets.Wavelet::new(wname)
    let f = @wavelets.central_frequency(
      w,
      precision=jint(jget(case, "precision")),
    )
    if tabs(f - jnum(jget(case, "f"))) > 1.0e-12 {
      fail("centrfrq[\{i}] \{wname}: \{f} != \{jnum(jget(case, "f"))}")
    }
  }
}

///|
test "integrate_wavelet matches pywt (5 cases)" {
  let cases = jarr(jget(golden_root, "intwave"))
  for i, case in cases {
    let wname = jstr(jget(case, "w"))
    let w = @wavelets.Wavelet::new(wname)
    let (int_psi, x) = @wavelets.integrate_wavelet(
      w,
      precision=jint(jget(case, "precision")),
    )
    if !arr_close(int_psi, jnums(jget(case, "int_psi")), 1.0e-9) ||
      !arr_close(x, jnums(jget(case, "x")), 1.0e-12) {
      fail("intwave[\{i}] \{wname}: mismatch")
    }
  }
}

///|
test "pad matches pywt (81 cases, all 9 modes)" {
  let cases = jarr(jget(golden_root, "pad"))
  for i, case in cases {
    let mname = jstr(jget(case, "mode"))
    let mode = @wavelets.Mode::parse(mname)
    let r = @wavelets.pad(
      jnums(jget(case, "x")),
      jint(jget(case, "left")),
      jint(jget(case, "right")),
      mode,
    )
    if !arr_close(r, jnums(jget(case, "r")), 1.0e-12) {
      fail("pad[\{i}] \{mname}: mismatch")
    }
  }
}

///|
test "threshold matches pywt (10 cases, boundary semantics)" {
  let cases = jarr(jget(golden_root, "threshold"))
  for i, case in cases {
    let mname = jstr(jget(case, "mode"))
    let mode = match mname {
      "soft" => @wavelets.Soft
      "hard" => @wavelets.Hard
      "greater" => @wavelets.Greater
      "less" => @wavelets.Less
      "garrote" => @wavelets.Garrote
      _ => panic()
    }
    let r = @wavelets.threshold(
      jnums(jget(case, "x")),
      jnum(jget(case, "t")),
      mode~,
    )
    if !arr_close(r, jnums(jget(case, "r")), 1.0e-12) {
      fail("threshold[\{i}] \{mname}: mismatch")
    }
  }
}

///|
test "dwt_max_level matches pywt (40 cases)" {
  let cases = jarr(jget(golden_root, "maxlevel"))
  for i, case in cases {
    let r = @wavelets.dwt_max_level(
      jint(jget(case, "n")),
      jint(jget(case, "f")),
    )
    if r != jint(jget(case, "r")) {
      fail("maxlevel[\{i}]: \{r} != \{jint(jget(case, "r"))}")
    }
  }
}

///|
test "dwt_coeff_len matches pywt (64 cases)" {
  let cases = jarr(jget(golden_root, "coefflen"))
  for i, case in cases {
    let mode = @wavelets.Mode::parse(jstr(jget(case, "mode")))
    let r = @wavelets.dwt_coeff_len(
      jint(jget(case, "n")),
      jint(jget(case, "f")),
      mode,
    )
    if r != jint(jget(case, "r")) {
      fail("coefflen[\{i}]: \{r} != \{jint(jget(case, "r"))}")
    }
  }
}

///|
test "qmf matches pywt" {
  let cases = jarr(jget(golden_root, "qmf"))
  for i, case in cases {
    let r = @wavelets.qmf(jnums(jget(case, "x")))
    if !arr_close(r, jnums(jget(case, "r")), 1.0e-15) {
      fail("qmf[\{i}]: mismatch")
    }
  }
}

///|
test "orthogonal_filter_bank matches pywt" {
  let cases = jarr(jget(golden_root, "orthfilt"))
  for i, case in cases {
    let (dl, dh, rl, rh) = @wavelets.orthogonal_filter_bank(
      jnums(jget(case, "x")),
    )
    if !arr_close(dl, jnums(jget(case, "dec_lo")), 1.0e-12) ||
      !arr_close(dh, jnums(jget(case, "dec_hi")), 1.0e-12) ||
      !arr_close(rl, jnums(jget(case, "rec_lo")), 1.0e-12) ||
      !arr_close(rh, jnums(jget(case, "rec_hi")), 1.0e-12) {
      fail("orthfilt[\{i}]: mismatch")
    }
  }
}
