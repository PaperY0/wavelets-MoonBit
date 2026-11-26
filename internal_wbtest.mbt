// 白盒测试：内部数值工具（同包可见，无前缀）。

///|
test "upsample2 inserts trailing zero" {
  assert_eq(upsample2([1.0, 2.0, 3.0]), [1.0, 0.0, 2.0, 0.0, 3.0, 0.0])
  assert_eq(upsample2([7.0]), [7.0, 0.0])
}

///|
test "conv_full matches naive convolution" {
  assert_eq(conv_full([1.0, 2.0], [3.0, 4.0]), [3.0, 10.0, 8.0])
  let a = [1.0, -2.0, 3.0, 0.5]
  let b = [0.25, 0.75, -1.0]
  let got = conv_full(a, b)
  let want = Array::makei(a.length() + b.length() - 1, n => {
    let mut acc = 0.0
    for k in 0..<a.length() {
      let j = n - k
      if j >= 0 && j < b.length() {
        acc += a[k] * b[j]
      }
    }
    acc
  })
  assert_eq(got, want)
}

///|
test "conv_down uses reversed filter and stride 2" {
  // 均值滤波：反转不变，stride 2 取相邻对。
  assert_eq(conv_down([1.0, 2.0, 3.0, 4.0], [0.5, 0.5], 2), [1.5, 3.5])
  // 非对称滤波器验证反转相位。
  assert_eq(conv_down([1.0, 2.0, 3.0, 4.0], [1.0, 0.0], 2), [2.0, 4.0])
}

///|
test "imod and ifloor_div handle negatives" {
  assert_eq(imod(-1, 5), 4)
  assert_eq(imod(-6, 5), 4)
  assert_eq(imod(7, 5), 2)
  assert_eq(ifloor_div(-3, 4), -1)
  assert_eq(ifloor_div(-4, 4), -1)
  assert_eq(ifloor_div(-5, 4), -2)
  assert_eq(ifloor_div(5, 4), 1)
}

///|
test "extend_sample spot checks across modes" {
  let x = [1.0, 2.0, 3.0]
  assert_eq(extend_sample(x, -1, Symmetric), 1.0)
  assert_eq(extend_sample(x, -1, Reflect), 2.0)
  assert_eq(extend_sample(x, -4, AntiSymmetric), 3.0)
  assert_eq(extend_sample(x, -1, AntiSymmetric), -1.0)
  assert_eq(extend_sample(x, -4, AntiReflect), -3.0)
  assert_eq(extend_sample(x, -1, AntiReflect), 0.0)
  assert_eq(extend_sample(x, 6, AntiReflect), 7.0)
  assert_eq(extend_sample(x, -4, Smooth), -3.0)
  assert_eq(extend_sample(x, 4, Smooth), 5.0)
  assert_eq(extend_sample(x, -1, Periodic), 3.0)
  assert_eq(extend_sample(x, 3, Constant), 3.0)
  assert_eq(extend_sample(x, -1, Zero), 0.0)
  assert_eq(extend_sample([5.0], -3, Reflect), 5.0)
  assert_eq(extend_sample([5.0], -3, Smooth), 5.0)
}

///|
test "fft matches naive DFT" {
  for n in [4, 8, 16] {
    let x = Array::makei(n, i => (i * 17 % 13).to_double() - 6.0)
    let (re, im) = fft_real(x)
    for k in 0..<n {
      let mut wr = 0.0
      let mut wi = 0.0
      for t in 0..<n {
        let ang = -2.0 * 3.141592653589793 * (k * t).to_double() / n.to_double()
        wr += x[t] * @math.cos(ang)
        wi += x[t] * @math.sin(ang)
      }
      assert_true(dabs(re[k] - wr) < 1.0e-9)
      assert_true(dabs(im[k] - wi) < 1.0e-9)
    }
  }
}

///|
test "hermite polynomials" {
  assert_eq(hermite(0, 1.5), 1.0)
  assert_eq(hermite(1, 1.5), 3.0)
  assert_eq(hermite(2, 1.5), 7.0) // 4t^2 - 2
  assert_eq(hermite(3, 1.0), -4.0) // 8t^3 - 12t
}

///|
test "linspace endpoints and spacing" {
  assert_eq(linspace(0.0, 1.0, 5), [0.0, 0.25, 0.5, 0.75, 1.0])
  assert_eq(linspace(-8.0, 8.0, 3), [-8.0, 0.0, 8.0])
  assert_eq(linspace(2.0, 9.0, 1), [2.0])
}

///|
test "upcoef_step drops the last full-convolution sample" {
  // up = [sqrt2, 0], filt = db2 rec_lo：输出前 4 项即 phi 级联首层。
  let s = Double::sqrt(2.0)
  let w = Wavelet::new("db2")
  let got = upcoef_step([s, 0.0], w.rec_lo)
  assert_eq(got.length(), 4)
  assert_true(dabs(got[0] - 0.6830127018922193) < 1.0e-12)
  assert_true(dabs(got[1] - 1.1830127018922194) < 1.0e-12)
}
