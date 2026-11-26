// 注册表与对象模型测试。

///|
test "wavelist counts match the supported set" {
  assert_eq(@wavelets.wavelist(kind=Some(@wavelets.Discrete)).length(), 106)
  assert_eq(@wavelets.wavelist(kind=Some(@wavelets.Continuous)).length(), 10)
  assert_eq(@wavelets.wavelist().length(), 116)
}

///|
test "families list" {
  assert_eq(@wavelets.families(), [
    "haar", "db", "sym", "coif", "bior", "rbio", "dmey", "gaus", "mexh", "morl",
  ])
}

///|
test "unknown or miscased names are rejected" {
  assert_true(try_new("nope") is None)
  assert_true(try_new("DB2") is None)
  assert_true(try_new("") is None)
  assert_true(try_new("db2") is Some(_))
}

///|
test "db1 is the Haar filter bank under the Daubechies name" {
  let db1 = @wavelets.Wavelet::new("db1")
  let haar = @wavelets.Wavelet::new("haar")
  assert_eq(db1.dec_lo, haar.dec_lo)
  assert_eq(db1.dec_hi, haar.dec_hi)
  assert_eq(db1.family_name, "Daubechies")
  assert_eq(haar.family_name, "Haar")
}

///|
test "filter lengths of representative wavelets" {
  let cases = [
    ("haar", 2),
    ("db2", 4),
    ("db38", 76),
    ("sym20", 40),
    ("coif1", 6),
    ("coif17", 102),
    ("dmey", 62),
    ("bior3.5", 12),
  ]
  for (name, len) in cases {
    assert_eq(@wavelets.Wavelet::new(name).dec_len(), len)
  }
}

///|
test "orthogonality and symmetry metadata" {
  let db2 = @wavelets.Wavelet::new("db2")
  assert_true(db2.orthogonal)
  assert_true(db2.biorthogonal)
  assert_eq(db2.symmetry, @wavelets.Asymmetric)
  assert_eq(db2.vm_psi, 2)
  let sym4 = @wavelets.Wavelet::new("sym4")
  assert_eq(sym4.symmetry, @wavelets.NearSymmetric)
  let bior = @wavelets.Wavelet::new("bior3.5")
  assert_true(!bior.orthogonal)
  assert_true(bior.biorthogonal)
  assert_eq(bior.symmetry, @wavelets.Symmetric)
  let haar = @wavelets.Wavelet::new("haar")
  assert_eq(haar.vm_psi, 1)
}

///|
test "continuous wavelets carry support bounds and empty filter banks" {
  let mexh = @wavelets.Wavelet::new("mexh")
  assert_eq(mexh.kind, @wavelets.Continuous)
  assert_eq(mexh.lower_bound, -8.0)
  assert_eq(mexh.upper_bound, 8.0)
  assert_eq(mexh.dec_len(), 0)
  let gaus1 = @wavelets.Wavelet::new("gaus1")
  assert_eq(gaus1.lower_bound, -5.0)
  assert_eq(gaus1.upper_bound, 5.0)
}

///|
test "mode parsing accepts pywt names and legacy aliases" {
  assert_eq(@wavelets.Mode::parse("per"), @wavelets.Periodization)
  assert_eq(@wavelets.Mode::parse("sym"), @wavelets.Symmetric)
  assert_eq(@wavelets.Mode::parse("zpd"), @wavelets.Zero)
  assert_eq(@wavelets.Mode::parse("ppd"), @wavelets.Periodic)
  assert_eq(@wavelets.Mode::parse("sp1"), @wavelets.Smooth)
  assert_eq(@wavelets.Mode::parse("cpd"), @wavelets.Constant)
  for mode in [@wavelets.Zero, @wavelets.Symmetric, @wavelets.Periodization] {
    assert_eq(@wavelets.Mode::parse(mode.name()), mode)
  }
}

///|
fn try_parse_mode(s : String) -> @wavelets.Mode? {
  Some(@wavelets.Mode::parse(s)) catch {
    _ => None
  }
}

///|
test "unknown mode is rejected" {
  assert_true(try_parse_mode("bogus") is None)
}

///|
fn dwt_rejects_continuous() -> Bool raise @wavelets.WaveletError {
  let w = @wavelets.Wavelet::new("mexh")
  let r = Some(@wavelets.dwt([1.0, 2.0], w)) catch { _ => None }
  r is None
}

///|
test "discrete-only APIs reject continuous wavelets" {
  assert_true(dwt_rejects_continuous())
}

///|
fn dwt_raises(
  data : Array[Double],
  mode : @wavelets.Mode,
) -> Bool raise @wavelets.WaveletError {
  let w = @wavelets.Wavelet::new("haar")
  let r = Some(@wavelets.dwt(data, w, mode~)) catch { _ => None }
  r is None
}

///|
fn pad_n1_reflect() -> Array[Double] raise @wavelets.WaveletError {
  @wavelets.pad([5.0], 2, 2, @wavelets.Mode::Reflect)
}

///|
test "dwt rejects n<2 for reflect/antireflect (upstream restriction)" {
  // 上游 dwt 对 [anti]reflect 要求 N>=2；pad 无此限制（n=1 返回常数延拓）。
  assert_true(dwt_raises([5.0], @wavelets.Mode::Reflect))
  assert_true(dwt_raises([5.0], @wavelets.Mode::AntiReflect))
  assert_true(!dwt_raises([5.0], @wavelets.Mode::AntiSymmetric))
  assert_true(!dwt_raises([5.0, 6.0], @wavelets.Mode::Reflect))
  assert_eq(pad_n1_reflect(), [5.0, 5.0, 5.0, 5.0, 5.0])
}
