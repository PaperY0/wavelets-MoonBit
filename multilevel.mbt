// 多层一维分解与重构。
//
// 移植来源：PyWavelets `pywt/_multilevel.py::wavedec/waverec`（MIT）。
// waverec 的裁剪规则与上游一致：每步 idwt 之前，若近似系数比细节系数
// 长 1，则截去末尾一个样本；最终输出长度由最细层细节长度决定
// （奇数长度输入可能比原信号长 1，属上游行为）。

///|
/// 多层分解，返回 [cA_n, cD_n, cD_{n-1}, ..., cD_1]。
///
/// `level` 省略时取 dwt_max_level；level 为 0 时返回 [data]。
/// 与上游的差异：level 超过最大有用层数时本库报错而非警告
/// （见 docs/known-deviations.md）。
pub fn wavedec(
  data : Array[Double],
  wavelet : Wavelet,
  mode~ : Mode = Symmetric,
  level~ : Int32? = None,
) -> Array[Array[Double]] raise WaveletError {
  guard wavelet.kind is Discrete else { raise NotDiscrete(wavelet.name) }
  if data.length() == 0 {
    raise InvalidArgument("wavedec requires a non-empty signal")
  }
  let max_lvl = dwt_max_level(data.length(), wavelet.dec_len())
  let lvl = match level {
    None => max_lvl
    Some(l) => l
  }
  if lvl < 0 {
    raise InvalidLevel("level must be >= 0")
  }
  if lvl > max_lvl {
    raise InvalidLevel(
      "level \{lvl} exceeds maximum useful level \{max_lvl} for this signal/wavelet",
    )
  }
  if lvl == 0 {
    return [data]
  }
  let mut a = data
  let details : Array[Array[Double]] = []
  for _ in 0..<lvl {
    let (ca, cd) = dwt(a, wavelet, mode~)
    a = ca
    details.push(cd)
  }
  let out : Array[Array[Double]] = [a]
  let n = details.length()
  for i in 0..<n {
    out.push(details[n - 1 - i])
  }
  out
}

///|
/// 多层重构，输入为 wavedec 的系数列表 [cA_n, cD_n, ..., cD_1]。
pub fn waverec(
  coeffs : Array[Array[Double]],
  wavelet : Wavelet,
  mode~ : Mode = Symmetric,
) -> Array[Double] raise WaveletError {
  guard wavelet.kind is Discrete else { raise NotDiscrete(wavelet.name) }
  if coeffs.length() == 0 {
    raise InvalidArgument("coefficient list must contain at least one array")
  }
  if coeffs.length() == 1 {
    return coeffs[0]
  }
  let mut a = coeffs[0]
  for i in 1..<coeffs.length() {
    let d = coeffs[i]
    if a.length() == d.length() + 1 {
      a = Array::makei(d.length(), j => a[j])
    } else if a.length() != d.length() {
      raise ShapeMismatch("approximation/detail length mismatch at level \{i}")
    }
    a = idwt(Some(a), Some(d), wavelet, mode~)
  }
  a
}

///|
/// wmaxlev：dwt_max_level 的小波对象便捷版（对应 pywt.wmaxlev）。
pub fn wmaxlev(data_len : Int32, wavelet : Wavelet) -> Int32 raise WaveletError {
  guard wavelet.kind is Discrete else { raise NotDiscrete(wavelet.name) }
  dwt_max_level(data_len, wavelet.dec_len())
}
