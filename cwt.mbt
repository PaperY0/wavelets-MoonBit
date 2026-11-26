// 连续小波变换与频率换算。
//
// 移植来源：PyWavelets `pywt/_cwt.py::cwt`（conv 方法）与
// `pywt/_functions.py`（integrate_wavelet / central_frequency /
// scale2frequency / frequency2scale）（MIT）。
// 积分网格精度 precision 默认 12（2^12 点），与上游一致。

///|
/// 矩形法累积积分 psi：int_psi[k] = step * sum_{i<=k} psi[i]。
/// 返回 (int_psi, x)。支持连续小波与正交离散小波；双正交离散小波
/// 上游返回两个积分，v0.1 暂不支持（见 known-deviations）。
pub fn integrate_wavelet(
  wavelet : Wavelet,
  precision~ : Int32 = 8,
) -> (Array[Double], Array[Double]) raise WaveletError {
  if wavelet.kind is Discrete && wavelet.biorthogonal && !wavelet.orthogonal {
    raise NotSupported(
      "integrate_wavelet for biorthogonal discrete wavelets is not in v0.1",
    )
  }
  let wf = wavelet.wavefun(level=precision)
  let psi = wf.psi
  let x = wf.x
  let step = x[1] - x[0]
  let mut acc = 0.0
  let int_psi = Array::makei(psi.length(), i => {
    acc += psi[i]
    acc * step
  })
  (int_psi, x)
}

///|
/// psi 的中心频率（归一化到采样频率）。v0.1 仅支持连续小波
/// （离散小波需要任意长度 FFT，见 known-deviations）。
pub fn central_frequency(
  wavelet : Wavelet,
  precision~ : Int32 = 8,
) -> Double raise WaveletError {
  guard wavelet.kind is Continuous else {
    raise NotSupported("central_frequency for discrete wavelets is not in v0.1")
  }
  let wf = wavelet.wavefun(level=precision)
  let psi = wf.psi
  let x = wf.x
  let domain = x[x.length() - 1] - x[0]
  let (re, im) = fft_real(psi)
  let n = psi.length()
  // pywt: index = argmax(|fft[1:]|) + 2；切片 argmax = best-1，故 index = best+1。
  let mut best = 0
  let mut best_mag = -1.0
  for k in 1..<n {
    let mag = re[k] * re[k] + im[k] * im[k]
    if mag > best_mag {
      best_mag = mag
      best = k
    }
  }
  let mut index = best + 1
  if index > n / 2 {
    index = n - index + 2
  }
  (index - 1).to_double() / domain
}

///|
/// 尺度 -> 归一化频率：central_frequency / scale。
pub fn scale2frequency(
  wavelet : Wavelet,
  scale : Double,
  precision~ : Int32 = 8,
) -> Double raise WaveletError {
  if scale <= 0.0 {
    raise InvalidArgument("scale must be positive")
  }
  central_frequency(wavelet, precision~) / scale
}

///|
/// 归一化频率 -> 尺度：central_frequency / freq。
pub fn frequency2scale(
  wavelet : Wavelet,
  freq : Double,
  precision~ : Int32 = 8,
) -> Double raise WaveletError {
  if freq <= 0.0 {
    raise InvalidArgument("frequency must be positive")
  }
  central_frequency(wavelet, precision~) / freq
}

///|
/// 一维连续小波变换（conv 方法）。
///
/// 返回 (coefs, frequencies)：coefs 为 scales.length() × data.length()
/// 矩阵；frequencies = central_frequency / (scale * sampling_period)。
/// 系数本身与 sampling_period 无关（与上游一致）。
pub fn cwt(
  data : Array[Double],
  scales : Array[Double],
  wavelet : Wavelet,
  sampling_period~ : Double = 1.0,
  precision~ : Int32 = 12,
) -> (Matrix, Array[Double]) raise WaveletError {
  guard wavelet.kind is Continuous else { raise NotContinuous(wavelet.name) }
  if data.length() == 0 {
    raise InvalidArgument("cwt requires a non-empty signal")
  }
  if scales.length() == 0 {
    raise InvalidArgument("cwt requires at least one scale")
  }
  for s in scales {
    if s <= 0.0 {
      raise InvalidArgument("scales must be positive")
    }
  }
  if sampling_period <= 0.0 {
    raise InvalidArgument("sampling_period must be positive")
  }
  let (int_psi, x) = integrate_wavelet(wavelet, precision~)
  let step = x[1] - x[0]
  let span = x[x.length() - 1] - x[0]
  let n = data.length()
  let rows : Array[Array[Double]] = []
  for scale in scales {
    // j = floor(arange(scale*span + 1) / (scale*step))，截取 < len(int_psi)。
    let count = @math.ceil(scale * span + 1.0).to_int()
    let denom = scale * step
    let selected : Array[Double] = []
    for k in 0..<count {
      let j = (k.to_double() / denom).to_int()
      if j < int_psi.length() {
        selected.push(int_psi[j])
      }
    }
    // 反转。
    let m = selected.length()
    let kernel = Array::makei(m, i => selected[m - 1 - i])
    let conv = conv_full(data, kernel)
    // coef = -sqrt(scale) * diff(conv)。
    let coef_len = conv.length() - 1
    let coef = Array::makei(coef_len, i => {
      -scale.sqrt() * (conv[i + 1] - conv[i])
    })
    let d = (coef_len - n).to_double() / 2.0
    let row = if d > 0.0 {
      let start = @math.floor(d).to_int()
      let end = coef_len - @math.ceil(d).to_int()
      Array::makei(end - start, i => coef[start + i])
    } else if d < 0.0 {
      raise InvalidArgument("selected scale \{scale} is too small")
    } else {
      coef
    }
    rows.push(row)
  }
  let cf = central_frequency(wavelet, precision~)
  let freqs = Array::makei(scales.length(), i => {
    cf / scales[i] / sampling_period
  })
  (Matrix::from_rows(rows), freqs)
}
