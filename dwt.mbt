// 单层离散小波变换与直接重构。
//
// 移植来源：PyWavelets `pywt/_dwt.py` + C 端 convolution.template.c（MIT）。
// 关键约定（已用 tools/probe.py 对上游逐一验证）：
// - 输出长度：非 periodization 为 floor((N + F - 1) / 2)，periodization 为 ceil(N/2)；
// - 左延拓 F-2 个样本，右延拓由输出长度反推；
// - 分析卷积为滤波器反转 + stride 2；
// - idwt 输出长度 2L - F + 2（periodization 为 2L，循环卷积）。

///|
/// 系数分支选择（upcoef/downcoef 用）。
#derive(Eq, Debug)
pub(all) enum CoeffPart {
  /// 近似（低通）分支
  Approx
  /// 细节（高通）分支
  Detail
}

///|
/// 单层离散小波变换，返回 (cA, cD)。
///
/// 长度约定与 pywt.dwt 一致：奇数长度输入在非 periodization 模式下
/// 可能产生“多一个”的系数（上游行为，忠实保留）。
pub fn dwt(
  data : Array[Double],
  wavelet : Wavelet,
  mode~ : Mode = Symmetric,
) -> (Array[Double], Array[Double]) raise WaveletError {
  guard wavelet.kind is Discrete else { raise NotDiscrete(wavelet.name) }
  let n = data.length()
  if n == 0 {
    raise InvalidArgument("dwt requires a non-empty signal")
  }
  // 上游限制：reflect/antireflect 的周期为 2N-2，要求 N>=2（pad/idwt 无此限制）。
  if n < 2 && (mode is Reflect || mode is AntiReflect) {
    raise InvalidArgument(
      "reflect/antireflect modes require at least 2 samples",
    )
  }
  let f = wavelet.dec_len()
  // periodization：奇数长度先复制末值补偶，再周期延拓。
  let src = if mode is Periodization && n % 2 == 1 {
    concat_arrays(data, [data[n - 1]])
  } else {
    data
  }
  let nn = src.length()
  let out_len = dwt_coeff_len(n, f, mode)
  // periodization：左右各 (F-2)/2 的对称周期延拓（所有离散小波 F 均为偶数）；
  // 其余模式：左 F-2、右补齐到 2*out_len。
  let (left, right) = if mode is Periodization {
    ((f - 2) / 2, (f - 2) / 2)
  } else {
    (f - 2, 2 * out_len - nn)
  }
  let eff = if mode is Periodization { Periodic } else { mode }
  let ext = pad(src, left, right, eff)
  (
    conv_down(ext, wavelet.dec_lo, out_len),
    conv_down(ext, wavelet.dec_hi, out_len),
  )
}

///|
/// 单层逆离散小波变换。`cA`/`cD` 之一可为 None（等价于补零分支），
/// 但不能同时为 None；两者长度必须一致（与 pywt.idwt 一致，
/// “cA 比 cD 长 1”的裁剪由 waverec 负责）。
pub fn idwt(
  cA : Array[Double]?,
  cD : Array[Double]?,
  wavelet : Wavelet,
  mode~ : Mode = Symmetric,
) -> Array[Double] raise WaveletError {
  guard wavelet.kind is Discrete else { raise NotDiscrete(wavelet.name) }
  let (a, d) = match (cA, cD) {
    (Some(a), Some(d)) => (a, d)
    (Some(a), None) => (a, Array::make(a.length(), 0.0))
    (None, Some(d)) => (Array::make(d.length(), 0.0), d)
    (None, None) =>
      raise InvalidArgument("at least one of cA/cD must be provided")
  }
  let l = a.length()
  if l != d.length() {
    raise ShapeMismatch("cA and cD must have the same length")
  }
  if l == 0 {
    raise InvalidArgument("idwt requires non-empty coefficients")
  }
  let f = wavelet.rec_len()
  let per = mode is Periodization
  let out_len = if per { 2 * l } else { 2 * l - f + 2 }
  if out_len <= 0 {
    raise InvalidArgument("coefficients too short for this wavelet")
  }
  let up_a = upsample2(a)
  let up_d = upsample2(d)
  let period = 2 * l
  // 周期模式窗口偏移 (F-2)/2，其余模式 F-2（与 pywt 相位一致）。
  let off = if per { (f - 2) / 2 } else { f - 2 }
  Array::makei(out_len, m => {
    let n = m + off
    let mut acc = 0.0
    for k in 0..<f {
      let idx = n - k
      let j = if per { imod(idx, period) } else { idx }
      if j >= 0 && j < period && j % 2 == 0 {
        acc += wavelet.rec_lo[k] * up_a[j] + wavelet.rec_hi[k] * up_d[j]
      }
    }
    acc
  })
}

///|
/// 由单个分支系数直接重构（不做边界裁剪），对应 pywt.upcoef。
/// `take > 0` 时取结果的中央 `take` 个样本。
pub fn upcoef(
  part : CoeffPart,
  coeffs : Array[Double],
  wavelet : Wavelet,
  level~ : Int32 = 1,
  take~ : Int32 = 0,
) -> Array[Double] raise WaveletError {
  guard wavelet.kind is Discrete else { raise NotDiscrete(wavelet.name) }
  if level < 1 {
    raise InvalidLevel("upcoef level must be >= 1")
  }
  if coeffs.length() == 0 {
    raise InvalidArgument("upcoef requires non-empty coefficients")
  }
  let mut a = coeffs
  for i in 0..<level {
    // 细节分支只在第一层走 rec_hi，其余层走 rec_lo（多层重构语义）；
    // 近似分支每层都走 rec_lo。与 pywt upcoef 一致。
    let filt = if i == 0 {
      match part {
        Approx => wavelet.rec_lo
        Detail => wavelet.rec_hi
      }
    } else {
      wavelet.rec_lo
    }
    a = upcoef_step(upsample2(a), filt)
  }
  if take == 0 {
    a
  } else {
    if take > a.length() {
      raise InvalidArgument("take exceeds reconstructed length")
    }
    let start = (a.length() - take) / 2
    Array::makei(take, i => a[start + i])
  }
}

///|
/// 只计算某一层的一个分支系数，对应 pywt.downcoef。
pub fn downcoef(
  part : CoeffPart,
  data : Array[Double],
  wavelet : Wavelet,
  mode~ : Mode = Symmetric,
  level~ : Int32 = 1,
) -> Array[Double] raise WaveletError {
  if level < 1 {
    raise InvalidLevel("downcoef level must be >= 1")
  }
  let mut a = data
  for i in 0..<level {
    let (ca, cd) = dwt(a, wavelet, mode~)
    // 中间层一律沿近似支下采样，仅在最后一层选取 part（与 pywt._downcoef 一致）。
    a = if i == level - 1 {
      match part {
        Approx => ca
        Detail => cd
      }
    } else {
      ca
    }
  }
  a
}

///|
/// 最大有用分解层数 floor(log2(N / (F - 1)))，下限 0。
/// 与 pywt.dwt_max_level 一致。
pub fn dwt_max_level(
  data_len : Int32,
  filter_len : Int32,
) -> Int32 raise WaveletError {
  if filter_len < 2 {
    raise InvalidArgument("filter length must be >= 2")
  }
  if data_len <= 0 {
    return 0
  }
  let ratio = data_len.to_double() / (filter_len - 1).to_double()
  let lvl = @math.floor(@math.log2(ratio)).to_int()
  if lvl < 0 {
    0
  } else {
    lvl
  }
}

///|
/// DWT 输出系数长度。非 periodization：floor((N + F - 1) / 2)；
/// periodization：ceil(N / 2)。与 pywt.dwt_coeff_len 一致。
pub fn dwt_coeff_len(
  data_len : Int32,
  filter_len : Int32,
  mode : Mode,
) -> Int32 {
  if mode is Periodization {
    (data_len + 1) / 2
  } else {
    (data_len + filter_len - 1) / 2
  }
}
