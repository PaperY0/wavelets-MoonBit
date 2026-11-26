// 阈值去噪与滤波器组工具。
//
// 移植来源：PyWavelets `pywt/_thresholding.py` 与
// `pywt/_functions.py::qmf/orthogonal_filter_bank`（MIT）。
// 边界约定（已对上游验证）：hard 含等号（|x| >= t 保留），
// greater 严格（|x| > t），less 含等号（|x| <= t），
// garrote 在 |x| >= t 时取 x*(1 - t^2/x^2)（t=0 且 x=0 时为 NaN，上游行为）。

///|
/// 阈值模式。
#derive(Eq, Debug)
pub(all) enum ThresholdMode {
  Soft
  Hard
  Greater
  Less
  Garrote
}

///|
/// 阈值处理。`value` 必须非负；`substitute` 为被置零分支的替代值
/// （soft/garrote 固定用 0，与上游一致）。
pub fn threshold(
  data : Array[Double],
  value : Double,
  mode~ : ThresholdMode = Soft,
  substitute~ : Double = 0.0,
) -> Array[Double] raise WaveletError {
  if value < 0.0 {
    raise InvalidArgument("threshold value must be >= 0")
  }
  Array::makei(data.length(), i => {
    let x = data[i]
    let mag = dabs(x)
    match mode {
      // soft/garrote 复刻 pywt 的 x*clip(1-t/|x|) 形式：
      // x=0 且 t=0 时 0/0 -> NaN；x<0 且截断为 0 时保留 -0.0。
      Soft =>
        if x == 0.0 {
          if value == 0.0 {
            0.0 / 0.0
          } else {
            0.0
          }
        } else {
          let f = 1.0 - value / mag
          x * (if f < 0.0 { 0.0 } else { f })
        }
      Hard => if mag >= value { x } else { substitute }
      Greater => if x >= value { x } else { substitute }
      Less => if x <= value { x } else { substitute }
      Garrote =>
        if x == 0.0 {
          if value == 0.0 {
            0.0 / 0.0
          } else {
            0.0
          }
        } else {
          let f = 1.0 - value * value / (x * x)
          x * (if f < 0.0 { 0.0 } else { f })
        }
    }
  })
}

///|
/// 正交镜像滤波器：反转后对奇数下标取负（对应 pywt.qmf）。
pub fn qmf(filt : Array[Double]) -> Array[Double] {
  let n = filt.length()
  Array::makei(n, i => {
    let v = filt[n - 1 - i]
    if i % 2 == 1 {
      -v
    } else {
      v
    }
  })
}

///|
/// 由尺度滤波器构造正交滤波器组 (dec_lo, dec_hi, rec_lo, rec_hi)。
/// 输入长度必须为偶数（对应 pywt.orthogonal_filter_bank）。
pub fn orthogonal_filter_bank(
  scaling_filter : Array[Double],
) -> (Array[Double], Array[Double], Array[Double], Array[Double]) raise WaveletError {
  let n = scaling_filter.length()
  if n == 0 || n % 2 != 0 {
    raise InvalidArgument("scaling filter length must be even and non-empty")
  }
  let mut sum = 0.0
  for v in scaling_filter {
    sum += v
  }
  if sum == 0.0 {
    raise InvalidArgument("scaling filter must have non-zero sum")
  }
  let scale = Double::sqrt(2.0) / sum
  let rec_lo = Array::makei(n, i => scaling_filter[i] * scale)
  let dec_lo = Array::makei(n, i => rec_lo[n - 1 - i])
  let rec_hi = qmf(rec_lo)
  let dec_hi = Array::makei(n, i => rec_hi[n - 1 - i])
  (dec_lo, dec_hi, rec_lo, rec_hi)
}
