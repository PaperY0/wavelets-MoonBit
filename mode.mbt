// 信号扩展模式与 pad 实现。
//
// 移植来源：PyWavelets `pywt/_dwt.py::pad` 与 C 端 convolution.template.c 的
// 边界语义（MIT）。九种模式与上游逐一对齐，包括 periodization 的
// “奇数长度先复制末值补偶、再周期延拓”约定。

///|
/// 信号扩展模式（对应 pywt.Modes）。
#derive(Eq, Debug)
pub(all) enum Mode {
  /// 零延拓
  Zero
  /// 边缘值延拓（等价 numpy 的 edge，注意不是补零）
  Constant
  /// 整样本对称（含边缘）
  Symmetric
  /// 周期延拓
  Periodic
  /// 一阶导数线性外插
  Smooth
  /// 周期化：输出长度 ceil(N/2)，奇数长度先复制末值补偶
  Periodization
  /// 半样本对称（不含边缘）
  Reflect
  /// 整样本反对称（分段变号）
  AntiSymmetric
  /// 半样本反对称（奇反射）
  AntiReflect
}

///|
/// 解析模式名；接受 pywt 的全名与旧别名（zpd/cpd/sym/ppd/sp1/per）。
pub fn Mode::parse(s : String) -> Mode raise WaveletError {
  match s {
    "zero" | "zpd" => Zero
    "constant" | "cpd" => Constant
    "symmetric" | "sym" => Symmetric
    "periodic" | "ppd" => Periodic
    "smooth" | "sp1" => Smooth
    "periodization" | "per" => Periodization
    "reflect" => Reflect
    "antisymmetric" => AntiSymmetric
    "antireflect" => AntiReflect
    _ => raise InvalidMode(s)
  }
}

///|
/// 模式的 pywt 规范名。
pub fn Mode::name(self : Mode) -> String {
  match self {
    Zero => "zero"
    Constant => "constant"
    Symmetric => "symmetric"
    Periodic => "periodic"
    Smooth => "smooth"
    Periodization => "periodization"
    Reflect => "reflect"
    AntiSymmetric => "antisymmetric"
    AntiReflect => "antireflect"
  }
}

///|
/// 按给定模式延拓一维信号：左侧补 `left` 个、右侧补 `right` 个样本。
/// 与 `pywt.pad(x, (left, right), mode)` 行为一致。
pub fn pad(
  x : Array[Double],
  left : Int32,
  right : Int32,
  mode : Mode,
) -> Array[Double] raise WaveletError {
  if left < 0 || right < 0 {
    raise InvalidArgument("pad widths must be >= 0")
  }
  let n = x.length()
  if n == 0 {
    raise InvalidArgument("cannot pad an empty signal")
  }
  // periodization：奇数长度先复制末值补偶，再做周期延拓。
  let src = if mode is Periodization && n % 2 == 1 {
    concat_arrays(x, [x[n - 1]])
  } else {
    x
  }
  let nn = src.length()
  let eff = if mode is Periodization { Periodic } else { mode }
  Array::makei(left + nn + right, i => extend_sample(src, i - left, eff))
}

///|
/// 计算延拓信号在相对位置 p（可为负或 >= n）处的样本值。
fn extend_sample(src : Array[Double], p : Int32, mode : Mode) -> Double {
  let n = src.length()
  if p >= 0 && p < n {
    return src[p]
  }
  match mode {
    Zero => 0.0
    Constant => if p < 0 { src[0] } else { src[n - 1] }
    Symmetric => {
      let t = imod(p, 2 * n)
      src[if t < n { t } else { 2 * n - 1 - t }]
    }
    Reflect =>
      if n == 1 {
        src[0]
      } else {
        let t = imod(p, 2 * n - 2)
        src[if t < n { t } else { 2 * n - 2 - t }]
      }
    Periodic | Periodization => src[imod(p, n)]
    Smooth =>
      if n == 1 {
        src[0]
      } else if p < 0 {
        src[0] + (-p).to_double() * (src[0] - src[1])
      } else {
        src[n - 1] + (p - n + 1).to_double() * (src[n - 1] - src[n - 2])
      }
    AntiSymmetric => {
      let t = imod(p, 2 * n)
      let idx = if t < n { t } else { 2 * n - 1 - t }
      let seg = if p < 0 { (-p - 1) / n } else { (p - n) / n }
      let sign = if seg % 2 == 0 { -1.0 } else { 1.0 }
      sign * src[idx]
    }
    AntiReflect =>
      if n == 1 {
        src[0]
      } else {
        // 奇反射：以 P = 2(n-1) 为反周期，f(p+P) = f(p) + 2*(x[n-1]-x[0])。
        let period = 2 * (n - 1)
        let q = ifloor_div(p, period)
        let r = p - q * period
        let base = if r < n {
          src[r]
        } else {
          2.0 * src[n - 1] - src[period - r]
        }
        base + 2.0 * q.to_double() * (src[n - 1] - src[0])
      }
  }
}
