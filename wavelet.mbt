// 小波对象与内置注册表。
//
// 移植来源：PyWavelets `pywt/_extensions/_pywt.pyx` 的 Wavelet /
// ContinuousWavelet 对象模型（MIT）。滤波器系数数据由 tools/gen_wavelets.py
// 从上游精确导出（wavelet_data.mbt，17 位有效数字，bit-identical）。

///|
/// 小波类别。
#derive(Eq, Debug, Show)
pub(all) enum WaveletKind {
  Discrete
  Continuous
}

///|
/// 对称性（对应 pywt 的 symmetry 属性）。
#derive(Eq, Debug)
pub(all) enum Symmetry {
  Asymmetric
  NearSymmetric
  Symmetric
  AntiSymmetric
}

///|
/// 对称性的 pywt 字符串表示。
pub fn Symmetry::name(self : Symmetry) -> String {
  match self {
    Asymmetric => "asymmetric"
    NearSymmetric => "near symmetric"
    Symmetric => "symmetric"
    AntiSymmetric => "anti-symmetric"
  }
}

///|
/// 一个小波：离散小波携带四个滤波器组，连续小波携带支撑区间。
///
/// 与 pywt 一致：`dec_lo/dec_hi` 为分解低通/高通，`rec_lo/rec_hi` 为重构
/// 低通/高通；连续小波的滤波器组为空数组。
#derive(Eq, Debug)
pub struct Wavelet {
  name : String
  family_name : String
  short_name : String
  kind : WaveletKind
  dec_lo : Array[Double]
  dec_hi : Array[Double]
  rec_lo : Array[Double]
  rec_hi : Array[Double]
  orthogonal : Bool
  biorthogonal : Bool
  symmetry : Symmetry
  /// psi 消失矩数目（0 表示未知/不适用）
  vm_psi : Int32
  /// phi 消失矩数目（0 表示未知/不适用）
  vm_phi : Int32
  lower_bound : Double
  upper_bound : Double
}

///|
pub extend Wavelet with Eq::{equal, not_equal}

///|
pub extend Wavelet with @debug.Debug::{to_repr}

///|
/// 按名称构造小波（如 "db4"、"bior3.5"、"mexh"）。名称区分大小写，
/// 与 pywt 一致；"db1" 与 "haar" 均为合法名称。
pub fn Wavelet::new(name : String) -> Wavelet raise WaveletError {
  for spec in discrete_specs {
    if spec.name == name {
      return {
        name: spec.name,
        family_name: spec.family_name,
        short_name: spec.short_name,
        kind: Discrete,
        dec_lo: spec.dec_lo,
        dec_hi: spec.dec_hi,
        rec_lo: spec.rec_lo,
        rec_hi: spec.rec_hi,
        orthogonal: spec.orthogonal,
        biorthogonal: spec.biorthogonal,
        symmetry: spec.symmetry,
        vm_psi: spec.vm_psi,
        vm_phi: spec.vm_phi,
        lower_bound: 0.0,
        upper_bound: 0.0,
      }
    }
  }
  for spec in continuous_specs {
    if spec.name == name {
      return {
        name: spec.name,
        family_name: spec.family_name,
        short_name: spec.short_name,
        kind: Continuous,
        dec_lo: [],
        dec_hi: [],
        rec_lo: [],
        rec_hi: [],
        orthogonal: false,
        biorthogonal: false,
        symmetry: spec.symmetry,
        vm_psi: 0,
        vm_phi: 0,
        lower_bound: spec.lower_bound,
        upper_bound: spec.upper_bound,
      }
    }
  }
  raise UnknownWavelet(name)
}

///|
/// 分解滤波器长度。
pub fn Wavelet::dec_len(self : Wavelet) -> Int32 {
  self.dec_lo.length()
}

///|
/// 重构滤波器长度。
pub fn Wavelet::rec_len(self : Wavelet) -> Int32 {
  self.rec_lo.length()
}

///|
/// 滤波器组 (dec_lo, dec_hi, rec_lo, rec_hi)。
pub fn Wavelet::filter_bank(
  self : Wavelet,
) -> (Array[Double], Array[Double], Array[Double], Array[Double]) {
  (self.dec_lo, self.dec_hi, self.rec_lo, self.rec_hi)
}

///|
/// 逆滤波器组 (rec_lo[::-1], rec_hi[::-1], dec_lo[::-1], dec_hi[::-1])，
/// 双正交小波的 wavefun 会用到（对应 pywt inverse_filter_bank）。
pub fn Wavelet::inverse_filter_bank(
  self : Wavelet,
) -> (Array[Double], Array[Double], Array[Double], Array[Double]) {
  (
    reverse_array(self.rec_lo),
    reverse_array(self.rec_hi),
    reverse_array(self.dec_lo),
    reverse_array(self.dec_hi),
  )
}

///|
fn reverse_array(a : Array[Double]) -> Array[Double] {
  let n = a.length()
  Array::makei(n, i => a[n - 1 - i])
}

///|
/// 内置小波族短名（v0.1 支持的实值小波族）。
pub fn families() -> Array[String] {
  ["haar", "db", "sym", "coif", "bior", "rbio", "dmey", "gaus", "mexh", "morl"]
}

///|
/// 列出内置小波名。`kind` 省略时返回全部（离散在前、连续在后），
/// 顺序与 pywt.wavelist 一致（复值小波族除外，见 known-deviations）。
pub fn wavelist(kind~ : WaveletKind? = None) -> Array[String] {
  let out : Array[String] = []
  match kind {
    None | Some(Discrete) =>
      for spec in discrete_specs {
        out.push(spec.name)
      }
    _ => ()
  }
  match kind {
    None | Some(Continuous) =>
      for spec in continuous_specs {
        out.push(spec.name)
      }
    _ => ()
  }
  out
}
