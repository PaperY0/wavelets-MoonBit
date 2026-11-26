// 尺度函数 / 小波函数的级联近似（wavefun）。
//
// 移植来源：PyWavelets `pywt/_extensions/_pywt.pyx::wavefun`（MIT）。
// 算法：以 (sqrt(2))^level 为单系数做 upcoef 级联，前端补一个 0，
// 按 keep_length 截断、右端补零到 output_length；双正交小波对
// 分解侧使用逆滤波器组，psi 的符号由消失矩数目决定。

///|
/// wavefun 结果。连续小波只有 psi；正交小波有 phi/psi；
/// 双正交小波另有重构侧 phi_r/psi_r（phi/psi 为分解侧）。
#derive(Eq, Debug)
pub struct Wavefun {
  phi : Array[Double]?
  psi : Array[Double]
  phi_r : Array[Double]?
  psi_r : Array[Double]?
  x : Array[Double]
}

///|
pub extend Wavefun with Eq::{equal, not_equal}

///|
pub extend Wavefun with @debug.Debug::{to_repr}

///|
/// 计算小波函数近似。level 为细化层数（网格 2^level），上限 20
/// 以防止内存爆炸（上游无硬上限，见 known-deviations）。
pub fn Wavelet::wavefun(
  self : Wavelet,
  level~ : Int32 = 8,
) -> Wavefun raise WaveletError {
  if level < 1 || level > 20 {
    raise InvalidLevel("wavefun level must be in 1..=20")
  }
  match self.kind {
    Continuous => {
      let n = 1 << level
      let x = linspace(self.lower_bound, self.upper_bound, n)
      let psi = Array::makei(n, i => continuous_psi(self.name, x[i]))
      { phi: None, psi, phi_r: None, psi_r: None, x, }
    }
    Discrete =>
      if self.orthogonal {
        let n = @math.pow(2.0, level.to_double() / 2.0)
        let (phi, psi, out_len) = wavefun_side(self, level, n, true)
        let p2 = (1 << level).to_double()
        let x = linspace(0.0, (out_len - 1).to_double() / p2, out_len)
        { phi: Some(phi), psi, phi_r: None, psi_r: None, x, }
      } else if self.biorthogonal {
        let n = @math.pow(2.0, level.to_double() / 2.0)
        let n_mul = if self.vm_psi % 4 != 1 { -n } else { n }
        // 分解侧：逆滤波器组构成的小波。
        let (idl, idh, irl, irh) = self.inverse_filter_bank()
        let other : Wavelet = {
          name: self.name,
          family_name: self.family_name,
          short_name: self.short_name,
          kind: Discrete,
          dec_lo: idl,
          dec_hi: idh,
          rec_lo: irl,
          rec_hi: irh,
          orthogonal: false,
          biorthogonal: true,
          symmetry: self.symmetry,
          vm_psi: self.vm_psi,
          vm_phi: self.vm_phi,
          lower_bound: 0.0,
          upper_bound: 0.0,
        }
        let (phi_d, _, _) = wavefun_side(other, level, n, false)
        let (_, psi_d_mul, _) = wavefun_side(other, level, n_mul, false)
        let (phi_r, _, out_len) = wavefun_side(self, level, n, false)
        let (_, psi_r_mul, _) = wavefun_side(self, level, n_mul, false)
        let p2 = (1 << level).to_double()
        let x = linspace(0.0, (out_len - 1).to_double() / p2, out_len)
        {
          phi: Some(phi_d),
          psi: psi_d_mul,
          phi_r: Some(phi_r),
          psi_r: Some(psi_r_mul),
          x,
        }
      } else {
        raise NotSupported(
          "wavefun requires an orthogonal or biorthogonal wavelet",
        )
      }
  }
}

///|
/// 单侧级联：返回 (phi, psi, output_length)。
/// plus_one 区分正交（(F-1)*2^level + 1）与双正交（(F-1)*2^level）。
fn wavefun_side(
  w : Wavelet,
  level : Int32,
  n_coef : Double,
  plus_one : Bool,
) -> (Array[Double], Array[Double], Int32) raise WaveletError {
  let p2 = 1 << level
  let f = w.dec_len()
  let mut output_length = (f - 1) * p2 + (if plus_one { 1 } else { 0 })
  let mut keep = 1
  for _ in 0..<level {
    keep = 2 * keep + (f - 2)
  }
  if output_length - keep - 2 < 0 {
    output_length = keep + 2
  }
  let right = output_length - keep - 1
  let phi_raw = upcoef(Approx, [n_coef], w, level~)
  let psi_raw = upcoef(Detail, [n_coef], w, level~)
  let phi = assemble_wavefun(phi_raw, keep, right)
  let psi = assemble_wavefun(psi_raw, keep, right)
  (phi, psi, output_length)
}

///|
/// [0] + raw[:keep] + zeros(right)。
fn assemble_wavefun(
  raw : Array[Double],
  keep : Int32,
  right : Int32,
) -> Array[Double] {
  let k = if keep < raw.length() { keep } else { raw.length() }
  Array::makei(1 + k + right, i => {
    if i == 0 {
      0.0
    } else if i <= k {
      raw[i - 1]
    } else {
      0.0
    }
  })
}
