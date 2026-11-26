// 二维离散小波变换（可分离：先行后列）。
//
// 移植来源：PyWavelets `pywt/_multidim.py::dwt2/idwt2` 与
// `_multilevel.py::wavedec2/waverec2`（MIT）。
// 方向约定（已对上游验证）：先沿行（axis=1）变换，再沿列（axis=0）；
// cA=LL, cH=LH（水平细节）, cV=HL（垂直细节）, cD=HH。

///|
/// 一层 2D 细节系数（三个方向）。
#derive(Eq, Debug)
pub struct Detail2 {
  horizontal : Matrix
  vertical : Matrix
  diagonal : Matrix
}

///|
pub extend Detail2 with Eq::{equal, not_equal}

///|
pub extend Detail2 with @debug.Debug::{to_repr}

///|
/// 单层 2D 变换结果；任一分量可为 None（idwt2 时等价补零）。
#derive(Eq, Debug)
pub struct Coeffs2 {
  approx : Matrix?
  horizontal : Matrix?
  vertical : Matrix?
  diagonal : Matrix?
}

///|
pub extend Coeffs2 with Eq::{equal, not_equal}

///|
pub extend Coeffs2 with @debug.Debug::{to_repr}

///|
/// 多层 2D 分解结果：最终近似 + 各层细节（从粗到细）。
#derive(Eq, Debug)
pub struct Decomposition2 {
  approx : Matrix
  details : Array[Detail2]
}

///|
pub extend Decomposition2 with Eq::{equal, not_equal}

///|
pub extend Decomposition2 with @debug.Debug::{to_repr}

///|
/// 单层 2D 离散小波变换。
pub fn dwt2(
  data : Matrix,
  wavelet : Wavelet,
  mode~ : Mode = Symmetric,
) -> Coeffs2 raise WaveletError {
  guard wavelet.kind is Discrete else { raise NotDiscrete(wavelet.name) }
  let (rows, cols) = data.shape()
  let f = wavelet.dec_len()
  let oc = dwt_coeff_len(cols, f, mode)
  let or_ = dwt_coeff_len(rows, f, mode)
  // 先沿行变换。
  let lo = Matrix::zeros(rows, oc)
  let hi = Matrix::zeros(rows, oc)
  for r in 0..<rows {
    let (a, d) = dwt(data.row(r), wavelet, mode~)
    for c in 0..<oc {
      lo.data[r * oc + c] = a[c]
      hi.data[r * oc + c] = d[c]
    }
  }
  // 再沿列变换。
  let ll = Matrix::zeros(or_, oc)
  let lh = Matrix::zeros(or_, oc)
  let hl = Matrix::zeros(or_, oc)
  let hh = Matrix::zeros(or_, oc)
  for c in 0..<oc {
    let (a, d) = dwt(lo.col(c), wavelet, mode~)
    let (b, e) = dwt(hi.col(c), wavelet, mode~)
    for r in 0..<or_ {
      ll.data[r * oc + c] = a[r]
      lh.data[r * oc + c] = d[r]
      hl.data[r * oc + c] = b[r]
      hh.data[r * oc + c] = e[r]
    }
  }
  {
    approx: Some(ll),
    horizontal: Some(lh),
    vertical: Some(hl),
    diagonal: Some(hh),
  }
}

///|
/// 单层 2D 逆离散小波变换。分量可为 None（补零），但给出的分量
/// 形状必须一致；重构顺序为先列后行。
pub fn idwt2(
  coeffs : Coeffs2,
  wavelet : Wavelet,
  mode~ : Mode = Symmetric,
) -> Matrix raise WaveletError {
  guard wavelet.kind is Discrete else { raise NotDiscrete(wavelet.name) }
  // 确定系数形状。
  let shape = match
    (coeffs.approx, coeffs.horizontal, coeffs.vertical, coeffs.diagonal) {
    (Some(m), _, _, _)
    | (_, Some(m), _, _)
    | (_, _, Some(m), _)
    | (_, _, _, Some(m)) => m.shape()
    _ => raise InvalidArgument("at least one coefficient matrix is required")
  }
  let (r0, c0) = shape
  let check = (m : Matrix?) => {
    match m {
      Some(mm) =>
        if mm.shape() != shape {
          raise ShapeMismatch("all coefficient matrices must share one shape")
        }
      None => ()
    }
  }
  check(coeffs.approx)
  check(coeffs.horizontal)
  check(coeffs.vertical)
  check(coeffs.diagonal)
  let zeros = Matrix::zeros(r0, c0)
  let ll = coeffs.approx.unwrap_or(zeros)
  let lh = coeffs.horizontal.unwrap_or(zeros)
  let hl = coeffs.vertical.unwrap_or(zeros)
  let hh = coeffs.diagonal.unwrap_or(zeros)
  let f = wavelet.rec_len()
  let per = mode is Periodization
  let out_rows = if per { 2 * r0 } else { 2 * r0 - f + 2 }
  let out_cols = if per { 2 * c0 } else { 2 * c0 - f + 2 }
  if out_rows <= 0 || out_cols <= 0 {
    raise InvalidArgument("coefficients too short for this wavelet")
  }
  // 先沿列（axis 0）重构：LL+LH -> lo（行低通支），HL+HH -> hi（行高通支）。
  let lo = Matrix::zeros(out_rows, c0)
  let hi = Matrix::zeros(out_rows, c0)
  for c in 0..<c0 {
    let a = idwt(Some(ll.col(c)), Some(lh.col(c)), wavelet, mode~)
    let b = idwt(Some(hl.col(c)), Some(hh.col(c)), wavelet, mode~)
    for r in 0..<out_rows {
      lo.data[r * c0 + c] = a[r]
      hi.data[r * c0 + c] = b[r]
    }
  }
  // 再沿行重构。
  let out = Matrix::zeros(out_rows, out_cols)
  for r in 0..<out_rows {
    let row = idwt(Some(lo.row(r)), Some(hi.row(r)), wavelet, mode~)
    for c in 0..<out_cols {
      out.data[r * out_cols + c] = row[c]
    }
  }
  out
}

///|
/// 多层 2D 分解。level 省略时取两个方向最大有用层数的较小值；
/// 超过最大值时报错（与 1D 一致的上游差异，见 known-deviations）。
pub fn wavedec2(
  data : Matrix,
  wavelet : Wavelet,
  mode~ : Mode = Symmetric,
  level~ : Int32? = None,
) -> Decomposition2 raise WaveletError {
  guard wavelet.kind is Discrete else { raise NotDiscrete(wavelet.name) }
  let (rows, cols) = data.shape()
  let f = wavelet.dec_len()
  let max_lvl = {
    let a = dwt_max_level(rows, f)
    let b = dwt_max_level(cols, f)
    if a < b {
      a
    } else {
      b
    }
  }
  let lvl = match level {
    None => max_lvl
    Some(l) => l
  }
  if lvl < 0 {
    raise InvalidLevel("level must be >= 0")
  }
  if lvl > max_lvl {
    raise InvalidLevel(
      "level \{lvl} exceeds maximum useful level \{max_lvl} for this matrix/wavelet",
    )
  }
  let mut a = data
  let details : Array[Detail2] = []
  for _ in 0..<lvl {
    let c = dwt2(a, wavelet, mode~)
    let d = {
      horizontal: c.horizontal.unwrap(),
      vertical: c.vertical.unwrap(),
      diagonal: c.diagonal.unwrap(),
    }
    details.push(d)
    a = c.approx.unwrap()
  }
  // details 目前从细到粗，反转为从粗到细。
  let n = details.length()
  let ordered : Array[Detail2] = Array::makei(n, i => details[n - 1 - i])
  { approx: a, details: ordered, }
}

///|
/// 多层 2D 重构。裁剪规则与 pywt.waverec2 一致：每步之前，
/// 若近似矩阵某维比细节矩阵长 1，则截去该维末尾一行/列。
pub fn waverec2(
  dec : Decomposition2,
  wavelet : Wavelet,
  mode~ : Mode = Symmetric,
) -> Matrix raise WaveletError {
  guard wavelet.kind is Discrete else { raise NotDiscrete(wavelet.name) }
  let mut a = dec.approx
  for d in dec.details {
    let (ar, ac) = a.shape()
    let (dr, dc) = d.horizontal.shape()
    if d.vertical.shape() != (dr, dc) || d.diagonal.shape() != (dr, dc) {
      raise ShapeMismatch("detail matrices must share one shape")
    }
    let rows = if ar == dr + 1 {
      dr
    } else if ar == dr {
      ar
    } else {
      raise ShapeMismatch("approximation/detail row mismatch")
    }
    let cols = if ac == dc + 1 {
      dc
    } else if ac == dc {
      ac
    } else {
      raise ShapeMismatch("approximation/detail column mismatch")
    }
    let trimmed = Matrix::new(
      rows,
      cols,
      Array::makei(rows * cols, i => a.get(i / cols, i % cols)),
    )
    a = idwt2(
      {
        approx: Some(trimmed),
        horizontal: Some(d.horizontal),
        vertical: Some(d.vertical),
        diagonal: Some(d.diagonal),
      },
      wavelet,
      mode~,
    )
  }
  a
}
