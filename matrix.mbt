// 二维数据载体：行主序矩阵。

///|
/// 行主序矩阵（rows × cols），用于 2D 小波变换。
#derive(Eq, Debug)
pub struct Matrix {
  rows : Int32
  cols : Int32
  data : Array[Double]
}

///|
pub extend Matrix with Eq::{equal, not_equal}

///|
pub extend Matrix with @debug.Debug::{to_repr}

///|
/// 由行主序数据构造矩阵；长度必须等于 rows * cols。
pub fn Matrix::new(
  rows : Int32,
  cols : Int32,
  data : Array[Double],
) -> Matrix raise WaveletError {
  if rows <= 0 || cols <= 0 {
    raise InvalidArgument("matrix dimensions must be positive")
  }
  if data.length() != rows * cols {
    raise ShapeMismatch("data length does not match rows * cols")
  }
  { rows, cols, data, }
}

///|
/// 全零矩阵。
pub fn Matrix::zeros(rows : Int32, cols : Int32) -> Matrix raise WaveletError {
  if rows <= 0 || cols <= 0 {
    raise InvalidArgument("matrix dimensions must be positive")
  }
  { rows, cols, data: Array::make(rows * cols, 0.0), }
}

///|
/// 由若干行构造矩阵（各行长度必须一致）。
pub fn Matrix::from_rows(
  rows : Array[Array[Double]],
) -> Matrix raise WaveletError {
  if rows.length() == 0 {
    raise InvalidArgument("matrix must have at least one row")
  }
  let cols = rows[0].length()
  if cols == 0 {
    raise InvalidArgument("matrix must have at least one column")
  }
  let data : Array[Double] = []
  for row in rows {
    if row.length() != cols {
      raise ShapeMismatch("all rows must have the same length")
    }
    for v in row {
      data.push(v)
    }
  }
  { rows: rows.length(), cols, data, }
}

///|
pub fn Matrix::get(self : Matrix, r : Int32, c : Int32) -> Double {
  self.data[r * self.cols + c]
}

///|
/// 提取第 r 行。
pub fn Matrix::row(self : Matrix, r : Int32) -> Array[Double] {
  Array::makei(self.cols, c => self.data[r * self.cols + c])
}

///|
/// 提取第 c 列。
pub fn Matrix::col(self : Matrix, c : Int32) -> Array[Double] {
  Array::makei(self.rows, r => self.data[r * self.cols + c])
}

///|
/// 拆成若干行。
pub fn Matrix::to_rows(self : Matrix) -> Array[Array[Double]] {
  Array::makei(self.rows, r => self.row(r))
}

///|
/// 形状 (rows, cols)。
pub fn Matrix::shape(self : Matrix) -> (Int32, Int32) {
  (self.rows, self.cols)
}
