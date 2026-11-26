// 最小可运行示例：小波软阈值去噪（1D）。
// 运行：moon run examples/denoise
//
// 构造含噪信号 -> wavedec 多层分解 -> 对细节系数软阈值 -> waverec 重构，
// 打印去噪前后相对干净信号的 MSE，验证去噪有效。

///|
fn main {
  let n = 128
  let pi = 3.141592653589793
  let clean = Array::makei(n, i => {
    let t = i.to_double() / n.to_double()
    @math.sin(2.0 * pi * 3.0 * t) + (if i == 60 { 1.5 } else { 0.0 })
  })
  // 确定性噪声（可复现，不依赖随机源）
  let noisy = Array::makei(n, i => {
    clean[i] +
    0.35 *
    @math.sin(i.to_double() * 12.9898) *
    @math.cos(i.to_double() * 78.233)
  })
  let w = @lib.Wavelet::new("db4")
  let coeffs = @lib.wavedec(noisy, w, level=Some(3))
  // 保留近似支（coeffs[0]），对各层细节支软阈值
  let shrunk : Array[Array[Double]] = Array::makei(coeffs.length(), k => {
    if k == 0 {
      coeffs[0]
    } else {
      @lib.threshold(coeffs[k], 0.3, mode=@lib.Soft)
    }
  })
  let rec = @lib.waverec(shrunk, w)
  let before = mse(noisy, clean, n)
  let after = mse(rec, clean, n)
  println("noisy   MSE = \{before}")
  println("denoise MSE = \{after}")
  println("improved    = \{after < before}")
}

///|
fn mse(a : Array[Double], b : Array[Double], n : Int32) -> Double {
  let mut acc = 0.0
  for i in 0..<n {
    let d = a[i] - b[i]
    acc += d * d
  }
  acc / n.to_double()
}
