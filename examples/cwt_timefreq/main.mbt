// 最小可运行示例：连续小波变换时频分析（CWT）。
// 运行：moon run examples/cwt_timefreq
//
// 对线性调频（chirp）信号做 mexh 连续小波变换，打印各尺度对应的频率与
// 该尺度下的最大幅值，展示能量随尺度（频率）的分布。

///|
fn main {
  let n = 64
  let pi = 3.141592653589793
  let x = Array::makei(n, i => {
    let t = i.to_double() / n.to_double()
    @math.sin(2.0 * pi * 5.0 * t) + @math.sin(2.0 * pi * 20.0 * t * t)
  })
  let w = @lib.Wavelet::new("mexh")
  let scales = [1.0, 2.0, 4.0, 8.0, 16.0]
  let (coefs, freqs) = @lib.cwt(x, scales, w)
  for s in 0..<scales.length() {
    let mut mx = 0.0
    for v in coefs.row(s) {
      let a = if v < 0.0 { -v } else { v }
      if a > mx {
        mx = a
      }
    }
    println("scale \{scales[s]}  freq \{freqs[s]}  max|coef| \{mx}")
  }
}
