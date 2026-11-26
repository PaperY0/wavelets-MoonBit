// 内部 radix-2 FFT（仅 2 的幂长度），供 central_frequency 使用。
// 不公开；白盒测试见 fft_wbtest.mbt。

///|
/// 原地迭代 FFT。re/im 长度必须为 2 的幂。
fn fft_inplace(re : Array[Double], im : Array[Double]) -> Unit {
  let n = re.length()
  // 位反转置换。
  let mut j = 0
  for i in 0..<(n - 1) {
    if i < j {
      let tr = re[i]
      re[i] = re[j]
      re[j] = tr
      let ti = im[i]
      im[i] = im[j]
      im[j] = ti
    }
    let mut k = n / 2
    while k <= j {
      j -= k
      k /= 2
    }
    j += k
  }
  // 蝶形运算。
  let mut len = 2
  while len <= n {
    let ang = -2.0 * 3.141592653589793 / len.to_double()
    let wr = @math.cos(ang)
    let wi = @math.sin(ang)
    let half = len / 2
    let mut start = 0
    while start < n {
      let mut cur_r = 1.0
      let mut cur_i = 0.0
      for k in 0..<half {
        let a = start + k
        let b = a + half
        let tr = cur_r * re[b] - cur_i * im[b]
        let ti = cur_r * im[b] + cur_i * re[b]
        re[b] = re[a] - tr
        im[b] = im[a] - ti
        re[a] = re[a] + tr
        im[a] = im[a] + ti
        let next_r = cur_r * wr - cur_i * wi
        cur_i = cur_r * wi + cur_i * wr
        cur_r = next_r
      }
      start += len
    }
    len *= 2
  }
}

///|
/// 实序列 FFT，返回 (re, im)。
fn fft_real(x : Array[Double]) -> (Array[Double], Array[Double]) {
  let n = x.length()
  let re = Array::makei(n, i => x[i])
  let im = Array::make(n, 0.0)
  fft_inplace(re, im)
  (re, im)
}
