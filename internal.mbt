// 内部数值工具：数组拼接、取模、卷积与上采样。
// 这些 helper 不公开；白盒测试见 internal_wbtest.mbt。

///|
/// 数组拼接（不修改输入）。
fn concat_arrays(a : Array[Double], b : Array[Double]) -> Array[Double] {
  let na = a.length()
  let nb = b.length()
  Array::makei(na + nb, i => if i < na { a[i] } else { b[i - na] })
}

///|
/// 绝对值。
fn dabs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

///|
/// 非负模：结果恒在 [0, m)。
fn imod(a : Int32, m : Int32) -> Int32 {
  (a % m + m) % m
}

///|
/// 向下取整除法（对负数也成立）。
fn ifloor_div(a : Int32, m : Int32) -> Int32 {
  (a - imod(a, m)) / m
}

///|
/// 上采样 2 倍：[a0, a1, ...] -> [a0, 0, a1, 0, ..., a_{L-1}, 0]，长度 2L。
/// 与 pywt 的 upsampling 约定一致（末尾补零）。
fn upsample2(a : Array[Double]) -> Array[Double] {
  let n = a.length()
  Array::makei(2 * n, i => if i % 2 == 0 { a[i / 2] } else { 0.0 })
}

///|
/// 下采样卷积（DWT 分析端）：
/// `out[i] = sum_j filt[F-1-j] * ext[2*i + j]`，i ∈ [0, out_len)。
/// 即滤波器反转后与扩展信号做 stride-2 的有效卷积，与 pywt 的
/// downsampling convolution 相位一致。
fn conv_down(
  ext : Array[Double],
  filt : Array[Double],
  out_len : Int32,
) -> Array[Double] {
  let f = filt.length()
  Array::makei(out_len, i => {
    let mut acc = 0.0
    for j in 0..<f {
      acc += filt[f - 1 - j] * ext[2 * i + j]
    }
    acc
  })
}

///|
/// 上采样卷积的单步（upcoef 级联用）：
/// `out[n] = sum_k filt[k] * up[n-k]`，n ∈ [0, 2L + F - 2)。
/// 即完整卷积去掉最后一个样本，与 pywt upcoef 的输出长度约定一致。
fn upcoef_step(up : Array[Double], filt : Array[Double]) -> Array[Double] {
  let lu = up.length()
  let f = filt.length()
  let out_len = lu + f - 2
  Array::makei(out_len, n => {
    let mut acc = 0.0
    for k in 0..<f {
      let idx = n - k
      if idx >= 0 && idx < lu {
        acc += filt[k] * up[idx]
      }
    }
    acc
  })
}

///|
/// 完整卷积：`out[n] = sum_k a[k] * b[n-k]`，长度 la + lb - 1。
fn conv_full(a : Array[Double], b : Array[Double]) -> Array[Double] {
  let la = a.length()
  let lb = b.length()
  Array::makei(la + lb - 1, n => {
    let mut acc = 0.0
    for k in 0..<la {
      let j = n - k
      if j >= 0 && j < lb {
        acc += a[k] * b[j]
      }
    }
    acc
  })
}

///|
/// 等距网格 linspace(start, stop, n)（含端点，与 numpy.linspace 一致）。
fn linspace(start : Double, stop : Double, n : Int32) -> Array[Double] {
  if n == 1 {
    return [start]
  }
  let step = (stop - start) / (n - 1).to_double()
  Array::makei(n, i => start + i.to_double() * step)
}
