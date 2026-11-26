// 实值连续小波的解析表达式。
//
// 移植来源：PyWavelets C 端 wavelets.c 的 mexh / morl / gaus 族（MIT）。
// gausN 为 exp(-t^2) 的 N 阶导数（物理学家 Hermite 多项式表示），
// 按连续 L2 范数归一化；归一化常数由 tools/probe.py 对上游 wavefun
// 精确拟合（17 位有效数字，与网格无关）。

///|
/// mexh 归一化常数 2 / (sqrt(3) * pi^(1/4))。
let mexh_norm : Double = 0.8673250705840776

///|
/// gaus1..gaus8 的归一化常数（含符号），拟合自 pywt wavefun。
let gaus_norms : Array[Double] = [
  0.8932438417380022, -0.51571457247941122, -0.23063456821024109, 0.087171673031294444,
  0.029057224343764818, -0.0087610827815862557, -0.0024298871690456858, 0.00062739416925846833,
]

///|
/// 物理学家 Hermite 多项式 H_n(t)：H_0 = 1, H_1 = 2t,
/// H_k = 2t H_{k-1} - 2(k-1) H_{k-2}。
fn hermite(n : Int32, t : Double) -> Double {
  if n == 0 {
    return 1.0
  }
  let mut h_prev = 1.0
  let mut h = 2.0 * t
  for k in 1..<n {
    let next = 2.0 * t * h - 2.0 * k.to_double() * h_prev
    h_prev = h
    h = next
  }
  h
}

///|
/// 计算实值连续小波在 t 处的 psi 值。name 必须是注册表中的连续小波。
fn continuous_psi(name : String, t : Double) -> Double raise WaveletError {
  match name {
    "mexh" => mexh_norm * (1.0 - t * t) * @math.exp(-t * t / 2.0)
    "morl" => @math.exp(-t * t / 2.0) * @math.cos(5.0 * t)
    _ =>
      if name.has_prefix("gaus") {
        let order = parse_gaus_order(name)
        let sign = if order % 2 == 0 { 1.0 } else { -1.0 }
        gaus_norms[order - 1] * sign * hermite(order, t) * @math.exp(-t * t)
      } else {
        raise UnknownWavelet(name)
      }
  }
}

///|
fn parse_gaus_order(name : String) -> Int32 raise WaveletError {
  let suffix = name[4:].to_owned()
  let mut value = 0
  if suffix.length() == 0 {
    raise UnknownWavelet(name)
  }
  for c in suffix {
    if !c.is_ascii_digit() {
      raise UnknownWavelet(name)
    }
    value = value * 10 + (c.to_int() - '0'.to_int())
  }
  if value < 1 || value > 8 {
    raise UnknownWavelet(name)
  }
  value
}
