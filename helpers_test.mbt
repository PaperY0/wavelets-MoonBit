// 黑盒测试辅助：黄金 JSON 的读取与浮点容差比较。
// 期望值中的 null 表示 NaN（上游 garrote 在 t=0、x=0 处的行为）。

///|
using @json {type Json}

///|
let golden_root : Json = @json.parse(golden_json[:]) catch { _ => panic() }

///|
fn jget(j : Json, key : String) -> Json {
  match j {
    Object(m) => m.get(key).unwrap()
    _ => panic()
  }
}

///|
fn jnum(j : Json) -> Double {
  match j {
    Number(v, ..) => v
    _ => panic()
  }
}

///|
fn jstr(j : Json) -> String {
  match j {
    String(s) => s
    _ => panic()
  }
}

///|
fn jarr(j : Json) -> Array[Json] {
  match j {
    Array(a) => a
    _ => panic()
  }
}

///|
fn jint(j : Json) -> Int32 {
  jnum(j).to_int()
}

///|
/// 数值数组；null -> NaN。
fn jnums(j : Json) -> Array[Double] {
  let a = jarr(j)
  Array::makei(a.length(), i => {
    match a[i] {
      Null => 0.0 / 0.0
      Number(v, ..) => v
      _ => panic()
    }
  })
}

///|
fn jmat(j : Json) -> Array[Array[Double]] {
  let rows = jarr(j)
  Array::makei(rows.length(), i => jnums(rows[i]))
}

///|
fn tabs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

///|
fn arr_close(a : Array[Double], b : Array[Double], eps : Double) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i in 0..<a.length() {
    let x = a[i]
    let y = b[i]
    if y.is_nan() {
      if !x.is_nan() {
        return false
      }
    } else if x.is_nan() || tabs(x - y) > eps {
      return false
    }
  }
  true
}

///|
fn mat_close(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
  eps : Double,
) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i in 0..<a.length() {
    if !arr_close(a[i], b[i], eps) {
      return false
    }
  }
  true
}

///|
fn try_new(name : String) -> @wavelets.Wavelet? {
  Some(@wavelets.Wavelet::new(name)) catch {
    _ => None
  }
}
