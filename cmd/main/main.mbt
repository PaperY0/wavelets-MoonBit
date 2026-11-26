// wavelets CLI：查询小波元数据、执行变换并输出 JSON。
//
// 用法：
// ```bash
// moon run cmd/main -- list [discrete|continuous]
// moon run cmd/main -- info db2
// moon run cmd/main -- dwt db2 symmetric 1 2 3 4 5 6 7 8
// moon run cmd/main -- wavedec db2 symmetric 2 1 2 3 4 5 6 7 8
// moon run cmd/main -- cwt mexh 1,2,4,8 1 2 3 4 5 6 7 8
// moon run cmd/main -- wavefun db2 5
// moon run cmd/main -- threshold soft 2 -3 -2 -1 0 1 2 3
// moon run cmd/main -- pad antisymmetric 4 4 1 2 3
// ```
//
// 数值结果一律输出 JSON，供 tools/diff_test.py 与 PyWavelets 差分对拍。

///|
using @json {type Json}

///|
fn main {
  run(cli_args()) catch {
    e => {
      eprintln("error: \{e}")
      panic()
    }
  }
}

///|
#cfg(target="js")
fn cli_args() -> Array[String] {
  @env.args()[1:].to_owned()
}

///|
#cfg(not(target="js"))
fn cli_args() -> Array[String] {
  @env.args()
}

///|
fn usage() -> String {
  "usage: moon run cmd/main -- <command> [args]\n" +
  "commands:\n" +
  "  list [discrete|continuous]        list built-in wavelet names\n" +
  "  info <wavelet>                    show metadata and filter bank\n" +
  "  dwt <wavelet> <mode> <v...>       single-level DWT -> {cA, cD}\n" +
  "  wavedec <wavelet> <mode> <level> <v...>  multilevel -> {coeffs}\n" +
  "  cwt <wavelet> <scales-csv> <v...> continuous WT -> {coefs, freqs}\n" +
  "  wavefun <wavelet> <level>         scaling/wavelet functions\n" +
  "  threshold <mode> <value> <v...>   soft|hard|greater|less|garrote\n" +
  "  pad <mode> <left> <right> <v...>  signal extension"
}

///|
fn need(
  args : Array[String],
  i : Int32,
  what : String,
) -> String raise @lib.WaveletError {
  if i < args.length() {
    args[i]
  } else {
    raise @lib.InvalidArgument("missing \{what}")
  }
}

///|
fn parse_num(s : String) -> Double raise @lib.WaveletError {
  @string.parse_double(s[:]) catch {
    _ => raise @lib.InvalidArgument("not a number: \{s}")
  }
}

///|
fn parse_int(s : String) -> Int32 raise @lib.WaveletError {
  let v = parse_num(s)
  if v != v.to_int().to_double() {
    raise @lib.InvalidArgument("expected an integer: \{s}")
  }
  v.to_int()
}

///|
fn values_from(
  args : Array[String],
  from : Int32,
) -> Array[Double] raise @lib.WaveletError {
  if from >= args.length() {
    raise @lib.InvalidArgument("at least one sample value is required")
  }
  Array::makei(args.length() - from, i => parse_num(args[from + i]))
}

///|
fn split_csv(s : String) -> Array[String] {
  let out : Array[String] = []
  let mut cur = ""
  for c in s {
    if c == ',' {
      out.push(cur)
      cur = ""
    } else {
      cur = cur + c.to_string()
    }
  }
  out.push(cur)
  out
}

///|
fn emit(j : Json) -> Unit {
  println(j.stringify())
}

///|
fn obj(pairs : Array[(String, Json)]) -> Json {
  Json::object(Map::from_array(pairs))
}

///|
fn run(args : Array[String]) -> Unit raise @lib.WaveletError {
  if args.length() < 2 {
    println(usage())
    return
  }
  match args[1] {
    "--help" | "-h" | "help" => println(usage())
    "list" => {
      let kind = if args.length() > 2 {
        match args[2] {
          "discrete" => Some(@lib.Discrete)
          "continuous" => Some(@lib.Continuous)
          _ =>
            raise @lib.InvalidArgument("list kind must be discrete|continuous")
        }
      } else {
        None
      }
      for name in @lib.wavelist(kind~) {
        println(name)
      }
    }
    "info" => print_info(@lib.Wavelet::new(need(args, 2, "wavelet name")))
    "dwt" => {
      let w = @lib.Wavelet::new(need(args, 2, "wavelet name"))
      let mode = @lib.Mode::parse(need(args, 3, "mode"))
      let (a, d) = @lib.dwt(values_from(args, 4), w, mode~)
      emit(obj([("cA", @json.to_json(a)), ("cD", @json.to_json(d))]))
    }
    "wavedec" => {
      let w = @lib.Wavelet::new(need(args, 2, "wavelet name"))
      let mode = @lib.Mode::parse(need(args, 3, "mode"))
      let level = parse_int(need(args, 4, "level"))
      let coeffs = @lib.wavedec(
        values_from(args, 5),
        w,
        mode~,
        level=Some(level),
      )
      emit(obj([("coeffs", @json.to_json(coeffs))]))
    }
    "cwt" => {
      let w = @lib.Wavelet::new(need(args, 2, "wavelet name"))
      let scale_strs = split_csv(need(args, 3, "scales (comma separated)"))
      let scales : Array[Double] = []
      for s in scale_strs {
        scales.push(parse_num(s))
      }
      let (coefs, freqs) = @lib.cwt(values_from(args, 4), scales, w)
      emit(
        obj([
          ("coefs", @json.to_json(coefs.to_rows())),
          ("freqs", @json.to_json(freqs)),
        ]),
      )
    }
    "wavefun" => {
      let w = @lib.Wavelet::new(need(args, 2, "wavelet name"))
      let level = parse_int(need(args, 3, "level"))
      let wf = w.wavefun(level~)
      let opt_json = (a : Array[Double]?) => {
        match a {
          Some(xs) => @json.to_json(xs)
          None => Json::null()
        }
      }
      emit(
        obj([
          ("phi", opt_json(wf.phi)),
          ("psi", @json.to_json(wf.psi)),
          ("phi_r", opt_json(wf.phi_r)),
          ("psi_r", opt_json(wf.psi_r)),
          ("x", @json.to_json(wf.x)),
        ]),
      )
    }
    "threshold" => {
      let mode = match need(args, 2, "threshold mode") {
        "soft" => @lib.Soft
        "hard" => @lib.Hard
        "greater" => @lib.Greater
        "less" => @lib.Less
        "garrote" => @lib.Garrote
        other => raise @lib.InvalidArgument("unknown threshold mode: \{other}")
      }
      let t = parse_num(need(args, 3, "threshold value"))
      let r = @lib.threshold(values_from(args, 4), t, mode~)
      emit(obj([("r", @json.to_json(r))]))
    }
    "pad" => {
      let mode = @lib.Mode::parse(need(args, 2, "mode"))
      let left = parse_int(need(args, 3, "left pad width"))
      let right = parse_int(need(args, 4, "right pad width"))
      let r = @lib.pad(values_from(args, 5), left, right, mode)
      emit(obj([("r", @json.to_json(r))]))
    }
    other => raise @lib.InvalidArgument("unknown command: \{other}\n\{usage()}")
  }
}

///|
fn print_info(w : @lib.Wavelet) -> Unit {
  println("name        : \{w.name}")
  println("family      : \{w.family_name} (\{w.short_name})")
  println("kind        : \{w.kind}")
  println("symmetry    : \{w.symmetry.name()}")
  match w.kind {
    @lib.Discrete => {
      println("orthogonal  : \{w.orthogonal}")
      println("biorthogonal: \{w.biorthogonal}")
      println("vm_psi      : \{w.vm_psi}")
      println("vm_phi      : \{w.vm_phi}")
      println("dec_len     : \{w.dec_len()}")
      println("rec_len     : \{w.rec_len()}")
      println("dec_lo      : \{@json.to_json(w.dec_lo).stringify()}")
      println("dec_hi      : \{@json.to_json(w.dec_hi).stringify()}")
      println("rec_lo      : \{@json.to_json(w.rec_lo).stringify()}")
      println("rec_hi      : \{@json.to_json(w.rec_hi).stringify()}")
    }
    @lib.Continuous =>
      println("support     : [\{w.lower_bound}, \{w.upper_bound}]")
  }
}
