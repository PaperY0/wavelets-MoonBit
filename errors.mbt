// 库级错误类型：所有公开 API 的失败路径都收口到 `WaveletError`。

///|
/// 小波库错误。
#derive(Debug, Show)
pub(all) suberror WaveletError {
  /// 未知小波名（不在内置注册表中）
  UnknownWavelet(String)
  /// 该操作需要离散小波，但传入了连续小波
  NotDiscrete(String)
  /// 该操作需要连续小波，但传入了离散小波
  NotContinuous(String)
  /// 未知扩展模式名
  InvalidMode(String)
  /// 参数非法（空信号、负 pad 宽度、尺度非正等）
  InvalidArgument(String)
  /// 分解层数非法（负数或超过最大有用层数）
  InvalidLevel(String)
  /// 系数/矩阵形状不匹配
  ShapeMismatch(String)
  /// v0.1 尚未支持的功能（见 docs/known-deviations.md）
  NotSupported(String)
}
