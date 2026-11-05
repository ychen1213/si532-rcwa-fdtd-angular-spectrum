# 程序包上传核对记录

核对日期：2026-11-06，Asia/Hong_Kong。

在新增本次讨论原文之前，公共仓库当前提交为 `af07c4f30f889537616d3f6fd38235582350caf8`，仓库工作区干净。原 D4C 程序包已经完整上传，本次无需重传相同的二进制数据。

## 已发布内容

| 项目 | 核对结果 |
| --- | --- |
| 原研究目录 Git 跟踪文件 | 695 个，包括 694 个公开清单条目及清单本身 |
| 原公开清单条目总字节 | 265,461,947 字节 |
| 源目录原有 Python 脚本 | 32 个，全部已跟踪、在清单内、与源文件字节一致 |
| 完成的 FDTD 工况 | 21 个，均有完成记录 |
| FDTD 原生工程 | 21 个 FSP |
| 原始近场数据 | 84 个 NPZ，即每工况 4 个监视面，保留复 E/H |
| RCWA 独立参考 | 9 组 |
| 报告 | Markdown、PDF、HTML、TeX 均已发布 |
| 清单校验 | 原有 694/694 文件 SHA-256 与 Git 索引 blob 字节检查通过 |

核对覆盖 `analysis/prepare_spectra.py`、`analysis/zero_order_models.py`、`analysis/run_analysis.py`、`analysis/physical_component_scan.py` 及全部其他原有 Python 文件；未发现遗漏或未发布的源脚本改动。

科学数组和原生工程按原字节保留。原公开包已声明：22 个求解日志仅移除了机器／许可证抬头信息；原始日志仍保留在本地研究目录。瞬时构建文件、缓存和未完成仿真的临时产物不作为完成结果发布。

## 本次新增

本次增加完整对话 Markdown、问答阅读版 Markdown、原文校验清单、归档导航、此核对记录，以及可复现的原文导出脚本。程序包索引新增归档链接，公开文件清单同步更新；既有科学结果和脚本不改写。

最终公开包的逐文件 SHA-256、文件数量及总字节数，以更新后的 [PUBLIC_COPY_MANIFEST.json](../../PUBLIC_COPY_MANIFEST.json) 为准。提交前再次运行 [check_public_index.py](../check_public_index.py)，核对公开副本及 Git 暂存文件与清单一致。
