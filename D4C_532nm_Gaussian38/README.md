# D4C 扫描曲线与有限 Gaussian FDTD 的统一参数对标

先用匹配材料的 RCWA 生成 D4C 目标，再用真实有限 9 周期 Gaussian FDTD 的完整复角谱积分计算扫描曲线；检查共同参数能否同时覆盖 AH53、50% 占空比 Mark，以及 0/100/150/200 nm SiO₂ 填平覆盖层。

- [教材报告 PDF](report/REPORT_zh.pdf)
- [教材报告 Markdown](report/REPORT_zh.md)
- [单文件离线 HTML](report/REPORT_zh.html)
- [用户原始公式](inputs/D4C_FDTD_user_formula_20261105.md)及[后续确认](docs/user_clarifications.md)
- [完整实验定义](experiment.json)
- [全部对话原文 Markdown（2026-10-31 至 2026-11-06）](docs/discussions/README.md)：含全部问答、进度消息、逐条原文校验及程序包核对；另提供不含进度消息的问答阅读版。

532 nm、二维 TE、16 μm 周期、300 nm Si 台阶、38 μm 强度 1/e² 直径光斑（w₀=19 μm）；四边 PML，真实 finite Mark，不使用横向周期边界。主比较采用用户已确认 Tₘ=|rₘ|²，K₂=.4947、K₄=.49135 为用户授权线性插值。SUM/DIFF 只比较原文同一干涉核下的 leakage 参数，不推定未知真实端口符号。

全部扫描在无限像面假设下使用角谱功率积分；没有再加像面接收面积。比较保留绝对值，不逐条曲线平移或归一化峰高。训练集固定为两种几何的 0/100 nm 膜厚，150/200 nm 锁定参数后验证。

另加用户指定的第九个弱扫描例：AH53 台阶 180 nm、SiO₂ 平坦覆盖 400 nm；D4C 的 SUM 峰值约 0.0018897、非零谷底约 0.000593687。它改变了台阶高度，不能归因于单纯膜厚变化；总反射率约 0.301，不能把扫描信号称作总反射率。该例不参与原参数训练，另有100/200/400 fs 时间核验。专用数据在 `outputs/extra_weak/`。

所有 D4C 谷底都保留零级残余和 leakage，均大于零。报告增加了独立谷底放大图与逐结构最小值表，避免全量程线性图使微小底值看起来等于零。

## 数据

`rcwa/cases/*/reference.npz` 为独立复振幅参考，`d4c_targets.csv` 为全扫描曲线。`fdtd/cases/*` 保留实际 FSP、原生 LSF、six-component 复 E/H、网格与 DONE。`outputs/spectra/` 为归一化前的共同参考面角谱，`outputs/comparison_fine/` 为固定参数和共同拟合结果。

`analysis/synthetic_control/` 明确为 RCWA frozen-q0 与有限 Gaussian 的公式控制，不是 Mark FDTD。`outputs/physical_component_fine/` 是用真实 FDTD 源模板系数构造的物理波包读出对照，不是额外电磁仿真。

## 复现

电磁求解需要自备 Lumerical v252 许可，安装路径在 `fdtd/runtime_support.py` 中配置。网格与实际结束时间以各 DONE 为准，不用请求参数代替实际参数。保留的复场足以在不启动商业求解器的情况下重做全部后处理。

```powershell
python -m pip install -r requirements.txt
python analysis/prepare_spectra.py ah53_cap000 --suffix _fine
# 对 experiment.json 中另外七个 variant 同样执行上行
python analysis/run_analysis.py
python analysis/physical_component_scan.py --suffix _fine --require-all-eight
python analysis/plot_results.py
python report/build_report.py
python qa/verify_delivery.py --report report/REPORT_zh.md
```

PDF 构建还需要 XeLaTeX（ctex 与中文字体），HTML/PDF 构建需要 Pandoc（pypandoc_binary 提供）。脚本只从保存的结果生成图表，不以解析曲线冒充 FDTD。
