# 最终交付核验

`delivery_audit.json` 为最终严格模式核验，overall complete/true，0 FAIL、0 PENDING。实际检查数以该JSON为准；构建期间存在的PDF临时Markdown也可能被链接扫描器计数，因此不把临时文件导致的检查总数变化当作物理验证范围变化。

这些检查包括数据完整性、原始复场及模型哈希、实际边界/网格/源、独立目标重算、共同参数冻结、零级会计恒等式、时间控制和报告链接。检查数量包含逐个分量/文件的完整性检查，不等同于1533项独立物理实验。

实际完成的FDTD引擎求解为21个。另有1个弱结构 finerY 100fs尝试按要求取消，保留CANCELLED记录，无DONE，未作为空间收敛证据。物理精度请看报告中的网格、时间、RCWA级数与Gaussian角依赖对照，而不是只看机械核验数量。

重新运行：

```powershell
python qa/verify_delivery.py --report report/REPORT_zh.md
```

公开副本仅移除运行日志中的机器/许可头部信息，原始scientific NPZ/FSP保留；`PUBLIC_COPY_MANIFEST.json`记录公开文件的精确字节哈希。项目.gitattributes禁用科学数据及元数据的自动换行转换，以保持哈希可复核。
