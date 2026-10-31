# RCWA 基准与复振幅相位审计

本目录包含真实 Ansys Lumerical v252 RCWA 计算结果、原生模型、原始数据、收敛试验和可重现脚本。**与有限长 FDTD 的偶数周期、中心沟槽结构比较时，请读取 `reference_for_analysis.npz` 中的 `m` 和 `r_native_origin`。** 这里的 `native_origin` 指 FDTD 原生坐标的横向原点，并不指未经处理的 RCWA 原生输出。

## 1. 物理结构与材料

采用时间约定 $\exp(-i\omega t)$。入射沿 canonical $+z$，TE 电场为 $E_y$；空气位于 $z<0$，Si 脊顶位于 $z=0$，沟槽底与基底顶位于 $z=100\,\mathrm{nm}$。周期 $\Lambda=2\,\mu\mathrm m$，Si 脊宽 $1\,\mu\mathrm m$，占空比 50%，波长 $532\,\mathrm{nm}$。RCWA 单胞为 $x\in[-1,+1]\,\mu\mathrm m$，脊位于 $[-0.5,+0.5]\,\mu\mathrm m$。

材料来自 `Si (Silicon) - Palik`，实际查询结果保存在 `material.json`：

| 查询 | 532 nm 的复折射率 |
|---|---|
| 数据库 `getindex` | $4.1501040188082685+0.04393515056151286i$ |
| FDTD 500–565 nm 拟合的 `getfdtdindex` | $4.149873311036015+0.04833494047986202i$ |

为使单频 RCWA 与 FDTD 使用相同材料，RCWA 将第二行的复折射率写入固定 `(n,k) Material`，并关闭再次拟合。`probe_material.py` 中另外通过 `getindex` 回读固定材料并检查数值一致。这避免把数据库插值与 FDTD 拟合的微小差别误认为角谱提取误差。

## 2. 求解器输出与电场系数

反射电场展开使用

$$E_y^r(x,z)=\sum_m r_m\exp(ik_{xm}x-ik_{zm}z),\qquad k_{xm}=\frac{2\pi m}{\Lambda},\quad k_{zm}=\sqrt{k_0^2-k_{xm}^2}.$$

参考面为脊顶 $z=0$，并令该平面的入射 $E_y=1$。空气中只有 $m=-3,\ldots,+3$ 为传播级。相应功率为 $R_m=\operatorname{Re}(k_{zm})|r_m|^2/k_0$。

应使用 `field_amplitudes.eky_bs` 提取电场幅值。`grating_characterization.Rss` 是功率归一化的散射系数，其模平方直接对应衍射效率；它与上述 $E_y$ 系数的归一化不同。有关原生结果的定义见 [Ansys RCWA solver 文档](https://optics.ansys.com/hc/en-us/articles/12959229278611-RCWA-Solver-Simulation-Object)。

本模型的第一界面为 $z_{\rm ref}=-200\,\mathrm{nm}$。沿 z 的去嵌入以及 TE 入射符号校正为

$$r_m^{\rm raw,z0}=-\mathrm{eky\_bs}_m\exp[i(k_0+k_{zm})z_{\rm ref}].$$

负号来自该正入射原生 s 基底的入射 $E_y=-1$。`flat_calibration` 使用真正的平坦 Si 界面验证此式：得到 Fresnel 电场反射系数 $(1-n)/(1+n)$，最大复数误差 $1.51\times10^{-14}$。这一步只能校准电场正负号与 z 参考面；零级对横向原点不敏感，所以仍需下面的独立空间审计。

## 3. 两个必须分开的横向相位问题

### 3.1 RCWA 原生横向参考从单胞左边界开始

原生 `field_amplitudes` 在本次求解中的横向相位基准是单胞左边界 $x_{\min}=-1\,\mu\mathrm m$。这不是根据 FDTD 结果调符号，而是由另外运行的 RCWA 实空间场监视器独立确定。

`audit_spatial_origin.py` 在空气中 $z=-100\,\mathrm{nm}$，对物理坐标 $x\in[-1,+1]\,\mu\mathrm m$ 取 1601 个场点；直接进行 $E_y(x)\exp(-ik_{xm}x)$ 的复数梯形积分，再做 z 去嵌入。结果满足

$$r_m^{\rm spatial}=r_m^{\rm raw,z0}\exp(-ik_{xm}x_{\min})=(-1)^m r_m^{\rm raw,z0},$$

最大复数误差为 $6.41\times10^{-15}$。完整原始场与积分结果分别保存在 `spatial_origin_Es.npz`、`spatial_origin_projection.npz` 和 `spatial_origin_projection.json`。场监视器的原生定义见 [Ansys RCWA field monitor 文档](https://optics.ansys.com/hc/en-us/articles/16277147107731-RCWA-Field-Monitor-Simulation-Object)。

### 3.2 有限 RCWA 网格留下精确的半个实际网格线性相位

样品几何和导出的实际折射率采样均以 $x=0$ 镜像对称，Si 采样坐标均值在 $10^{-13}\,\mathrm{nm}$ 量级。但未校正的原生场系数满足

$$\frac{r_m^{\rm raw,z0}}{r_{-m}^{\rm raw,z0}}=\exp(i mG\Delta x_{\rm actual}),\qquad G=2\pi/\Lambda.$$

此关系在六组截断/网格设置、$m=1,2,3$ 中的最大复数残差为 $2.05\times10^{-13}$，表明它是线性横向相位偏置。重参照到镜像中心所需的因子为 $\exp(-ik_{xm}\Delta x_{\rm actual}/2)$。校正后 $r_m=r_{-m}$ 的最大残差为 $1.47\times10^{-14}$。

这项校正不改变任何 $R_m$，也不是声称几何中心发生了位移。它是本版本、二维 TE 求解设置中用精确镜像对称性校准的有限网格相位偏置。`sampled_index.npz`、`actual_mesh.json` 和 `phase_origin_audit.csv` 保留全部证据。更细网格试验还验证了未校正偏置随实际网格缩小，校正后的场系数稳定。

### 3.3 与本项目 FDTD 几何的对应

FDTD 有限偶数周期结构采用 $x_j=(j-(N-1)/2)\Lambda$，因此 $x=0$ 为沟槽中心；它相对 RCWA 的物理结构平移了半个周期。这个物理平移也乘以 $\exp(-im\pi)$。它恰好抵消从 RCWA 左边界转到全局坐标的 $\exp(+im\pi)$。

因此，最终与 FDTD 同原点的系数为

$$\boxed{r_m^{\rm FDTD\ origin}=-\mathrm{eky\_bs}_m\exp[i(k_0+k_{zm})z_{\rm ref}]\exp(-ik_{xm}\Delta x_{\rm actual}/2).}$$

若要输出本目录模型中脊中心为 $x=0$ 的物理系数，则再乘 $(-1)^m$。这两种系数都保存在最终 NPZ 和 CSV 中。**不要在 `r_native_origin` 上再额外乘一次 $(-1)^m$。**

## 4. 收敛与最终结果

计算按最大保留级数 20、40、80，随后减小实际网格，再增加到 160 级确认。v252 二维模型的 `max number k vectors=H` 实际保留 $2H+1$ 个谐波。申请 5 nm 网格的 80 级模型被求解器自动细化成 640 个单胞网格、实际 $3.125\,\mathrm{nm}$；文中的相位校正始终使用实际网格。

| 设置变化 | 校正后最大 $|\Delta r_m|$ | 最大 $|\Delta R_m|$ |
|---|---:|---:|
| 20 → 40，实际 dx=5 nm | $9.737\times10^{-4}$ | $1.751\times10^{-4}$ |
| 40/dx5 → 80/dx3.125 nm | $1.192\times10^{-4}$ | $1.379\times10^{-5}$ |
| 80/dx3.125 → 80/dx2.5 nm | $1.029\times10^{-5}$ | $5.304\times10^{-6}$ |
| 80/dx2.5 → 80/dx1.25 nm | $1.371\times10^{-5}$ | $7.065\times10^{-6}$ |
| 80 → 160，实际 dx=1.25 nm | $1.457\times10^{-5}$ | $1.936\times10^{-6}$ |

最终参考使用 321 个谐波、实际 $1.25\,\mathrm{nm}$ 网格。全部逐级数据、实际运行耗时和求解器计时见各子目录 `summary.json`。最终反射总功率为 $0.36136091512762$，逐级功率之和与求解器原生 `total_energy.Rs` 一致到浮点精度。

| 级次 | 角度 | FDTD 沟槽原点的电场复系数 | $R_m$ |
|---|---:|---|---:|
| 0 | 0° | $-0.141144589-0.236559705i$ | 0.0758822890 |
| ±1 | ±15.426381° | $0.333462073-0.148842777i$ | 0.1285468895，各一级 |
| ±2 | ±32.140687° | $0.063222182+0.013902463i$ | 0.0035481321，各一级 |
| ±3 | ±52.939539° | $-0.092036922+0.095872119i$ | 0.0106442914，各一级 |

## 5. 重现顺序和文件约定

一条命令即可从没有任何求解结果的干净源码目录重建全部六组光栅、平坦校准、网格审计、最终适配器和独立空间相位审计。在具备 v252 许可的 Windows 机器上，先按实际安装位置设置 `LUMERICAL_ROOT`，在仓库根目录用 PowerShell 运行：

```powershell
$env:LUMERICAL_ROOT = 'D:\Ansys Inc\v252\Lumerical'
& "$env:LUMERICAL_ROOT\python\python.exe" .\rcwa\reproduce_reference.py
```

`reproduce_reference.py` 默认包含所有六组光栅：20/dx5、40/dx5、80/dx5、80/dx2.5、80/dx1.25 和 160/dx1.25 nm；不是只运行前三档截断。若从包含已归档结果的公共仓库 clone，完整求解器输出会被复用，只重新生成与验证适配器，不需要重新运行 Maxwell 求解。可添加 `--dry-run` 查看所需步骤而不执行。

若需要逐步运行，等价的完整顺序为（使用上面同一个 Python 可执行文件）：

```text
python probe_material.py
python run_rcwa.py
python inspect_mesh.py
python finalize_reference.py
python audit_spatial_origin.py
```

每个真实计算子目录保留 `model.fsp`、独立可运行的 `build_run_export.lsf`、`native_datasets.mat`、各矩阵数据集的 NPZ/JSON、原生级次 CSV、模型输入和运行统计。原生 `index` 数据因 v252 Python 翻译器对界面维度处理异常，另用 Lumerical 脚本逐个取出原始坐标与 `index_z` 数组保存；完整原生对象仍在 MATLAB 文件中。所有原始求解器输出保留，不由校正结果覆盖。

最终适配文件 `reference_for_analysis.npz` 的关键字段：

- `m`：$-3,\ldots,+3$。
- `r_native_origin`：与本项目 FDTD 沟槽中心原点一致的最终电场复系数。
- `r_center_origin`：RCWA 几何脊中心为物理 $x=0$ 时的最终电场复系数。
- `r_rcwa_cell_origin`：去除半网格相位后的原生单胞参考系数；在本几何中数值等于 `r_native_origin`。
- `native_eky_bs`：原始 RCWA 电场幅值；`r_before_lateral_correction`：仅校正 TE 正负号与 z 参考面的系数。
- `reference_z_factor`、`lateral_half_mesh_factor`、`rcwa_xmin_to_physical_origin_factor`、`fdtd_groove_origin_factor`：每一步复数因子，便于逐项复算。
- `R`、`kx`、`kz`：对应功率与波矢。

最终完整数值表是 `final_reference.csv`；相位审计、收敛和参考定义分别见 `phase_origin_audit.csv`、`convergence_corrected.csv`、`reference_metadata.json`。
