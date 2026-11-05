"""Build a detailed Chinese textbook from completed, traceable comparisons."""
from pathlib import Path
import argparse,csv,json,re,subprocess,sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1];HERE=ROOT/'report'
def js(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def rows(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def n(x):return f'{float(x):.7g}'
def table(headers,records):
    return '\n| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+''.join('| '+' | '.join(map(str,r))+' |\n' for r in records)+'\n'

def main(allow_partial=False):
    comp=ROOT/'outputs/comparison_fine';summary=js(comp/'summary.json')
    if not summary['all_eight_complete'] and not allow_partial:raise RuntimeError('Eight real FDTD variants must complete before final report')
    extra_path=ROOT/'outputs/extra_weak/summary.json'
    if not extra_path.exists() and not allow_partial:raise RuntimeError('User-requested extra weak-signal holdout must complete')
    if extra_path.exists() and not allow_partial and not js(extra_path).get('controls_complete'):raise RuntimeError('Weak-signal time/mesh controls must complete')
    if extra_path.exists() and not allow_partial and not js(extra_path).get('time_control_passed'):raise RuntimeError('Weak-signal production time must pass its own span criterion')
    er=rows(comp/'candidate_errors.csv');fe=rows(comp/'shared_fit_errors.csv');zp=rows(comp/'zero_power_audit.csv')
    fit=js(comp/'shared_fits.json');req=rows(comp/'required_T0_diagnostic.csv')
    physical=ROOT/'outputs/physical_component_fine'
    pfit=js(physical/'shared_fits.json');per=rows(physical/'shared_fit_errors.csv')
    pmethod='actual_nearest_band_NA0.60'
    gainkey=pmethod+'__nonzero_only__gain_ac'
    awarekey=pmethod+'__nonzero_only_baseline_aware__gain_ac+zero_transmission'
    pgain=pfit[gainkey]['parameters'];paware=pfit[awarekey]['parameters']
    rc=js(ROOT/'rcwa/reference_summary.json');syn=js(ROOT/'analysis/synthetic_control/summary.json')
    ids=summary['available_structures'];lines=[]
    def add(x):lines.append(x.strip()+'\n')
    def figure(path,caption,purpose,discussion):
        add('**图的目的：** '+purpose+'\n\n!['+caption+']('+path+')\n\n**读图与结论：** '+discussion)
    def err(name,model,ch='SUM'):
        return next(r for r in er if r['id']==name and r['model']==model and r['channel']==ch)
    def fitted(name,key,ch='SUM'):
        return next(r for r in fe if r['id']==name and r['fit']==key and r['channel']==ch)
    def recommended(name,ch='SUM',aware=False):
        return next(r for r in per if r['id']==name and r['channel']==ch and r['method']==pmethod and
            r['gain_scope']=='nonzero_only' and r['free']==('gain_ac+zero_transmission' if aware else 'gain_ac') and
            r['objective']==('curve_and_s0_equal_weight' if aware else 'curve_only'))
    bandkey='decband_K_over_gamma_gain+zero_transmission';modekey='dec2_K_over_gamma_gain+zero_transmission'
    bf=fit[bandkey]['parameters'];mf=fit[modekey]['parameters']
    hold=[r for r in fe if r['fit']==modekey and r['split']=='heldout']
    maxhold=max([float(r['relative_L2_pct']) for r in hold],default=float('nan'))
    add('# 实验概览')
    add('532 nm；AH53 与 50% 占空比 Mark；9 周期；38 μm 光斑；SiO₂ 填平覆盖层 0/100/150/200 nm。')
    if extra_path.exists():add('另加用户指定弱信号压力测试：AH53 台阶 180 nm、SiO₂ 覆盖 400 nm，使扫描纵坐标处于 0.00x 量级；它不参与原八组标定。')
    add('本报告按“先生成目标 → 独立求场 → 统一公式 → 控制误差 → 跨结构验证 → 回答应修改哪些参数”的顺序组织。所有扫描曲线保留绝对纵坐标；没有对每条曲线分别平移、归一化峰值或拟合独立增益。报告中的真实 FDTD、RCWA 和人为构造的公式控制分别注明来源。')
    add('# 先明确这次实际要完成的问题')
    add(r'''用户的目标是：先使用已经给定的 D4C 算法，从多个 Mark 的 RCWA 复振幅得到目标扫描曲线；再使用有限 Mark、Gaussian 入射的真实 FDTD 角谱，代入 FDTD 的 B 方法公式；检查是否存在一套共用参数，使所有 Mark 的绝对扫描曲线都一致。特别要检查零级挡板与 $T_0P_0$ 是否被重复计算。

这里有三个不同的核验层次。第一，代码是否忠实计算指定公式；第二，两套电磁求解在相同连续物理条件下是否收敛；第三，给定的有限角谱算符是否等价于 D4C 的离散级次算符。第一层通过不等于后三者自动一致。本实验分别构造对照，不把一个拟合增益解释为所有差异的物理原因。

全流程只使用角谱积分。假设像面无限大，光在像面上的传播积分由 Parseval 关系直接消去；没有再算有限像面窗口，也没有给图像另选接收面积。''')
    add('整体路线为：①匹配材料的 RCWA 生成 8 条 D4C 目标；②全 PML 的有限 9 周期 FDTD 保存复 E/H；③空域 FFT 与共同参考面；④比较级次权重、零级定义、挡板与标量尾项；⑤只用 0/100 nm 训练；⑥锁定参数预测 150/200 nm。')
    add('# 目前最重要的实际结果')
    add(f'完整数据集包含 {len(ids)} 组真实有限 Mark FDTD 对比和 8 组独立 RCWA 目标。补上正确 K/γ、零级按用户给定 T₀ 只计一次之后，固定参数结果见下表。随后只用两种 Mark 的无膜与 100 nm 膜厚拟合一个共用非零级 gain；不修改原 T₀，也不缩放零级项；150/200 nm 为留出数据。')
    add(table(['Mark','膜厚 / nm','固定参数 L2 / %','共用非零 gain 后 L2 / %','s=0 基线差 / %','用途'],[
        [name.split('_')[0],name[-3:],n(err(name,'weighted_nostop_mode2_scalar')['relative_L2_pct']),
         n(recommended(name)['relative_L2_pct']),n(recommended(name)['baseline_relative_pct']),
         '训练' if name.endswith(('000','100')) else '留出'] for name in ids]))
    add(f'这套更简洁、保护零级的参数是 g_ac={n(pgain["gain_ac"])}，T₀ 保持 0.000629728295091705，两种 leakage 保持原值。它只校正当前网格下约 1% 的非零级信号尺度，不给每个 Mark 独立常数。是否需要这约 1% 的修正，要结合后文更细网格控制判断，不能直接视为硬件标定。')
    add('**必须同时核查 s=0 的零级基线。** 整条曲线 L2 主要由较强的非零级调制决定；降低这个误差可能牺牲弱零级基线。所以上表不是“零级已经匹配”的充分证据。图 13 和专门的基线表保留固定 T₀、全曲线拟合 T₀ 的绝对基线，另外检查优先保护零级的拟合。')
    add('# 第一阶段：先用 RCWA 建立 D4C 目标')
    notes=(HERE/'rcwa_reference_notes.md').read_text(encoding='utf8')
    notes=re.sub(r'^#[^\n]+\n','',notes,count=1)
    notes=re.sub(r'^(#{2,4})\s+\d+(?:\.\d+)*\.?\s+',r'\1 ',notes,flags=re.M)
    # The delegated chapter has its own heading hierarchy and nearby figure
    # interpretation. Keep its full scientific content.
    add(notes)
    add('## 关于“图中最小值为什么看起来是 0”的核对')
    add(r'''用户指出：无论 SUM 还是 DIFF，都有漏光项，曲线最小值不应为零。这个判断适用于本组结构。逐个核查原始绘图数组后，16 条目标的最小值全部严格大于零，且分别等于保存的零级残余加 leakage。原全量程线性图中，最小值比峰值小几个数量级，视觉上贴近 0；不能将这种显示效果解释为算法把漏光漏掉或把底值归零。

本次对称结构在 $s=0$ 有 $1-\cos=0$，因此只消掉非零级干涉项，零级尾项仍为

$$S_{\min,\chi}=S_\chi(0)=[T_0+L_\chi(1-T_0)]|r_0|^2>0.$$

为消除误读，下面专门给出谷底局部线性放大图：每个纵轴仍从 0 开始，蓝、橙点与 0 有明确间隔，并标出实际最小值。局部更密位移网格只是对同一 D4C 公式重新求值，没有增加或冒充新的 FDTD 求解。''')
    figure('../figures/16_d4c_positive_minima.png','图 16：SUM 与 DIFF leakage 版本的非零谷底','直接核查漏光项是否保留，并让读者能从图上看见曲线没有到零。',
        '蓝色为 SUM，橙色为文档同一干涉核下采用 DIFF leakage 的版本。每个面板使用适合该结构的局部位移与纵轴范围，横轴为 nm；两种最小值在点旁明确给出。最弱的 AH53 cap150 仍有 SUM=9.78432×10⁻⁶、DIFF=5.16505×10⁻⁶，均不是 0。不能把不同面板的自动纵轴范围当作相同信号尺度，完整绝对量级仍以目标总览和下表为准。')
    minima=rows(ROOT/'figures/data/16_d4c_minima_table.csv')
    add(table(['结构','D4C SUM 最小值','D4C DIFF leakage 最小值'],[[r['id'],n(r['SUM_min']),n(r['DIFF_same_kernel_min'])] for r in minima]))
    add('# 第二阶段：真正的有限 Mark FDTD 如何取得角谱')
    add(r'''## 几何、入射、边界与材料的一致性

每个有限 Mark 沿 x 只有 9 个周期，周期 16 μm，总长度 144 μm；z 方向不变，计算二维 TE，电场分量为 $E_z$。Si 台阶高度 300 nm，下面是 Si 衬底。AH53 的一个周期内有三条宽 1.6 μm 的 Si 线，中心在 $-3.2,0,+3.2$ μm；二元 Mark 的 Si 线宽为 8 μm，即 50% 占空比。

无膜组是裸 Si/空气结构，不能把它误解为已经填平但 cap 厚度为零。其余三组使用 SiO₂ 填满沟槽，并在 Si 台阶最高点上方分别保留 100、150、200 nm；覆盖层顶面平坦。有限 Mark 外围也保留同一平坦覆盖层，从而不引入位于受照区域的额外膜边缘。

532 nm 入射 Gaussian 束的**强度 $1/e^2$ 直径为 38 μm**，因此场为 $E_i(x,0)\propto\exp(-x^2/w_0^2)$，$w_0=19$ μm，腰面 y=0。源在 y=1.35 μm，Mark 不移动；本次复现的是用户指定 B 方法的相位扫描，而非每个位移都移动几何重新求解。

x 计算宽度 192 μm，复场监视器宽 184 μm。y 范围为 −0.8 到 +2.5 μm，四个方向全部是 PML；没有横向周期边界。反射监视器位于 y=1.8、2.0 μm，位于源上方；空域入射参考在 y=0.8 μm。原始 six-component E/H、原生 LSF、模型 FSP、实际 mesh 与结束记录都保留。

生产网格为全域 dx=10 nm、dy=5 nm；对照网格为 dx=20 nm、dy=10 nm。统一网格延伸穿过空气、Si 与边界，避免旧算例中源区粗细网格交界带来的额外回波。材料使用本次 FDTD 的实际拟合复折射率，RCWA 使用完全相同数值；没有分别取材料数据库近似值。''')
    add('Si 的实际拟合折射率为 4.149971915645845+0.048064634885764825i，SiO₂ 为 1.4601075923261915。裸 Si 平面解析反射率为 0.37416823853。网格及有限时间误差都保留在结果中，未乘以 Fresnel 修正因子。')
    add(r'''## 从总复场到同一参考面的反射角谱

采用时间约定 $e^{-i\omega t}$ 与傅里叶变换

$$A(k;y)=\int E_z(x,y)e^{-ikx}\,dx.$$

电场单位 V/m，积分后 $A$ 的单位为 V。反射监视器在源上方，因此主提取使用同一源、同一网格的空结构计算，先作**复电场相减**，再传播到 yᵣ=0.85 μm：

$$A_r(k;y_r)=[A_{\rm mark,R}(k)-A_{\rm air,R}(k)]\exp[i\beta(k)(y_r-y_R)],\quad\beta(k)=\sqrt{k_0^2-k^2}.$$

该相减不是减平坦 Si 背景，也不是对强度相减。减掉平坦 Si 会改变总零级，直接破坏这里要比较的 $I_0$。独立 E/H 方向分离作为交叉核验；连续介质 E/H 导纳与离散 Yee 格的细小差异也单独报告。

空域参考的向下入射角谱使用

$$A_i^{\downarrow}(k)=\tfrac12\left[A_E(k)-\frac{Z_0}{\gamma(k)}A_H(k)\right],\qquad\gamma(k)=\sqrt{1-(k/k_0)^2}.$$

实际入射功率由它计算：

$$P_{\rm in}=\int\frac{\gamma(k)}{4\pi Z_0}|A_i^{\downarrow}(k)|^2\,dk.$$

同一 $P_{\rm in}$ 同时归一化非零级、零级残余与 leakage；没有把峰高当入射功率，也没有重复除以入射量。FFT 零填充只细化 dk，物理接收窗口仍为 184 μm。''')
    add('## 数值误差先独立检查，避免用零级参数补偿')
    figure('../figures/04_numerical_accuracy.png','图 4：平坦 Si 与扫描曲线的网格核验','确定电磁离散误差是否已经足够小，并检查读出参数是否被用来掩盖网格误差。',
        '左图统一扣除空域复场后，平坦 Si 反射率误差从 dx20/dy10 nm 的约 −3.072% 降至 dx10/dy5 nm 的约 −0.750%，再到 dx10/dy2.5 nm 的约 −0.187%。右图使用完全相同的读出公式，展示无膜 AH53 扫描残差随网格改变；因此这部分变化确实来自电磁求解，而非修改 T₀。该控制不能给所有复杂膜层结构提供同一个误差修正量。')
    add('生产仿真时间使用真实短时/长时对照决定。无膜 AH53 长时基准实际计算到约 428.93 fs，达到 9.04×10⁻¹⁰ auto-shutoff；100 fs 对照虽然全域剩余能量未达到该阈值，但 NA=0.6 的扫描曲线最大差仅约目标 span 的 0.0225%。总复场相对 L2 差约 0.0874%，非零角带复场相对差约 0.2338%。因此 100 fs 对这条扫描曲线已足够，但不能把它描述为所有弱衍射级复振幅都达到 0.1%。膜层结构另外做时间控制，具体真实记录见 solver manifest 与 comparison JSON。')
    add('## 150 nm 膜层的时间与更细垂直网格控制')
    control=js(ROOT/'outputs/controls/mesh_time_summary.json')
    figure('../figures/15_film_mesh_and_time.png','图 15：同一膜层结构的时间和网格对照','判别约 1% 的曲线差是否应该由修改零级参数来解决。',
        '三个计算使用同一 K/γ、同一原 T₀ 与 leakage、同一零级分解，没有增益拟合。把时间从 100 fs 延长到 200 fs，仅使最大曲线变化达到 D4C span 的约 0.0270%；把 dy 从 5 nm 降到 2.5 nm，整曲线相对 L2 则从 1.15687% 降到 0.276862%，s=0 基线偏差从 −1.47702% 降到 −0.309295%。这显著支持剩余误差主要来自数值离散，而不是需要把 T₀ 改成另一个很不同的数。')
    add(table(['设置','实际时间 / fs','整曲线相对 L2 / %','s0 基线相对差 / %'],[
        [r['suffix'],n(r['actual_time_fs']),n(r['relative_L2_pct']),n(r['baseline_s0_relative_error_pct'])] for r in control['scan_rows']]))
    add('更细控制使用实际 Custom PML128，标准控制为 Standard PML64；已读回的 kappa/sigma/alpha 及多项式参数相同，y 方向实际 PML 厚度保持 320 nm，x 方向变厚。因此这个控制不应冒称只有一个参数发生改变。空域与双反射监视器核验均保留；平坦 Si 的同步收敛和时间对照共同支持数值解释。初始未扣空场的通量审计给出 fine 约 −0.726%、finerY 约 −0.163%；本报告图 4/15 与扫描主流程统一使用空场复场相减，因小背景的相干项，数值略变为 −0.750%/−0.187%，两种定义都没有通过拟合修正。')
    add('# 第三阶段：将 FDTD 角谱写成与 D4C 可比较的扫描公式')
    add(r'''## B 方法完整积分的对象

令 $G=2\pi/p$。最近级次映射记为 $M(k)$，边界是相邻 $mG$ 的中点。用户指定的 B 方法把同一角带的所有 k 点赋予相同扫描相位 $e^{-iM(k)Gs}$。它不同于真实刚性平移的连续相位 $e^{-iks}$。本次保持用户指定 B 方法，没有偷偷改成另一条扫描定义。

对 $k>0$ 的第 m 个角带 $B_m$，两支的相干积分应写为

$$S_{m}^{\rm EM}(s)=\int_{B_m,k>0}\frac{\gamma(k)}{4\pi Z_0P_{\rm in}}
\left|U(k)e^{-imGs}-U(-k)e^{+imGs}\right|^2dk.$$

可以压缩成两个特征：

$$D_m=\int_{B_m,k>0}w(k)[|U(k)|^2+|U(-k)|^2]dk,$$

$$C_m=\int_{B_m,k>0}w(k)U(k)U(-k)^*dk,\qquad w(k)=\frac{\gamma(k)}{4\pi Z_0P_{\rm in}},$$

$$S_m^{\rm EM}(s)=D_m-2\operatorname{Re}[C_m e^{-2imGs}].$$

这里 $D_m$ 与 $C_m$ 已归一化，单位都是相对于入射功率的无量纲数；$C_m$ 保留相位。只有 $+m/-m$ 等幅同相时，才可把它写成单一 $1-\cos$。本次八组结构镜面对称，RCWA 的配对相位基本为零，故相位正负约定不会影响这些目标曲线。不能由此推断任意非对称 Mark 都可省略相位。

## 为什么必须区分 Kₘ 和 γ(k)

对于非常窄的、彼此分离的第 m 级波包，电磁积分趋向

$$S_m^{\rm EM}\simeq2\gamma_m|r_m|^2[1-\cos(2mGs)],\qquad\gamma_m=\sqrt{1-(m\lambda/p)^2}.$$

用户已确认 D4C 用 $T_m=|r_m|^2$，其目标则是 $2K_m|r_m|^2[1-\cos(2mGs)]$。因此两套强度约定之间的固定转换为

$$\boxed{W_m=K_m/\gamma_m.}$$

最终非零级读出为 $\sum_{m=1}^{18}W_m[D_m-2\operatorname{Re}(C_m e^{-2imGs})]$。这不是逐级自由拟合，而是从两套既定公式的归一化推导的固定权重。连续通量中的 γ(k) 仍保留，不能删掉后又将它等同于 K。若用户选择的是功率效率而非 $|r|^2$，转换式会不同；本次这个歧义已经由用户明确消除。

低级 K 约为 0.49，而 γ 约为 1，因此若直接把纯电磁 B 曲线当 D4C 曲线，会出现近两倍的调制增益。仅拟合一个整体 0.5 可以暂时缓解，却不能严格代替随级次变化的已知 Kₘ，也可能把零级一并错误缩放。''')
    add('# 第四阶段：零级究竟应如何处理')
    add(r'''## 三种零级公式的含义

先定义同一个零级场在挡板之前与之后的功率：

$$P_{0,\rm raw}=\int_{\rm pupil}w(k)|A_0(k)|^2dk,\quad
P_{0,\rm pass}=\int_{\rm pupil}w(k)|T_{\rm stop}(k)A_0(k)|^2dk,$$

$$P_{0,\rm block}=P_{0,\rm raw}-P_{0,\rm pass},\qquad T_{0,\rm eff}=P_{0,\rm pass}/P_{0,\rm raw}.$$

此处 w 已含 $1/P_{\rm in}$，所以三个 P 在程序中是归一化功率。标量模型应只加一次

$$\boxed{Z_{\rm scalar,\chi}=[T_0+L_\chi(1-T_0)]P_{0,\rm raw}.}$$

若真实角谱挡板已经负责计算残余，则应使用

$$\boxed{Z_{\rm physical,\chi}=P_{0,\rm pass}+L_\chi P_{0,\rm block}.}$$

两者在 $T_{0,\rm eff}=T_0$ 时完全等价。相反，以下做法一般重复计入已经通过的那部分：

$$Z_{\rm repeated,\chi}=P_{0,\rm pass}+[T_0+L_\chi(1-T_0)]P_{0,\rm raw}.$$

其与标量分支之差恰好为 $P_{0,\rm pass}$。若固定同一非零级曲线，它只会改变上下位置，不改变扫描 span 或形状。因此某一组的 RMSE 变小，并不能证明重复计数有物理依据：它可能抵消了另一项错误。

原始 D4C 主式明确把干涉求和限定为 $m>0$，再加标量零级。不能先把零级放进理想相干奇偶端口，再机械加一次同一个标量项。SUM/DIFF 在本报告中沿用用户文档给定同一干涉核，仅替换 leakage；这里没有虚构真实 DIFF 端口的光路正负号。

**仅仅在 U 中使用了挡板 mask，并不自动构成重复计数。** 如果干涉求和本来只含 $m>0$，或者对称零级在理想暗口相消，那么已经通过挡板的零级仍可能没有作为 DC 读出，此时标量尾项仍有其算法作用。重复与否必须沿实际读出逐项核对：同一份 $P_{0,\rm pass}$ 是否已经以功率加到结果中。不能仅因公式同时出现 $T_{\rm stop}$ 与 $T_0P_0$ 就断言错误。

## 原挡板半宽并不对应给定 T₀

理想 $w_0=19$ μm Gaussian 的零级角谱功率密度正比于 $\exp[-w_0^2k^2/2]$。不透光、中心对称挡板的半宽为 h 时，忽略极小的 γ 变化有

$$T_{0,\rm eff}\simeq\operatorname{erfc}(w_0h/\sqrt2).$$

用户公式给出 $T_0=0.000629728295091705$，对应 $h=0.1799197989504$ rad/μm；旧半宽 0.2 rad/μm 对应 $T_{0,\rm eff}=0.000144696087850240$，两者相差约 4.35 倍。挡板内部材料的局部功率透射率 τ 也不是整束光的 T₀：即使内部全不透光，外侧 Gaussian 尾部仍通过。

本次直接测量实际 FDTD 源分布的挡板残余，并保留 TE 通量因子。离散硬边界采用边界 cell 与挡板的交叠比例积分，避免 h 在相邻 dk 样点之间跳动所产生的一阶假变化；这是积分处理，不是把物理挡板改成软边。''')
    figure('../figures/05_zero_angular_spectrum.png','图 5：真实 FDTD 的中央角谱与零级估计','看清物理零级波包、最近级次中央角带及挡板边界的关系。',
        '左图来自无膜 AH53 的真实细网格复场。绿色是仅使用空域入射模板估计的零级，橙色为总场减去该估计的余场。右图放大 B0 边界与挡板边界：它们不在同一位置。零级的尾部跨出 B0，邻级的尾部也可进入 B0。因此把中央角带功率直接命名为完整 P₀，是一项近似。模板场是独立场分解估计，并不从 D4C 目标扫描曲线反推。')
    add(r'''## 两种 P₀ 提取都要实际比较

第一种沿用旧 B0 定义，即 $|k|<G/2$ 中的总场功率。这种方法简单，但角带不是独立物理模式。

第二种使用空场实测 $A_i(q)$ 生成移位源模板 $A_i(k-mG)\exp[i\beta(k)y_r]$，联合拟合总复场。每级允许 q 的二阶响应项，从而估计完整 $A_0(k)$；再将 $A_0$ 从总场中减去后构造非零级干涉。这个分解没有输入 RCWA 的 $r_m$，也没有输入目标扫描曲线。保留常数响应与二阶响应的比较，以检查结果是否过分依赖模板自由度。

需要明确：总场 $A=A_0+A_n$ 是相干和，其总功率还含 $2\operatorname{Re}\int w A_0A_n^*dk$。程序记录这项交叉功率，不把分量功率可加当作恒等式。D4C 的标量零级读出本身定义了一个离散级次分离算符；源模板版本是在有限谱上尝试实现这个算符。它不是所有有限结构都存在唯一严格“物理零级”的证明。

对于只输入零级 Gaussian 的控制，最近角带 B 方法仍会把跨越 $G/2$ 的零级尾部赋予 m=1 的扫描相位，产生伪调制。无膜 AH53 对应量级约 $5.802\times10^{-5}$；显式提取并移除完整 $A_0$ 后，纯零级的扫描 range 为零。该反证说明不能仅因中央 B0 不随位移变化，就认为零级所有尾部都不随位移变化。''')
    figure('../analysis/synthetic_control/zero_transmission_definition.png','图 6a：同一挡板对不同零级定义的透过率','用完全相同的 RCWA 复振幅与有限源波包，单独检查 P₀ 定义，不夹带 FDTD 求解误差。',
        '这是 RCWA frozen-q0+Gaussian 合成控制，不能称作有限 Mark FDTD。实测空域源相同、挡板半宽也相同；已知完整 A₀ 时，各结构有效透过率一致。使用 B0 角带后，有效透过率约在 0.0002269 到 0.0016390 之间变化。改变 Mark 就需要不同 T₀ 的表象，可能源于角带混合，而不是仪器真的改变了挡板。')
    figure('../analysis/synthetic_control/absolute_curve_errors.png','图 6b：有限角谱公式控制的绝对误差','判断有限谱与离散级次的差异在不经过 FDTD 时是否已经存在。',
        '每个柱都对照相同 D4C 目标，保留绝对纵坐标。逐物理级积分在 NA=0.6 下的曲线相对差约 0.00585%–0.11641%；第 18 级中心 sinθ=0.5985 靠近 NA=0.6，有限谱被额外截断。固定同一个 K₁₈、将接收扩大至能完整覆盖选中 18 级的 NA=0.65 后，逐级控制误差降至约 0.000153%–0.000237%。这证明边缘截断产生额外作用，但没有证明 K₁₈ 已经包含了哪种真实仪器损失。')
    add('# 第五阶段：八条真实 FDTD 曲线是否都能匹配')
    figure('../figures/06_all_absolute_scans.png','图 7：八组绝对扫描曲线','直接查看真实有限 Gaussian FDTD 和 D4C 的全扫描范围、峰值、基线及波形。',
        '黑线是先计算好的 D4C 目标；橙线保留 B0 零级定义，绿色采用独立源模板零级估计。所有固定参数相同，gain=1，没有逐曲线对齐。AH53 展现多级干涉细节，二元结构更接近低级主导；这两类不同频谱是检验普适参数的必要条件。线条看似重合时，需要图 8 的残差检查绝对误差。')
    figure('../figures/07_scan_residuals.png','图 8：保留直流与形状信息的残差','把肉眼难分辨的差异展开，并检验重复零级项是否只改变曲线的整体上下位置。',
        '纵轴除以每条目标的 span，横轴仍是相同位移。紫线比绿色多加一次已通过的零级功率，因此两者差是一条常量；这不会修复局部峰形或各扫描谐波。如果某一结构因此更接近，也应同时检查其他膜厚是否偏离。表中相对 L2 使用目标整曲线范数，图中 span 归一仅用于显示，不用于重新缩放预测。')
    figure('../figures/08_model_comparison.png','图 9：各零级与权重分支的误差矩阵','把“补 K 权重”“改变 P₀ 定义”“实际挡板”“重复残余项”的效果逐项比较。',
        '每个格给出 gain=1 的绝对曲线相对 L2 百分数，颜色使用对数尺度。原纯电磁公式与 D4C 的权重不同，通常出现近倍数偏差；补上固定 K/γ 后，才有意义继续研究较小的零级基线差。不能在这些步骤之间随意加入独立每条曲线的纵向平移。')
    add(table(['结构','B0+标量 / %','模式+标量 / %','模式+实际挡板一次 / %','模式+重复项 / %'],[
        [name,n(err(name,'weighted_nostop_band_scalar')['relative_L2_pct']),n(err(name,'weighted_nostop_mode2_scalar')['relative_L2_pct']),
         n(err(name,'weighted_matchedstop_mode2_physical')['relative_L2_pct']),n(err(name,'weighted_matchedstop_mode2_double')['relative_L2_pct'])] for name in ids]))
    figure('../figures/09_zero_power_and_transmission.png','图 10：真实 FDTD 中提取的零级功率与透过率','检查剩余基线差究竟来自 P₀ 本身，还是挡板透过率。',
        '左图直接比较 FDTD 提取 P₀ 与 RCWA 的 I₀；即使将 T₀ 设为同一个数，输入零级功率的差也会传递到基线。右图使用同一 h=0.17992 rad/μm，比较通过比例；B0 与源模板的差体现定义影响。膜层改变零级时，不能用固定一个绝对加法常数补偿所有结构，因为正确零级项应随 P₀ 缩放。')
    add('# 第六阶段：只改哪几个共用参数，能达到什么程度')
    add(r'''## 标定与留出设计

训练集固定为 AH53 与 binary50 的无膜、100 nm 组；留出集固定为两种几何的 150、200 nm 组。拟合目标同时包含文档同干涉核下的 SUM/DIFF leakage 版本。每条曲线按其目标均方根做记录平衡，避免强信号曲线独占目标函数，但**预测曲线没有按记录独立缩放**。周期重复端点仅计算一次。

候选共同参数为一个总功率 gain、一个零级 $T_0$；$L_{\rm SUM}$、$L_{\rm DIFF}$ 固定。另一个物理挡板分支只搜索一个共同半宽 h，并只拟合一个共同 gain；其零级项使用实测 passed/blocked，不再拥有独立的标量 T₀。

为什么不同时拟合 T₀ 和两个 leakage？两通道只给出 $c_\chi=T_0+L_\chi(1-T_0)$ 两个组合。若把 T₀、两个 L 与 gain 全部放开，Jacobian 的秩不足；有界最小二乘返回一组数字也不代表这些参数已被辨识。实现会拒绝这种退化拟合。''')
    add(table(['零级定义','共同 gain','共同 T₀','固定 LSUM','固定 LDIFF'],[
        ['B0 角带',n(bf['gain']),n(bf['zero_transmission']),'.00138','.000431'],
        ['源模板 A₀',n(mf['gain']),n(mf['zero_transmission']),'.00138','.000431']]))
    figure('../figures/10_global_fit_validation.png','图 11：共同参数在留出结构上的预测','直接检验只改 T₀ 是否足够，及允许共同 gain 后是否能跨结构改进。',
        '背景有阴影的 150/200 nm 结构没有参与拟合。左图比较绝对曲线相对 L2，右图比较平均上下偏差；二者应一起看。T₀ 只提供与各自 P₀ 成比例的常数，无法改变扫描 span、各谐波比例及局部形状。gain 能改变全曲线尺度，但也无法任意修复各级响应的结构依赖。此图是检验“所有 Mark 共用参数”的主证据。')
    figure('../figures/13_zero_baseline_validation.png','图 13：零位移基线单独核验','防止把整曲线 RMSE 降低误判成用户最关心的零级基线已正确。',
        '绿色保持用户给定 T₀，紫色来自只按全曲线误差优化的共同 gain+T₀，橙色多加了一次已通过的零级。如果紫色基线误差变大，即使其全曲线 L2 更小，也不能把它推荐成正确的零级修正。SUM/DIFF 的泄漏底不同，相同 T₀ 调整造成的相对基线误差也不同；这正是需要保留两种 leakage 控制的原因。')
    add(table(['结构','固定 T₀ 基线差 / %','全曲线拟合基线差 SUM / %','全曲线拟合基线差 DIFF / %'],[
        [name,n(err(name,'weighted_nostop_mode2_scalar')['baseline_s0_relative_error_pct']),
         n(fitted(name,modekey)['baseline_s0_relative_error_pct']),n(fitted(name,modekey,'DIFF')['baseline_s0_relative_error_pct'])] for name in ids]))
    add('这里的基线误差以该结构、该 leakage 版本的 D4C 零位移值为分母。它与周期平均误差不是同一个量：非零级本身也贡献周期平均值。因此用 T₀ 修补周期平均误差时，可能把真正的 s=0 零级底调错。')
    add('## 同时保护整曲线与零级基线的拟合，以及更简单的固定 T₀ 方案')
    add(r'''为避免让强调制淹没零级底，另定义固定的拟合目标：每条训练曲线既提供一个按 RMS 归一化的整曲线残差块，也提供一个按目标零位移基线归一化的 s=0 残差块。两个块的权重均固定为 1，不为各 Mark 单独调权。拟合 $S=g_{\rm ac}N(s)+[T_0+L_\chi(1-T_0)]P_0$，leakage 不变，gain 只作用于整个非零级贡献，包含其自身的周期平均值。''')
    add(f'该基线约束版本得到 g_ac={n(paware["gain_ac"])}、T₀={n(paware["zero_transmission"])}。同时保留更简单的一参数版本：固定用户 T₀，只拟合 g_ac={n(pgain["gain_ac"])}。下表比较两者，不因增加一个参数就预设它更好。')
    add(table(['结构','固定 T₀+g_ac L2 / %','同式 s0 差 / %','保护基线 joint L2 / %','同式 s0 差 SUM / %','同式 s0 差 DIFF / %'],[
        [name,n(recommended(name)['relative_L2_pct']),n(recommended(name)['baseline_relative_pct']),
         n(recommended(name,aware=True)['relative_L2_pct']),n(recommended(name,aware=True)['baseline_relative_pct']),
         n(recommended(name,'DIFF',True)['baseline_relative_pct'])] for name in ids]))
    add('固定 T₀ 版本的零位移偏差直接反映 FDTD 零级功率与 RCWA I₀ 的差，最大约 1.48%，不需要把 T₀ 大幅改到约 0.0003。基线约束拟合只对 T₀ 做小幅调整；其额外自由度并没有自动消除电磁与有限角谱模型的残差。')
    figure('../figures/14_recommended_shared_readout.png','图 14：推荐的简洁共同参数与基线保护对照','直接回答是否存在一套可跨八结构使用的少量参数，并避免以零级正确性换取较小的整曲线误差。',
        f'蓝色采用同一个非零级功率 gain={n(pgain["gain_ac"])}，原 T₀ 与两种 leakage 全部保留。左图显示各结构相对曲线误差；阴影为未参与拟合的膜厚。右图单独显示 s=0 底值，固定 T₀ 时非零级 gain 不会改变这个底值。紫色/橙色是额外保护基线的联合拟合；它未全面优于更简单的固定 T₀ 方案，因此没有必要仅为增加参数而增加参数。')
    add('下面把“假如每个 Mark 只允许改 T₀，需取多少”作为诊断列出。这些数没有用于画预测曲线，也不是逐 Mark 标定方案。若它们不一致，说明只调一个共用 T₀ 不足以修复全部基线；若数值为负或很大，还说明其在补偿零级以外的误差。')
    add(table(['结构','B0 所需 T₀（gain=1）','源模板所需 T₀（gain=1）'],[
        [name,n(next(r['needed_T0_at_gain1'] for r in req if r['id']==name and r['channel']=='SUM' and r['decomposition']=='band' and r['convention']=='K_over_gamma')),
        n(next(r['needed_T0_at_gain1'] for r in req if r['id']==name and r['channel']=='SUM' and r['decomposition']=='2' and r['convention']=='K_over_gamma'))] for name in ids]))
    figure('../figures/11_shared_stop_width_sweep.png','图 12：共同实际挡板半宽扫描','检验实际零级挡板能否由一个共用 h 对应到 D4C 的标量 T₀。',
        '该分支只计一次 passed 零级与 blocked leakage，gain 只在训练结构拟合。虚线为理想 Gaussian 对应给定 T₀ 的 h，而不是预先强行指定的最优点。曲线显示实际 FDTD 与 D4C 目标之间的总模型残差会影响最优 h；不能将通过最小化总曲线误差得到的 h 直接当作真实硬件尺寸。')
    best=summary['best_shared_physical_stop']
    add(f'搜索范围 h=0.14–0.225 rad/μm，步长 0.0025；以训练集为准的共同最优 h={n(best["h_rad_um"])} rad/μm，gain={n(best["gain"])}，训练相对曲线误差 RMS={n(best["train_mean_relative_L2_pct"])}%。该结果与标量共用参数版本一起保留，未针对留出结构重新搜索。')
    if extra_path.exists():
        extra=js(extra_path);eerrors=extra['errors'];et=extra['target_summary'];weakid=extra['id']
        add('# 用户新增的第九组：较厚膜层与 0.00x 弱扫描曲线')
        add('用户明确这里指扫描曲线的纵坐标，不是零级反射率或总反射率。为保持原八组结论可检验，第九组只作额外留出，所有既有参数冻结；它不参与选择新的 gain、T₀ 或 leakage。')
        add('先对原 AH53 的 300 nm 台阶只扫描 SiO₂ cap 厚度：300–1000 nm、间隔 10 nm，共 71 点；这组实际筛选的最小 SUM 扫描峰值约为 0.0881549（cap400），仍没有达到 0.00x。最小零级 I₀ 出现在另一膜厚，不能把它偷换成扫描曲线的峰值。')
        add('因此新增的结构同时采用 180 nm Si 台阶与 400 nm 平坦 SiO₂ cap，顶面位于 y=0.58 μm。周期、线宽、三线位置、9 周期范围、38 μm Gaussian、边界和监视器均与原 AH53 相同。这里明确改变了台阶高度，不能把低扫描幅值单独归因于膜层变厚；其完整输入保存在 inputs/extra_weak_mark.json。')
        figure('../figures/extra_weak_mark_target.png','图 19：新增弱目标与原结构单独增厚膜层的筛选','展示为何需要新增一个不同台阶高度的测试结构，以及低扫描信号与低总反射率的差别。',
            '左图把 D4C 总扫描、非零级贡献与零级常数分开；右图的厚膜扫描始终保持原 300 nm 台阶，星号才是新增 180 nm 台阶。厚膜网格最低峰仍约 0.0882，新结构峰约 0.0018897；新结构总反射率却约 0.3013。因此不能用低扫描峰值冒充低总反射率，也不能把新高度的星号当作同一条单参数厚度曲线的一点。')
        add(table(['目标量','SUM','DIFF leakage 版本'],[
            ['实际扫描最小值',n(et['SUM']['minimum']),n(et['DIFF']['minimum'])],
            ['实际扫描最大值',n(et['SUM']['maximum']),n(et['DIFF']['maximum'])],
            ['扫描 span',n(et['SUM']['span']),n(et['DIFF']['span'])]]))
        figure('../figures/17_extra_weak_scan.png','图 17：第九个弱信号结构的冻结参数预测','在扫描本身只有 0.00x、零级底占比更高时，检验共同参数及零级会计是否仍成立。',
            '黑线是新结构独立 RCWA 的 D4C 目标，不是从 FDTD 拟合所得。绿色保留原 T₀ 并只计一次零级，蓝色使用在原四个训练结构得到的非零级 gain；橙色故意重复已通过零级，紫色使用原全曲线最优化但错误牺牲零级底的有效 T₀。下排单独放大非零谷底。这个弱信号例能把在强曲线上被相对范数淹没的绝对 DC 错误明显显示出来。')
        add(table(['模型（参数全部来自原训练）','SUM L2 / %','SUM 基线差 / %','DIFF L2 / %','DIFF 基线差 / %'],[
            [label]+[n(next(r[key] for r in eerrors if r['model']==model and r['channel']==ch)) for ch,key in
             [('SUM','relative_L2_pct'),('SUM','baseline_s0_relative_error_pct'),('DIFF','relative_L2_pct'),('DIFF','baseline_s0_relative_error_pct')]]
             for model,label in [('fixed_T0_once','固定原 T₀、零级一次'),('B0_band_scalar','原 B0 角带定义、标量一次'),('physical_stop_once','实际挡板、零级一次'),
             ('frozen_nonzero_gain','原共同非零 gain、原 T₀'),('frozen_baseline_aware','原保护基线 joint 参数'),
             ('repeated_passed_zero','故意重复 passed 零级'),('frozen_curve_only_gain_T0','原全曲线最小化 gain+T₀')]]))
        add('这些误差独立列出，没有把新增第九组重新混入原训练集，使其看起来又被“统一调准”。弱信号主结果采用经过额外时间核验的 200 fs；100 fs 的未收敛结果只作为失败时间控制保留，不能用它来宣称共享参数失效或计算已收敛模型误差。')
        add('可以用 span 直接判断“一个 gain 加零级常数能否精确拟合当前所有曲线”。T₀ 与 leakage 只改常数，不改 span；若只允许一个共同非零级 gain，那么每条曲线要求的 gain 都必须等于目标 span/FDTD span。下面只是代数诊断，不拿各自 gain 来生成预测：')
        add(table(['结构','固定原参数 FDTD span / 目标 span','单独匹配 span 所需 gain'],[
            [name,n(err(name,'weighted_nostop_mode2_scalar')['span_ratio']),n(1/float(err(name,'weighted_nostop_mode2_scalar')['span_ratio']))] for name in ids]+
            [[weakid,n(next(r['span_ratio'] for r in eerrors if r['model']=='fixed_T0_once' and r['channel']=='SUM')),
              n(1/float(next(r['span_ratio'] for r in eerrors if r['model']=='fixed_T0_once' and r['channel']=='SUM')))]]))
        add('如果旧结构需要大于 1 的 gain，而新弱结构需要小于 1，就不存在一个共同常数 gain 能把这批当前网格曲线的 span 同时精确对齐；再调零级常数也改变不了这个矛盾。这个判断只针对已声明的读出参数族与当前数据，不将有限的参数测试夸大为所有可能光学模型的不可能性结论。')
        if extra.get('controls_complete'):
            figure('../figures/18_extra_weak_controls.png','图 18：弱信号本身的 100/200/400 fs 时间控制','用弱曲线自己的小 span 评价时间误差，避免沿用强信号的宽松绝对容差。',
                '所有曲线保持原 T₀、leakage 和固定 K/γ，不使用经验 gain。左上、右上分别为两种 leakage 的绝对曲线；下排使用第九组自己的 D4C span 归一化残差。100 fs 对此例未收敛，不能将其偏差解释为共享参数失效；主弱信号对比采用 200 fs，并与独立 400 fs 计算核验。')
            add(table(['网格 / 时间设置','SUM L2 / %','SUM s0 差 / %','DIFF L2 / %','DIFF s0 差 / %'],[
                [suffix]+[n(next(r[key] for r in extra['time_mesh_control_rows'] if r['suffix']==suffix and r['channel']==ch)) for ch,key in
                    [('SUM','relative_L2_pct'),('SUM','baseline_s0_relative_error_pct'),('DIFF','relative_L2_pct'),('DIFF','baseline_s0_relative_error_pct')]]
                    for suffix in ('_fine','_fine_t200','_fine_t400')]))
            add('100 fs/200 fs 的最大差以这条弱曲线自身 span 为分母：'+', '.join(ch+'='+n(value)+'%' for ch,value in extra['weak_time_max_difference_over_own_span_pct'].items())+'。所有完成记录和复场保留，不能将强信号时间控制的百分数直接套用到这里。')
            add('用于验证主 200 fs 结果的 200 fs/400 fs 最大差为：'+', '.join(ch+'='+n(value)+'%' for ch,value in extra['weak_200_vs_400_max_difference_over_own_span_pct'].items())+'。这项比较与100fs失败记录分开报告。')
        weak_notes=(HERE/'extra_weak_mark_notes.md').read_text(encoding='utf8')
        pos=weak_notes.find('## 3.')
        if pos>=0:
            suffix_notes=weak_notes[pos:]
            suffix_notes=re.sub(r'^(#{2,4})\s+\d+(?:\.\d+)*\.?\s+',r'\1 ',suffix_notes,flags=re.M)
            add(suffix_notes)
    add('# 附加对照：对完整物理波包逐级积分')
    pn=(HERE/'physical_component_notes.md').read_text(encoding='utf8')
    pn=re.sub(r'^#[^\n]+\n','',pn,count=1)
    pn=re.sub(r'^(#{2,4})\s+\d+(?:\.\d+)*\.?\s+',r'\1 ',pn,flags=re.M)
    add(pn)
    prec=rows(physical/'candidate_errors.csv')
    add('以下为完整八结构实际运行的未拟合结果，说明逐物理波包积分并不会自动全面改善。保留这些更好或更差的结果，避免只选择有利分支。')
    add(table(['结构','实际场 band L2 / %','重构场 band L2 / %','完整物理 pairs L2 / %'],[
        [name]+[n(next(r['relative_L2_pct'] for r in prec if r['id']==name and r['channel']=='SUM' and r['method']==met)) for met in
        ['actual_nearest_band_NA0.60','projected_nearest_band_NA0.60','projected_physical_pairs_NA0.60']] for name in ids]))
    add('# 为什么仍不应该承诺“任意 Mark 完全一样”')
    add(r'''首先，连续 Maxwell 解的两种照明并不完全相同。目标输入是正入射单一 q=0 的无限周期 RCWA；38 μm Gaussian 具有一段入射角谱，真实场取决于 $r_m(q)$。9 周期已足以把受照边缘推到很远，却不会把 Gaussian 变成平面波。$x=\pm72$ μm 的理想场幅约 $5.8\times10^{-7}$，Mark 外入射功率约 $3.5\times10^{-14}$，说明此例主要差异不是照到端部；固定光斑所产生的角响应平均仍然存在。

其次，FDTD 的空间离散误差仍有限。即使平坦 Si 的反射率误差已小于 1%，它不能保证所有深度、膜厚、相位响应都有同一个乘法误差。将读出 gain 调到某个值，只是一个模型标定，不是把 FDTD 变成精确连续 Maxwell 解。

第三，有限谱与离散级次的读出算符有可量化的差异：B0 并非完整零级；最近级次相位会给零级尾部分配非零扫描频率；NA 边缘会剪切有限第 18 级谱。它们不可能都由一个常量 T₀ 改变形状。前述合成对照已在没有任何 Mark FDTD 求解误差的情况下显示这些作用。

因此，本报告结论限定为：实际测试了给定参数族是否能在本组结构上统一拟合，并给出未被拟合抹掉的残差；没有通过试验有限参数点就宣称所有可能模型都不可能匹配。要将精度继续推进，应先使输入照明与接收算符相同，再做 Maxwell 收敛，之后才辨识硬件的剩余参数。''')
    add('# 最后结论：直接回答用户的问题')
    add(r'''**需要先修改的不是每个 Mark 的独立纵向偏移，而是两套公式的共同定义。** 用户确认 $T_m=|r_m|^2$ 后，非零级 FDTD 通量积分要配合固定的 $K_m/\gamma_m$ 权重；单纯 γ 加权的电磁功率不是 D4C 的经验加权信号。成对复场公式中的两支与因子 2 也必须一致。

**零级只能按一种会计规则计入一次。** 若用标量 D4C 尾项，就用挡板前的完整零级功率并加 $[T_0+L_\chi(1-T_0)]P_0$；若用实际挡板，直接加 $P_{0,\rm pass}+L_\chi(P_{0,\rm raw}-P_{0,\rm pass})$，不再额外加同一个 $T_0P_0$。两者等价条件是实际 passed/raw 等于给定 T₀。38 μm 光斑下，理想硬挡板半宽约 0.17992 rad/μm 对应该 T₀，旧 0.2 rad/μm 本来就不是同一个参数。

**P₀ 的定义也需要统一。** 中央角带并不严格等于完整零级；若目标是复现 D4C 的独立级次处理，应提取完整零级分量并审计其尾部，而不能要求一个重新拟合的 T₀ 吸收所有 Mark 的角带混合。源模板分解是本报告测试的可实施办法，其拟合残差与交叉功率均保留。

**只调零级参数不能修复非零级波形误差。** T₀ 与 leakage 仅增加与 P₀ 成比例的常数，不能改变 span 或扫描谐波。最终是否匹配所有 Mark，必须看锁定参数后的留出曲线与残差表；不能因为某一个 Mark 多加一次零级项后 RMSE 更小，就认定这种做法正确。''')
    maxrec=max(float(recommended(name,ch)['relative_L2_pct']) for name in ids for ch in ('SUM','DIFF'))
    add(f'本次建议保留原 T₀ 和 leakage，使用正确 K/γ 与零级一次规则；若在当前生产网格上需要共同读出校正，只拟合一个非零级 g_ac={n(pgain["gain_ac"])} 即可。它在全部八结构、两种 leakage 版本上的整曲线相对 L2 最大差为 {n(maxrec)}%，零级基线最大相对差约 1.48%。只按整曲线误差把 T₀ 拟到 {n(mf["zero_transmission"])} 虽可让总 L2 更小，却会把零级基线调错，因此不推荐。若要求严格为零的残差，当前数据不支持“再调一个零级挡板参数就能做到”；应继续使网格、入射角平均和接收算符一致。')
    if extra_path.exists():
        ex=js(extra_path)
        before=next(r for r in ex['errors'] if r['channel']=='SUM' and r['model']=='fixed_T0_once')
        after=next(r for r in ex['errors'] if r['channel']=='SUM' and r['model']=='frozen_nonzero_gain')
        add(f'**第九个弱信号例使用完全冻结的原参数，也得到改善。** 在原生产网格、经过时间核验的弱信号记录上，保留原 T₀、零级一次时 SUM 曲线相对 L2 为 {n(before["relative_L2_pct"])}%，套用原八组 gain 后为 {n(after["relative_L2_pct"])}%；其零级基线差为 {n(before["baseline_s0_relative_error_pct"])}%。这支持该共同读出校正在本批数据上的实用性，但不等于各条曲线、各扫描谐波已严格相同；零级常数仍不能消除 span 与形状差。报告保留未收敛100fs、收敛时间核验及入射角控制，不通过重新拟合第九组来改变验证结果。')
    add('# 文件、数据和复现')
    add(table(['文件 / 目录','含义'],[
        ['report/REPORT_zh.md / .pdf / .html','教材全文；HTML 为含公式与内嵌图片的单文件'],
        ['inputs/D4C_FDTD_user_formula_20261105.md','用户提供的原始算法，完整保持原文件'],
        ['docs/user_clarifications.md','Tm 确认、K₂/K₄ 授权插值及膜层选择'],
        ['rcwa/cases/*/reference.npz','独立 RCWA 各级复系数及参考面'],
        ['rcwa/cases/*/d4c_targets.csv','每组结构完整目标曲线'],
        ['fdtd/cases/*/nearfield_*.npz','真实有限结构 six E/H 原始复场'],
        ['fdtd/cases/*/DONE.json','实际网格、时间、终止与求解记录'],
        ['outputs/spectra/*.npz','角谱、实测 Pin、源模板零级分量'],
        ['outputs/comparison_fine/scan_curves.npz','所有固定参数和共同拟合预测曲线'],
        ['outputs/comparison_fine/shared_fits.json','共同参数与训练/留出误差'],
        ['analysis/synthetic_control/','明确标记为合成的有限谱公式控制'],
        ['analysis/prepare_spectra.py / run_analysis.py','从原始场提谱、比较与全局拟合'],
        ['analysis/plot_results.py / report/build_report.py','重建实际图表及本报告'],
        ['qa/','独立科学与交付核验']]))
    add(r'''复现顺序：先按 `experiment.json` 和 `fdtd/run_fdtd.py` 生成实际模型；保存原始场后，逐组运行 `analysis/prepare_spectra.py <variant> --suffix _fine`；再运行 `analysis/run_analysis.py`、`analysis/plot_results.py` 和 `report/build_report.py`。FDTD 重算需要 Lumerical 许可；从保存的 NPZ 重新分析不需要启动电磁求解器。RCWA 求解器、材料与级数收敛保存在 rcwa 目录。

本文的功率公式和会计等式由所列 Fourier 约定直接推导，并以原始场验证。软件含义查阅了 Ansys 的 [FDTD convergence testing](https://optics.ansys.com/hc/en-us/articles/360034915833-Convergence-testing-process-for-FDTD-simulations)、[finite array / assume periodic distinction](https://optics.ansys.com/hc/en-us/articles/360034914793-Understanding-the-Assume-Periodic-option-in-far-field-projections)、[grating projections](https://optics.ansys.com/hc/en-us/articles/360034394354-Grating-projections-in-FDTD-overview) 与 [numerical dispersion and mesh](https://optics.ansys.com/hc/en-us/articles/360034382634-Understanding-the-non-uniform-mesh-in-FDTD)。这些官方资料支持方法边界；本报告的数值结论来自本次真实文件，不能用文献替代。''')
    add('# 补充讨论与用户原文确认')
    add(re.sub(r'^#[^\n]+\n','',(ROOT/'docs/user_clarifications.md').read_text(encoding='utf8'),count=1))
    add('原始完整算法文件随交付保存为 [D4C_FDTD_user_formula_20261105.md](../inputs/D4C_FDTD_user_formula_20261105.md)。这里不把原文尚未确认的事项静默改写；Tm 与 K₂/K₄ 的新确认在上一节独立追加。')
    md='\n'.join(lines)
    # Portable MathML and XeLaTeX glyphs; mathematical content is unchanged.
    md=re.sub(r'\\rm\s+([A-Za-z]+)',lambda m:chr(92)+'mathrm{'+m.group(1)+'}',md)
    md=md.replace('𝜒','χ')
    for i,glyph in enumerate('①②③④⑤⑥',1):md=md.replace(glyph,f'({i})')
    (HERE/'REPORT_zh.md').write_text(md,encoding='utf8')
    metadata='''title: "D4C 扫描曲线与有限 Gaussian FDTD 的统一参数对标"
subtitle: "8 组主验证 + 厚膜弱信号压力测试 · 完整角谱积分 · 零级会计"
author: "可复现数值实验报告"
date: "2026-11-05（Asia/Hong_Kong）"
lang: zh-CN
documentclass: ctexart
classoption: [UTF8, a4paper, 11pt]
geometry: [margin=21mm]
mainfont: "Times New Roman"
monofont: "Consolas"
toc: true
toc-depth: 2
numbersections: true
colorlinks: true
linkcolor: blue
urlcolor: blue
header-includes:
  - \\usepackage{float,needspace,etoolbox}
  - \\floatplacement{figure}{H}
  - \\BeforeBeginEnvironment{longtable}{\\Needspace{9\\baselineskip}}
  - \\setlength{\\emergencystretch}{3em}
  - \\setlength{\\parskip}{0.35em}
'''
    (HERE/'metadata.yaml').write_text(metadata,encoding='utf8')
    (HERE/'style.css').write_text('body{max-width:1050px;margin:3em auto;padding:0 2em;font:17px/1.8 "Microsoft YaHei",sans-serif;color:#243346}h1{border-bottom:2px solid #2b7773;padding-bottom:.35em;margin-top:2em}h2{color:#216b69}img{max-width:100%;height:auto}table{border-collapse:collapse;width:100%;font-size:.88em}th,td{border:1px solid #cad5db;padding:8px}th{background:#e9f1f3}pre{overflow:auto;background:#f1f4f6;padding:1em}math[display=block]{overflow:auto;margin:1.5em 0}a{color:#176a98}blockquote{border-left:4px solid #8aaba9;padding-left:1em}',encoding='utf8')
    try:
        import pypandoc
    except ModuleNotFoundError:
        # Optional local cache; ordinary reproduction uses requirements.txt.
        sys.path.insert(0,str(ROOT.parent/'Si532_RCWA_FDTD_AngularSpectrum_20261031/.builddeps'))
        import pypandoc
    pandoc=pypandoc.get_pandoc_path()
    args=[pandoc,'REPORT_zh.md','-f','markdown+tex_math_dollars+tex_math_single_backslash','-s','--metadata-file=metadata.yaml','--resource-path=.;..']
    from PIL import Image
    def reserve_figure(match):
        block=match.group(0)
        file=re.search(r'!\[[^\]]*\]\(([^)]+)\)',block)
        w,h=Image.open(HERE/file.group(1)).size
        reserve=min(38,max(20,int(11+30*h/w)))
        return chr(92)+f'Needspace{{{reserve}'+chr(92)+'baselineskip}\n\n'+block
    pdf_md=re.sub(r'\*\*图的目的：\*\*[^\n]+\n\n!\[[^\]]*\]\([^)]+\)',reserve_figure,md)
    pdf_md=pdf_md.replace('# 补充讨论与用户原文确认',chr(92)+'clearpage\n\n# 补充讨论与用户原文确认\n\n'+chr(92)+'begingroup'+chr(92)+'small\n')
    pdf_md+='\n'+chr(92)+'endgroup\n'
    pdf_source=HERE/'.report_pdf_build.md';pdf_source.write_text(pdf_md,encoding='utf8')
    pdfargs=[pandoc,pdf_source.name,*args[2:]]
    try:subprocess.run(pdfargs+['--pdf-engine=xelatex','-o','REPORT_zh.pdf'],cwd=HERE,check=True)
    finally:pdf_source.unlink(missing_ok=True)
    subprocess.run(args+['--toc','--mathml','--embed-resources','--css=style.css','-o','REPORT_zh.html'],cwd=HERE,check=True)
    subprocess.run(args+['-o','REPORT_zh.tex'],cwd=HERE,check=True)
    print(json.dumps(dict(characters=len(md),real_structures=len(ids),report=str(HERE/'REPORT_zh.pdf')),ensure_ascii=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--allow-partial',action='store_true');a=p.parse_args();main(a.allow_partial)
