"""Additional held-out weak-scan target and the film-only screening evidence."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from te_planarized_rcwa import dump
from run_reference import csvwrite

ROOT=Path(__file__).resolve().parents[1];R=ROOT/'rcwa';F=ROOT/'figures'


def main():
    c=json.loads((R/'thick_film_candidate.json').read_text(encoding='utf8'))
    sweep=json.loads((R/'thick_film_sweep.json').read_text(encoding='utf8'))
    stats=c['stats'];rows=sweep['all_sweep_rows']
    with np.load(ROOT/c['d4c_targets']) as z:
        s=z['s_um'];y=z['S_SUM'];ac=z['AC'];zero=float(z['zero_SUM']);d=z['S_DIFF_same_kernel']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
    fig,aa=plt.subplots(1,2,figsize=(9.7,4.7),layout='constrained')
    aa[0].plot(s,y,color='#235d8c',label='Extra mark: complete D4C SUM',lw=1.7)
    aa[0].plot(s,ac,color='#bd6b26',ls='--',label='Nonzero interference contribution',lw=1)
    aa[0].axhline(zero,color='#7e4167',ls=':',label='Scalar zero term counted once')
    aa[0].set(xlim=(-4,4),xlabel='s / um',ylabel='D4C model signal',
        title=f"(a) Extra held-out target: h180 / cap400 nm; N{stats['max_order']}")
    aa[0].ticklabel_format(axis='y',style='plain',useOffset=False)
    aa[1].plot([r['cap_nm'] for r in rows],[r['S_SUM_peak'] for r in rows],color='#235d8c',label='Original h300 nm; cap-only sweep N160')
    aa[1].axhspan(.001,.009,color='#b4d8ba',alpha=.5,label='Requested scan peak: 0.001 to 0.009')
    aa[1].plot([400],[stats['S_SUM_peak']],'*',color='#8b425a',ms=10,label='New height h180 nm (not film-only)')
    minimum=sweep['min_scan_peak']
    aa[1].annotate(f"h300 grid minimum: {minimum['S_SUM_peak']:.5f}",xy=(minimum['cap_nm'],minimum['S_SUM_peak']),
        xytext=(465,.067),arrowprops=dict(arrowstyle='->',color='#485968'),fontsize=9)
    aa[1].set(xlim=(300,1000),ylim=(0,.13),xlabel='Planar SiO2 cap / nm',ylabel='Maximum D4C SUM scan signal',
        title='(b) Film-only scan does not reach the requested scale')
    for ax in aa:ax.grid(alpha=.2)
    hh=[];ll=[]
    for ax in aa:
        h,l=ax.get_legend_handles_labels();hh.extend(h);ll.extend(l)
    fig.legend(hh,ll,loc='outside lower center',ncol=2,fontsize=8.5)
    fig.suptitle('Weak scan amplitude is distinct from low zero-order or low total reflection\n'+
        f"Extra target peak = {stats['S_SUM_peak']:.8f}; I0 = {stats['I0']:.5f}; total R = {stats['R']:.5f}",fontsize=11.5)
    for ext in ('.png','.svg'):fig.savefig(F/('extra_weak_mark_target'+ext),dpi=210,bbox_inches='tight')
    plt.close(fig)
    csvwrite(F/'extra_weak_mark_target.csv',[dict(s_um=float(a),S_SUM=float(b),AC=float(q),zero_SUM=zero,S_DIFF_same_kernel=float(dd)) for a,b,q,dd in zip(s,y,ac,d)])
    dump(F/'extra_weak_mark_target.json',dict(candidate=c,original_h300_sweep_source='rcwa/thick_film_sweep.csv',
        figure_scope='RCWA D4C target and cap-only screen; no FDTD result or fitted extra-mark parameters in this figure.'))
    conv='\n'.join('| '+str(row['max_order'])+' | '+f"{row['S_SUM_peak']:.10f}"+' | '+f"{row['I0']:.9f}"+' | '+f"{row['R']:.9f}"+' | '+
        ('—' if 'relative_complex_L2' not in row else f"{100*row['relative_complex_L2']:.6f}%")+' | '+
        ('—' if 'nonzero_relative_complex_L2' not in row else f"{100*row['nonzero_relative_complex_L2']:.6f}%")+' |' for row in c['convergence'])
    text=r'''# 额外留出：扫描纵坐标为 0.00x 的弱信号 Mark

本次新增要求指的是**扫描曲线纵坐标的量级**，而不是零级 $|r_0|^2$ 或所有反射级相加的总反射率 $R$。因此候选筛选直接使用完整 D4C $S_{\rm SUM}(s)$ 的扫描最大值，要求其处于 0.001–0.009。$I_0$ 与 $R$ 仍同时保存，用于防止把三个量混为一谈。

## 1. 先检查只增厚原结构的膜是否足够

保持原 AH53 的 300 nm Si 台阶、16 µm 母周期、1.6 µm 凸台宽和三个局部中心不变，令平坦化 SiO2 cap 从 300 扫到 1000 nm，间隔 10 nm，共 71 点。每点都计算完整各级复振幅并代入同一 D4C 公式。N160 的图形层输入导纳在厚度间复用，均匀 cap 用散射递推接入，不重新计算相同本征层。递推与原完整模态求解在已有 cap100/150/200 三点的最大逐复数差均小于 $10^{-10}$。

这份厚度网格上最低扫描峰是 **@MINPEAK@**，位于 cap = @MINCAP@ nm。未得到 0.001–0.009。与此同时，最低零级 $I_0$ 是 **@MINI0@**（cap@MINI0CAP@ nm），但其扫描峰为 **@MINI0PEAK@**；最低总反射率为 **@MINR@**（cap@MINRCAP@ nm）。因此，挑一个低 $I_0$ 的厚膜点不能回答“扫描幅值 0.00x”的要求。

这些结论针对已计算的 10 nm 厚度网格；不是对连续厚度或所有可能材料层的全局最小值证明。原始厚度网格、复数递推校验和实际扫描指标在 `rcwa/thick_film_sweep.csv/.json` 中。

## 2. 额外改变台阶高度，并将它作为第九个留出结构

为获得所需弱扫描信号，额外结构采用：

| 参数 | 新结构 ah53_h180_cap400 |
|---|---|
| 几何 | AH53，九个 16 µm 母周期，凸台宽 1.6 µm，局部中心 −3.2/0/+3.2 µm |
| Si 台阶高度 | **180 nm**，与原八组的 300 nm 不同 |
| SiO2 | 平坦化填槽，并在凸台上覆盖 **400 nm cap** |
| 最上表面 | $y=0.58\,\mu\mathrm m$，兼容已有 $y_{\rm ref}=0.85\,\mu\mathrm m$ 与源面 1.35 µm |
| 材料 | 与原八组完全相同的 FDTD 532 nm 实际拟合值 |
| 数据用途 | 第九个额外留出；不进入原训练集，不重新拟合 gain、$T_0$ 或结构专用偏置 |

这里既改变了台阶高度，也改变了 cap，不能把获得弱扫描信号归因于单纯膜层增厚。180 nm 接近 SiO2 中 532 nm 光的半波厚度（约 182 nm），提供了可降低局部凸台/沟槽反射相位差的候选思路；实际结果以完整电磁求解为准，侧壁散射与多次耦合没有被局部相位模型代替。

![额外弱扫描目标，以及原高度不变的厚膜筛选](../figures/extra_weak_mark_target.png)

**本图目的：** 直接显示新增结构的原始 D4C 目标幅值，同时展示为什么原台阶只改 cap 的已扫网格没有达到目标。左图分开完整扫描、非零级干涉贡献和零级常数；右图蓝线始终是 h300 nm，星号明确属于改变到 h180 nm 的新结构。

**本图结论：** 最终 N@N@ 的新增扫描峰为 **@PEAK@**、零位移底值为 **@ZERO@**，符合所需 0.00x 扫描量级。其 $I_0=\mathbf{@I0@}$、总 $R=\mathbf{@R@}$，均不在 0.00x 量级。左图来自 RCWA 与固定 D4C 系数，没有用 FDTD 输出或第九结构拟合参数；是否被 FDTD 预测正确仍要看冻结参数的独立留出对比。

## 3. 弱衍射参考须单独控制截断

| 最大级次 N | SUM 扫描峰 | $I_0$ | 总 R | 最后完整复向量变化 | 最后非零级复向量变化 |
|---:|---:|---:|---:|---:|---:|
@CONV@

弱衍射时，总反射复向量可被较强零级主导，因此同时报告非零级子向量的变化，不能只看包含零级的整体百分比。各档真实 RCWA 数据、完整目标数组及与原八组相同的字段位于 `rcwa/cases/ah53_h180_cap400/`。`inputs/extra_weak_mark.json` 独立记录第九结构，原 `experiment.json` 的八个 variants 与原训练/留出划分保持不变。
'''
    repl={'MINPEAK':f"{minimum['S_SUM_peak']:.9f}",'MINCAP':f"{minimum['cap_nm']:.0f}",
        'MINI0':f"{sweep['min_I0']['I0']:.9f}",'MINI0CAP':f"{sweep['min_I0']['cap_nm']:.0f}",'MINI0PEAK':f"{sweep['min_I0']['S_SUM_peak']:.9f}",
        'MINR':f"{sweep['min_total_R']['R']:.9f}",'MINRCAP':f"{sweep['min_total_R']['cap_nm']:.0f}",
        'N':str(stats['max_order']),'PEAK':f"{stats['S_SUM_peak']:.10f}",'ZERO':f"{zero:.10f}",
        'I0':f"{stats['I0']:.9f}",'R':f"{stats['R']:.9f}",'CONV':conv}
    for key,val in repl.items():text=text.replace('@'+key+'@',val)
    anglepath=R/'extra_gaussian_angle_diagnostic.json'
    if anglepath.exists():
        angle=json.loads(anglepath.read_text(encoding='utf8'))
        table='\n'.join(f"| {x['m']} | {x['T_q0']:.7g} | {x['T_Gaussian_mean']:.7g} | {100*x['relative_change']:+.4f}% |" for x in angle['per_order'])
        angular_text=r'''

## 4. 弱衍射是否使有限入射角宽产生数个百分点差异

这一项使用**同一 N160 截断**，把 $q=0$ 系数与 $q=-250,\ldots,250\times10^3\,\mathrm{rad/m}$、间隔 $25\times10^3\,\mathrm{rad/m}$ 的 21 点 Gaussian 入射角平均比较。主参考的 N640 正入射目标保持不变；这里没有把 N640 与 N160 的差混称角宽效应，也没有在 $K_m$ 外重新施加一个 NA 截断。

各级强度平均为 $\langle T_m\rangle=\int W(q)|r_m(q)|^2dq/\int W(q)dq$，$W=\gamma_{\rm in}(q)\exp[-(qw_0)^2/2]$。镜像对称结构满足 $r_{-m}(-q)=r_m(q)$；对于对称 Gaussian 包络，正负出射包的 parity 配对交叉项因而可用这个强度矩作诊断。该控制没有用角平均后的结果重新拟合或替换 D4C 主目标。

| m | 同 N160 的 $T_m(q=0)$ | Gaussian 平均 $\langle T_m\rangle$ | 相对变化 |
|---:|---:|---:|---:|
@ANGLE_TABLE@

弱级 m2、m12 分别有约 −3.01%、−3.42% 的相对变化，不能概括为“所有级变化都小于 1%”。但按固定 $K_m$ 和扫描核相加后，整条非零级 AC 的 L2 变化只有 **@ANGLE_L2@%**，去均值形状变化为 **@ANGLE_SHAPE@%**，peak/span 增加 **@ANGLE_PEAK@%**；零级平均相对 q0 增加 **@ANGLE_I0@%**。角域外遗漏的理想入射功率约 $2.03\times10^{-6}$。

因此，这份已完成的角敏感度控制支持“小角度依赖对弱级有可见相对影响”，但不能单独支持把数个百分点的曲线形状差异全部归因于角宽。尤其应先按弱信号自己的扫描范围核对时间收敛；100 fs 的初步差异包含已发现的时间未收敛，不能直接当作参数泛化失效。FDTD 网格、实际角域光路和本例近似 B 重建仍需分别核对；不能为了让曲线对齐而事后换一个角平均目标。全部逐级系数、两条同 N160 曲线和限定条件见 `rcwa/extra_gaussian_angle_diagnostic.json` 及新结构目录的 `Gaussian_angle_diagnostic_*`。
'''
        for key,val in {'ANGLE_TABLE':table,'ANGLE_L2':f"{100*angle['AC_relative_L2']:.6f}",
            'ANGLE_SHAPE':f"{100*angle['centered_AC_relative_L2']:.6f}",
            'ANGLE_PEAK':f"{100*angle['AC_peak_relative_change']:.6f}",
            'ANGLE_I0':f"{100*angle['I0_relative_change']:.6f}"}.items():angular_text=angular_text.replace('@'+key+'@',val)
        text+=angular_text
    (ROOT/'report/extra_weak_mark_notes.md').write_text(text,encoding='utf8')
    print('Extra weak-mark target figure/data and report notes written.')


if __name__=='__main__':main()
