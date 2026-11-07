# BF mark 与 D4C 偏振矩阵扩展：对话原文（问答阅读版）

> 讨论时间：2026-11-06 至 2026-11-07；归档日期：2026-11-07；时区：Asia/Hong_Kong（UTC+08:00）。

本文件从本聊天完整会话记录直接提取，范围从最初误称 CB mark 的提问，到本次上传原文的请求。保留更名 BF、几何澄清、用户图片、只讨论 RCWA→D4C 的范围修正、偏振矩阵推导及后续确认。原文不摘要、不改写；公式、代码、措辞、空格、换行与链接照原记录保留。编号、角色、时间及导航为归档附加排版。

用户交互回复和图片消息的原始文字封装用代码块显示。图片以原始 PNG 字节另存，并在该消息下方展示；每个原始文字块、图片及其内容顺序均在校验清单中记录。交互提问保留工具实际展示的原问题及全部选项，仅添加列表排版。

早期关于整体几何的描述随后由图片澄清，均按当时原文保留，请按顺序阅读。本次只归档讨论，没有新增 RCWA/FDTD 仿真，也不覆盖此前报告或 2026-11-06 的聊天归档。

本版共 12 条：用户 8 条、助手正式回答 3 条、进度消息 0 条、交互澄清 1 组（共 2 问）。

[归档导航](README.md) · [逐条原文与图片校验清单](20261107_bf_polarization_provenance.json)

## 正式回答导航

- [2026-11-06T22:41:42+08:00 · 第 1 次正式回答（原文）](#bf-message-0011)
- [2026-11-07T09:14:19+08:00 · 第 2 次正式回答（原文）](#bf-message-0013)
- [2026-11-07T09:16:30+08:00 · 第 3 次正式回答（原文）](#bf-message-0015)

---

<a id="bf-message-0001"></a>

## bf-message-0001 · 用户 · 2026-11-06T22:27:06+08:00

<!-- MESSAGE_BODY_BEGIN:bf-message-0001 -->
你知道CB mark这种alignment mark吗 我们已经有了RCWA复刻扫描曲线的公式 但是是对简单的mark来说 对于CBmark，横向pitch 3.2um 间隔1.6um 但是mark整体是倾斜45°的 我不知道这种情况如何求扫描曲线 应该是在一个坐标系下求正常的衍射级振幅，然后再转换到实际的坐标系下 R矩阵会有一个变化 不再象是之前那种纯单偏振情况了 你看看这个扫描曲线怎么处理 是否需要将公式拓展为偏振矩阵形式

<!-- MESSAGE_BODY_END:bf-message-0001 -->

---

<a id="bf-message-0003"></a>

## bf-message-0003 · 助手 · 交互提问 · 2026-11-06T22:29:33+08:00

<!-- MESSAGE_BODY_BEGIN:bf-message-0003 -->
你说的 CB mark 具体是哪种图形？

选项（按原顺序）：
- 二维棋盘格，方块交替排列
- 一维平行条纹，整组旋转45°
- 交叉或分段条纹组成的其他CB图形

3.2 μm 的 pitch 是沿哪个方向量的？这决定扫描相位中是否出现 √2。

选项（按原顺序）：
- 沿全局横向 x 方向测量
- 沿旋转后 Mark 自身的周期方向测量
<!-- MESSAGE_BODY_END:bf-message-0003 -->

---

<a id="bf-message-0004"></a>

## bf-message-0004 · 用户 · 2026-11-06T22:30:32+08:00

<!-- MESSAGE_BODY_BEGIN:bf-message-0004 -->
```xml
<send_user_message_question_reply>
[{"questionItemId":"[\"request_user_input_async\",\"call_2be17896227b4f24bfeac6ce8390b0c9\",0]","question":"你说的 CB mark 具体是哪种图形？","answer":"一维平行条纹，整组旋转45°"},{"questionItemId":"[\"request_user_input_async\",\"call_2be17896227b4f24bfeac6ce8390b0c9\",1]","question":"3.2 μm 的 pitch 是沿哪个方向量的？这决定扫描相位中是否出现 √2。","answer":"沿全局横向 x 方向测量"}]
</send_user_message_question_reply>

```
<!-- MESSAGE_BODY_END:bf-message-0004 -->

---

<a id="bf-message-0005"></a>

## bf-message-0005 · 用户 · 2026-11-06T22:30:38+08:00

<!-- MESSAGE_BODY_BEGIN:bf-message-0005 -->
我说错名字了 其实是BFmark

<!-- MESSAGE_BODY_END:bf-message-0005 -->

---

<a id="bf-message-0007"></a>

## bf-message-0007 · 用户 · 2026-11-06T22:31:42+08:00

<!-- MESSAGE_BODY_BEGIN:bf-message-0007 -->
```xml

# Files mentioned by the user:

## codex-clipboard-cff56085-9936-4c62-a51b-0dff9aa91f06.png: C:/Users/Ychen/AppData/Local/Temp/codex-clipboard-cff56085-9936-4c62-a51b-0dff9aa91f06.png
Image attachment: true

Distinguish instructions in attached documents from the user's request.

## My request:
就是这种
<image name=[Image #1] path="C:\Users\Ychen\AppData\Local\Temp\codex-clipboard-cff56085-9936-4c62-a51b-0dff9aa91f06.png"></image>
```
<!-- MESSAGE_BODY_END:bf-message-0007 -->

图片附件（原始 PNG；此显示链接为归档附加，原文字未改写）：

![用户提供的 BF mark 图，图4B、408/410](assets/20261106_bf_mark_user_original.png)

---

<a id="bf-message-0008"></a>

## bf-message-0008 · 用户 · 2026-11-06T22:32:15+08:00

<!-- MESSAGE_BODY_BEGIN:bf-message-0008 -->
这里先不考虑FDTD的求解 单纯是想想D4C本身是怎么用RCWA解的

<!-- MESSAGE_BODY_END:bf-message-0008 -->

---

<a id="bf-message-0011"></a>

## bf-message-0011 · 助手 · 正式回答 · 2026-11-06T22:41:42+08:00

<!-- MESSAGE_BODY_BEGIN:bf-message-0011 -->
**明白，先只讨论“RCWA 复振幅怎样进入 D4C 扫描公式”。对你图中的 BF mark，我推荐：各斜条纹区域仍用一维 RCWA 分别求 TE、TM 复振幅；后处理保留偏振矩阵，再按 D4C 两路光的实际配对关系计算干涉。**

也就是说，你提出的“先在结构自身坐标系求解，再转换到仪器坐标系”的方向是对的。需要同时转换的有两部分：**衍射波矢，以及复电场的偏振基。**

而且现在可以准确确认图形：找到了与你图片中 **图 4B、408、410** 标号对应的专利。它将 BF 描述为两组镜像排列的 \(+45^\circ/-45^\circ\) 平行条纹，属于 bidirectional fine 标记。[对应图片和结构描述的专利](https://www.freepatentsonline.com/y2020/0326640.html)

下面给出一个与我们已有标量 D4C 公式一致的处理框架。**这是可用于实现的矢量扩展；仅凭现有标量公式，还不能唯一确定你这版 D4C 内部的偏振光路。**

**第一步：把用于 RCWA 的周期设对。**

你已经确认：

\[
p_x=3.2\,\mu\mathrm m
\]

是沿全局 \(x\) 方向量到的周期。

但一维 RCWA 使用的周期，应当沿着**条纹法向**测量。因此，两组斜条纹的本体周期均为

\[
\boxed{
p_u=p_x\cos45^\circ
=\frac{3.2}{\sqrt2}
=2.262742\,\mu\mathrm m.
}
\]

如果你说的 \(1.6\,\mu\mathrm m\) 是沿 \(x\) 量到的槽宽，那么对应的法向槽宽为

\[
w_u=\frac{1.6}{\sqrt2}
=1.131371\,\mu\mathrm m,
\]

占空比仍为 50%。

取

\[
g_x=\frac{2\pi}{p_x},
\]

两组条纹的第 \(m\) 级横向波矢可以分别表示为

\[
\boxed{
\mathbf q_{+,m}=m g_x(1,1),\qquad
\mathbf q_{-,m}=m g_x(1,-1).
}
\]

这里的正负标签表示两种取向，不影响左右区域的实际命名。

如果某一路按通常的正负级次配对，即配对 \(\mathbf q_{\sigma,m}\) 和 \(-\mathbf q_{\sigma,m}\)，相对位移 \(\mathbf s=(s_x,s_y)\) 产生的干涉相位为

\[
\boxed{
\Phi_{\sigma,m}
=2\mathbf q_{\sigma,m}\cdot\mathbf s
=\frac{4\pi m}{p_x}(s_x+\sigma s_y).
}
\]

所以，**沿全局 \(x\) 扫描、保持 \(s_y\) 不变时，原公式仍可使用**

\[
\boxed{
\cos\!\left(\frac{4\pi m s_x}{3.2\,\mu\mathrm m}+\phi_{\sigma,m}\right).
}
\]

一级成对干涉的周期仍是 **\(1.6\,\mu\mathrm m\)**。不需要再给这个扫描周期额外乘一次或除一次 \(\sqrt2\)。

这里的两个 \(\sqrt2\) 已经在“法向周期”和“位移投影”中抵消了。

**第二步：每个斜条纹区域求两套复振幅，而不是只求一套。**

在各自的本体坐标中，设：

- \(u\)：条纹法向，也就是周期方向；
- \(v\)：沿条纹方向；
- \(z\)：表面法向。

对于正常入射、各向同性材料、沿条纹方向不变的一维结构，局部 TE、TM 仍然解耦。可以用一维 RCWA 分别得到

\[
r_{\mathrm{TM},m},\qquad r_{\mathrm{TE},m},
\]

并组成

\[
\boxed{
R_m^{\mathrm{loc}}
=
\begin{pmatrix}
r_{\mathrm{TM},m}&0\\
0&r_{\mathrm{TE},m}
\end{pmatrix}.
}
\]

如果左右两组具有完全相同的横截面、膜层和材料，那么通常可以复用这一组本体 TE/TM 结果，再分别处理两种取向。

**45° 旋转本身并不意味着本体结构出现了新的 TE↔TM 耦合。** 它意味着：实验室中的一个固定线偏振，在本体坐标中通常同时含有 TE 和 TM。

例如，全局 \(x\) 线偏振照射到 45° 条纹时，两个本体偏振的入射功率权重均为 \(1/2\)，但二者仍是相干的复电场分量。

求解时必须保留复数，包括相位。只保存两组衍射效率，后面就无法完整恢复偏振干涉。

**第三步：把本体偏振响应转换到仪器使用的偏振基。**

一般应写成

\[
\boxed{
R_{j,m}^{\mathrm{lab}}
=
C_{j,m}^{\mathrm{out}}\,
R_m^{\mathrm{loc}}\,
C_j^{\mathrm{in}},
}
\]

其中：

- \(j\) 表示左、右两组条纹；
- \(C_j^{\mathrm{in}}\) 把全局入射偏振转换成本体 TE/TM；
- \(C_{j,m}^{\mathrm{out}}\) 把第 \(m\) 级的出射偏振转换到仪器使用的输出基。

在近轴、已经统一输入输出横向基和反射符号的示意中，可以写成

\[
\widetilde R_m(\alpha)
=
Q(\alpha)
\begin{pmatrix}
r_{\mathrm{TM},m}&0\\
0&r_{\mathrm{TE},m}
\end{pmatrix}
Q^{T}(\alpha),
\]

\[
Q(\alpha)=
\begin{pmatrix}
\cos\alpha&-\sin\alpha\\
\sin\alpha&\cos\alpha
\end{pmatrix}.
\]

当周期法向角度取 \(\alpha=\pm45^\circ\) 时，

\[
\boxed{
\widetilde R_{m,\pm}
=
\frac12
\begin{pmatrix}
r_{\mathrm{TM},m}+r_{\mathrm{TE},m}
&
\pm(r_{\mathrm{TM},m}-r_{\mathrm{TE},m})
\\
\pm(r_{\mathrm{TM},m}-r_{\mathrm{TE},m})
&
r_{\mathrm{TM},m}+r_{\mathrm{TE},m}
\end{pmatrix}.
}
\]

这就直接说明了你关心的现象：

\[
\boxed{
\text{本体 }R_m\text{ 是对角矩阵，转换到全局偏振基后一般不再对角。}
}
\]

不过，实际实现时，不能对所有级次都盲目使用同一个二维旋转矩阵。各阶出射方向不同，其 s/p 基是三维、随方向变化的；TM 还具有相应的表面法向电场分量。最稳妥的做法是使用 RCWA 输出的完整复电场，或使用其明确规定的各阶 s/p 基，再完成统一投影。[Ansys 的 RCWA 偏振基定义](https://optics.ansys.com/hc/en-us/articles/4414575008787-RCWA-Solver-Introduction)

**第四步：D4C 的标量干涉项，扩展成两路偏振矢量之间的干涉。**

原来的公式为

\[
S_m(s)=
2K_m|r_m|^2
[1-\cos(2mGs+\phi_m)].
\]

这个式子已经隐含了较强的条件：参与干涉的两路场具有特定的幅度关系、偏振关系和相干性。

对 BF，更一般的做法是先确定：**D4C 的两路光分别将哪两个衍射分量送到了同一个探测模式。**

用 \(a,b\) 表示这一对分量，它们可以包含“左／右区域”和“正／负级次”的标签。定义两路到达共同探测偏振基后的场为

\[
\mathbf a
=
J_{1a}R_a^{\mathrm{lab}}\mathbf e_{\mathrm{in}},
\]

\[
\mathbf b
=
J_{2b}R_b^{\mathrm{lab}}\mathbf e_{\mathrm{in}}.
\]

这里 \(J_1,J_2\) 表示两路光的偏振传递，包含需要保留的分束、波片、分析器等作用。对应的差场为

\[
\mathbf D_{ab}(\mathbf s)
=
\mathbf a\,e^{-i\mathbf q_a\cdot\mathbf s}
-
\mathbf b\,e^{-i\mathbf q_b\cdot\mathbf s}.
\]

于是，一对分量的读出可以写成

\[
\boxed{
S_{ab,\chi}(\mathbf s)
=
\mathbf D_{ab}^{\dagger}
W_{ab,\chi}
\mathbf D_{ab}.
}
\]

\(W\) 表示该检测通道的功率／偏振读出权重。为了与原 D4C 标量归一化衔接，可以将原来的 \(K_m\) 放在这里，例如标量极限下取 \(W=K_mI\)。已经计入光路矩阵的权重不应再重复计算。

展开后，这个式子仍然很简单：

\[
\boxed{
S_{ab,\chi}(\mathbf s)
=
A+B
-
2\operatorname{Re}
\left[
C\,e^{i(\mathbf q_a-\mathbf q_b)\cdot\mathbf s}
\right],
}
\]

其中

\[
A=\mathbf a^\dagger W\mathbf a,\qquad
B=\mathbf b^\dagger W\mathbf b,\qquad
C=\mathbf a^\dagger W\mathbf b.
\]

也可以写成

\[
S_{ab,\chi}
=
A+B-2|C|
\cos\!\left[
(\mathbf q_a-\mathbf q_b)\cdot\mathbf s+\arg C
\right].
\]

**所以，偏振矩阵形式并不会让扫描曲线的形式变得很复杂。最终每一对干涉分量仍然只对应三个量：两路功率 \(A,B\)，以及复干涉系数 \(C\)。**

当两路场等幅、偏振一致，并满足原标量模型的条件时，

\[
A=B=K_m|r_m|^2,\qquad
|C|=K_m|r_m|^2,
\]

就恢复我们已有的

\[
2K_m|r_m|^2[1-\cos(\cdots)].
\]

这也解释了为什么不能简单将矩阵换成一个“总反射强度”：**总强度只能告诉你 \(A,B\)，不能告诉你决定干涉幅度和相位的 \(C\)。**

**第五步：对图中的左右两组条纹，要保留两组标签，不能直接假定它们互相干涉。**

具体配对由干涉仪的空间变换决定。例如：

- 若两路相对是 \(180^\circ\) 旋转，通常涉及 \(\mathbf q\) 与 \(-\mathbf q\) 的配对。
- 若两路相对是某个镜像变换，可能配对来自不同取向区域的衍射光。

\(180^\circ\) 旋转不会把 \(+45^\circ\) 条纹方向变成 \(-45^\circ\)。因此，不能仅凭左右镜像排列，就默认“左边的一级与右边的一级相干相加”。

公开的 BF 相关描述也明确区分了标记结构和传感器读出；一些传感器使用两路相对旋转 \(180^\circ\) 的自参考干涉，并分别形成和场、差场。[BF 与自参考干涉读出的公开描述](https://patents.google.com/patent/US20260118787A1/en)

两组信号只有在独立检测、输出模式互不重叠，或相应交叉项确实为零时，才能直接加功率。组合前还要统一各单元的空间原点和参考面，保留左右光栅的注册相位。

**对于你当前想做的最简 RCWA—D4C 扩展，我建议按这个顺序实现：**

1. **按 \(2.262742\,\mu\mathrm m\) 的法向周期建立一维 RCWA 单元。**
2. **分别计算 TE、TM，并保留各正负级次的复振幅。**
3. **为左右两种取向分别构造偏振响应矩阵和实际衍射波矢。**
4. **按 D4C 的两路映射关系，计算每对分量的 \(A,B,C\)。**
5. **求和得到扫描曲线，再按原规则加入零级残余和漏光项。**

如果进一步确认两路不混合本体 TE/TM、检测不分析偏振，而且两种偏振使用相同的标量权重，那么可以简化成

\[
S_j
=
|e_{\mathrm{TM},j}|^2S_{\mathrm{TM},j}
+
|e_{\mathrm{TE},j}|^2S_{\mathrm{TE},j}.
\]

对全局 \(x\) 线偏振和 45° 条纹，就是

\[
S_j=\frac12S_{\mathrm{TM},j}+\frac12S_{\mathrm{TE},j}.
\]

**这个可以作为最简单的对照模型，但它成立需要上述光路条件；建议程序内部先保留 Jones 矩阵，不要一开始就丢掉 TE/TM 之间的复相位。**

另外，若仍使用 532 nm、NA＝0.6，新的周期下几何上只收集到正负一级和二级。此前 16 μm 周期的 \(K_m\) 表对应不同的衍射角与瞳位置，不能仅凭相同的级次编号直接搬过来。这里需要使用 BF 对应的读出权重，或明确建立随瞳位置、偏振变化的权重模型。
<!-- MESSAGE_BODY_END:bf-message-0011 -->

---

<a id="bf-message-0012"></a>

## bf-message-0012 · 用户 · 2026-11-07T09:13:23+08:00

<!-- MESSAGE_BODY_BEGIN:bf-message-0012 -->
所以实际上是不是和一个本身入射椭圆偏振的正常x轴周期的mark的得到R矩阵后的处理方式基本相似

<!-- MESSAGE_BODY_END:bf-message-0012 -->

---

<a id="bf-message-0013"></a>

## bf-message-0013 · 助手 · 正式回答 · 2026-11-07T09:14:19+08:00

<!-- MESSAGE_BODY_BEGIN:bf-message-0013 -->
**是的。如果你指的是“RCWA 已经给出每个衍射级的复 \(R\) 矩阵之后，怎样计算两路光的干涉和扫描信号”，两者基本就是同一套 Jones 矢量处理。** 可以共用一个后处理程序，不需要为旋转 mark 和椭圆偏振入射分别写两套扫描公式。

不过，需要区分：**处理方式相同，几何和入射状态并不因此等价。**

对于沿正常 \(x\) 方向周期的 mark，任意入射偏振都可以写成

\[
\mathbf e_{\mathrm{in}}
=
\begin{pmatrix}
e_x\\e_y
\end{pmatrix},
\]

其中 \(e_x,e_y\) 可以是带相对相位的复数，因此也包含椭圆偏振。第 \(m\) 级输出为

\[
\boxed{
\mathbf E_m=R_m\mathbf e_{\mathrm{in}}.
}
\]

对旋转后的条纹，同样有

\[
\boxed{
\mathbf E_m^{\mathrm{lab}}
=
R_m^{\mathrm{lab}}\mathbf e_{\mathrm{in}}^{\mathrm{lab}}
=
C_m^{\mathrm{out}}R_m^{\mathrm{loc}}
\underbrace{C^{\mathrm{in}}\mathbf e_{\mathrm{in}}^{\mathrm{lab}}}
_{\mathbf e_{\mathrm{in}}^{\mathrm{loc}}}.
}
\]

所以，你可以把它理解成：

> 在条纹本体坐标中，使用转换后的入射 Jones 矢量求出衍射场，再把出射场送入共同的仪器偏振坐标；接下来的干涉计算，与任意椭圆偏振入射时使用的计算完全相同。

**但 45° 坐标旋转本身不会把线偏振变成椭圆偏振。**

例如，实验室中的 \(x\) 线偏振为

\[
\mathbf e_{\mathrm{in}}^{\mathrm{lab}}
=
\begin{pmatrix}1\\0\end{pmatrix}.
\]

转换到某一组 45° 条纹的本体坐标后，可以得到

\[
\mathbf e_{\mathrm{in}}^{\mathrm{loc}}
=
\frac1{\sqrt2}
\begin{pmatrix}1\\-1\end{pmatrix},
\]

它仍是线偏振，只是同时激励了局部 TM 和 TE。经过结构反射后，

\[
\mathbf E_m^{\mathrm{loc}}
=
\frac1{\sqrt2}
\begin{pmatrix}
r_{\mathrm{TM},m}\\
-r_{\mathrm{TE},m}
\end{pmatrix}.
\]

如果 \(r_{\mathrm{TM},m}\) 和 \(r_{\mathrm{TE},m}\) 具有不同的复相位，**出射光才可能变为椭圆偏振**。因此，这里要保留两分量的相对复相位，不能提前只留下两份强度。

拿到这些输出场以后，两种情况都可以使用同一个扫描计算模块。设经过两路光学传递后，参与干涉的复矢量为 \(\mathbf a,\mathbf b\)，则

\[
\boxed{
S(\mathbf s)
=
A+B
-
2\operatorname{Re}
\left[
C\,e^{i(\mathbf q_a-\mathbf q_b)\cdot\mathbf s}
\right]
+Z,
}
\]

其中

\[
A=\mathbf a^\dagger W\mathbf a,\qquad
B=\mathbf b^\dagger W\mathbf b,\qquad
C=\mathbf a^\dagger W\mathbf b.
\]

**无论偏振的两分量来自椭圆入射，还是来自旋转条纹对 TE/TM 的共同激励，最终都是计算这三个量。** 偏振对比度与附加相位，都通过复干涉系数 \(C\) 体现。

对你的 BF mark，实现时还应保留三点：

- **左右两组分别转换。** 它们的取向是 \(+45^\circ/-45^\circ\)，对应的入射基转换、出射基转换和波矢不同。
- **扫描相位使用实际波矢差。** 偏振矩阵处理可以共用，几何相位仍按 \((\mathbf q_a-\mathbf q_b)\cdot\mathbf s\) 计算。
- **坐标只转换一次。** 如果 \(R_m\) 已经转换到全局坐标，就直接乘全局入射矢量；若使用本体 \(R_m\)，才先转换入射矢量，并在输出端转换回来。

因此，你的这个类比很适合用来组织程序：**先把“RCWA＋坐标转换”统一成各阶的复矢量输出，再复用“任意偏振两路干涉”的 D4C 后处理。**
<!-- MESSAGE_BODY_END:bf-message-0013 -->

---

<a id="bf-message-0014"></a>

## bf-message-0014 · 用户 · 2026-11-07T09:16:01+08:00

<!-- MESSAGE_BODY_BEGIN:bf-message-0014 -->
所以我们实际上要做的就是把公式扩展到偏振矩阵而非单一偏振的形式

<!-- MESSAGE_BODY_END:bf-message-0014 -->

---

<a id="bf-message-0015"></a>

## bf-message-0015 · 助手 · 正式回答 · 2026-11-07T09:16:30+08:00

<!-- MESSAGE_BODY_BEGIN:bf-message-0015 -->
**对。就你现在“用 RCWA 复振幅生成 D4C 扫描曲线”的目标，核心改动就是：把单一偏振的标量振幅公式，扩展成保留两偏振相干关系的 Jones 矩阵／矢量公式。**

具体对应关系是：

| 原来的单偏振处理 | 扩展后的处理 |
|---|---|
| 复振幅 \(r_m\) | 每级 Jones 矩阵 \(R_m\) |
| 一个入射振幅 | 入射偏振矢量 \(\mathbf e_{\mathrm{in}}\) |
| 两个标量场相减 | 两个经过光路变换的偏振矢量场相减 |
| 对标量场取模平方 | 对矢量场按检测偏振和权重计算功率 |
| 相位 \(2mGs\) | 实际波矢差产生的相位 \((\mathbf q_a-\mathbf q_b)\cdot\mathbf s\) |

一对参与干涉的衍射分量，可以统一写成

\[
\boxed{
\mathbf D(\mathbf s)
=
J_1R_a\mathbf e_{\mathrm{in}}\,
e^{-i\mathbf q_a\cdot\mathbf s}
-
J_2R_b\mathbf e_{\mathrm{in}}\,
e^{-i\mathbf q_b\cdot\mathbf s}.
}
\]

然后计算

\[
\boxed{
S_\chi(\mathbf s)
=
\sum_{\text{参与干涉的各对}}
\mathbf D^\dagger(\mathbf s)\,
W_\chi\,
\mathbf D(\mathbf s)
+
Z_\chi.
}
\]

这里 \(J_1,J_2\) 表示两路的偏振传递，\(W_\chi\) 表示检测通道的读出权重，\(Z_\chi\) 保留对应的零级残余与漏光贡献。

**要注意，不能只把原公式中的 \(|r_m|^2\) 换成一个“矩阵的模平方”。** 必须先形成两路复矢量场并相干叠加，再计算功率；这样才能保留 TE/TM 的相对相位，以及两路偏振是否匹配对干涉对比度的影响。

在单偏振、两路等幅且偏振一致的极限下，取相应权重 \(W=K_mI\)，上式会恢复原来的

\[
2K_m|r_m|^2[1-\cos(2mGs+\phi_m)].
\]

对你的 BF，**RCWA 本体求解仍可先用一维 TE、TM 两套复振幅；主要需要扩展的是坐标变换和 D4C 干涉读出部分。** 左右 \(+45^\circ/-45^\circ\) 两组分别构造响应，再按实际光路配对。
<!-- MESSAGE_BODY_END:bf-message-0015 -->

---

<a id="bf-message-0016"></a>

## bf-message-0016 · 用户 · 2026-11-07T09:17:39+08:00

<!-- MESSAGE_BODY_BEGIN:bf-message-0016 -->
把以上几段对话原文上传至github

<!-- MESSAGE_BODY_END:bf-message-0016 -->

---

