
# RT-qPCR Protocol

## 一、First Strand Synthesis（RNA反转录）

### Step 1：RNA定量与配制

1. 使用NanoDrop测量RNA浓度（ng/μL），一般是取1μL的RNA溶液到NanoDrop里。
2. 每个样本取 **2 μg RNA（2000 ng）**。
3. 使用DEPC水将RNA补足至 **12 μL**。

RNA取样体积计算：

RNA体积（μL）= 2000 / RNA浓度（ng/μL）

例如，假设通过NanoDrop测的RNA浓度为500 ng/μL：
- RNA：4 μL
- DEPC水：8 μL
- 总体积：12 μL

这一步要拿出笔记本，算具体的浓度数值。具体来说，先记录RNA-aq的浓度c1，然后2000/c2得到需要加的RNA-aq的体积v1，然后再用12-v1得到水的体积v2。

然后算完了之后，要把ep管里的RNA溶液根据算出的体积转移到pcr管里，并用水补齐12μL体系。

加液体的时候先大后小。

注意：
RNase-free的枪头要拿外面的，每次测量完浓度都要擦一下上下两个头。

### Step 2：去除基因组DNA（gDNA Digestion）

1. 向上述12 μL RNA体系中加入：
   - RNA溶液：12 μL
   - 5× Yeason gDNA Digestion Mix：3 μL
   - **总体积：15 μL**

2. 轻轻混匀，放入PCR仪。

3. 设置反应程序：
   - 热盖温度：60℃
   - 反应体积：15 μL
   - 42℃，2 min

### Step 3：反转录合成cDNA

1. 向上述15 μL体系中加入：
   - gDNA消化后的RNA溶液：15 μL
   - 4× Hifair III SuperMix plus：5 μL
   - **总体积：20 μL**

2. 轻轻混匀，放入PCR仪。

3. 设置反应程序：
   - 热盖温度：105℃
   - 反应体积：20 μL

| 温度 | 时间 |
|---|---|
| 50℃ | 5 min |
| 85℃ | 5 s |

4. 反应结束后，获得cDNA，用于后续qPCR。

---

## 二、qPCR（实时荧光定量PCR）

### Step 1：配制qPCR反应体系

每孔总反应体积为 **15 μL**。

| 试剂 | 每孔体积 |
|---|---:|
| 2× SYBR Green qPCR Mix | 7.5 μL |
| Primer Mix（Forward + Reverse，1:1） | 0.5 μL |
| ddH₂O | 6 μL |
| cDNA Template | 1 μL |
| **Total** | **15 μL** |

### Step 2：加入试剂

1. 根据样本数量和检测基因数量计算所需试剂。
2. 将SYBR Green qPCR Mix、Primer Mix和ddH₂O混合，配制Master Mix。
3. 每孔加入 **14 μL Master Mix**。
4. 每孔加入 **1 μL对应样本的cDNA**。
5. 轻轻混匀，避免产生气泡。
6. 封板，准备放入qPCR仪。

### Step 3：qPCR上机检测

1. 将反应板放入实时荧光定量PCR仪。
2. 选择与SYBR Green体系匹配的检测程序。
3. 进行PCR扩增及实时荧光检测。
4. 反应结束后，查看扩增曲线、Ct/Cq值及熔解曲线。

---

## 三、注意事项

1. RNA操作尽量使用无RNase污染的耗材。
2. 不同样本尽量保持相同的RNA反转录投入量。
3. Primer Mix的实际浓度需要根据实验室引物工作液确认。
4. 配制多孔qPCR反应时，建议使用Master Mix，减少移液误差。
5. 建议设置技术重复和阴性对照（NTC）；必要时设置无逆转录对照（No-RT）。
6. 原Protocol未说明cDNA是否需要稀释，正式实验前应确认。
