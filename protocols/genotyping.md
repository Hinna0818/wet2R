
# Mouse Genotyping Protocol（小鼠基因型鉴定）

## 一、基因组DNA提取（Genomic DNA Extraction）

### Step 1：组织处理

1. 将已采集的小鼠脚趾或尾部组织放入 **1.5 mL EP管**。
2. 在EP管上标记对应的小鼠编号。
3. 每管加入 **200 μL 50 mM NaOH溶液**。

### Step 2：组织裂解

4. 将金属浴（Dry Bath）预热至 **100℃**。
5. 将EP管放入金属浴，加热至少 **15 min**，避免加热时间过长。
6. 加热结束后，短暂离心，使管盖上的液滴回到管底。

### Step 3：中和与DNA收集

7. 每管加入 **80 μL Tris-HCl（pH 8.1）**。
8. 充分混匀，得到含有基因组DNA的粗提取液。
9. DNA可直接用于后续PCR，或暂时保存于 **4℃**。

---

## 二、PCR扩增（Genotyping PCR）

### Step 1：配制PCR反应体系

每管总反应体积为 **15 μL**。

| 试剂 | 每管体积 |
|---|---:|
| 2× PCR Mix | 7.5 μL |
| Primer Mix（Forward + Reverse，1:1） | 1 μL |
| ddH₂O | 5 μL |
| Genomic DNA Template | 1.5 μL |
| **Total** | **15 μL** |

### Step 2：配置PCR反应

1. 根据需要鉴定的小鼠数量及基因型，准备对应的PCR管或PCR板。
2. 按照检测目标选择相应的Forward和Reverse引物。
3. 将PCR Mix、Primer Mix及ddH₂O混合，配制Master Mix。
4. 每管加入 **13.5 μL Master Mix**。
5. 每管加入 **1.5 μL对应小鼠的基因组DNA**。
6. 轻轻混匀，短暂离心，使液体集中在管底。

### Step 3：PCR扩增

1. 将PCR管或PCR板放入PCR仪。
2. 根据所用引物及PCR Mix设置扩增程序。
3. 反应结束后，取PCR产物进行基因型检测。

---

## 三、琼脂糖凝胶电泳（Agarose Gel Electrophoresis）

> 本部分为补充的常规检测流程，原始Protocol未提供具体条件。

1. 根据预期PCR产物大小，选择合适浓度的琼脂糖凝胶。
2. 将PCR产物与适当的Loading Buffer混合（若PCR Mix已含上样染料，可按说明书操作）。
3. 将DNA Ladder（分子量标准）及PCR产物分别加入凝胶孔。
4. 进行电泳，使不同长度的DNA片段分离。
5. 使用凝胶成像系统观察DNA条带。
6. 根据预期条带大小，结合阳性、阴性对照判断小鼠基因型。

例如，对于经过验证的Flox鉴定引物，可根据WT和Flox等位基因的不同扩增片段区分：

- WT/WT：野生型纯合。
- WT/Flox：杂合。
- Flox/Flox：Flox纯合。

具体条带大小必须根据对应小鼠品系和引物确定。

---

## 四、注意事项

1. 每个EP管必须标记小鼠编号，避免样本混淆。
2. NaOH用于组织裂解，Tris-HCl用于中和裂解液。
3. DNA粗提液中可能含有组织碎片，取PCR模板时尽量避免吸入大量沉淀。
4. 不同基因型（如Flox、Cre）可能需要分别配置PCR引物。
5. 建议设置已知基因型的阳性对照和无模板对照（NTC）。
6. PCR扩增程序及预期条带大小应使用经过验证的品系专用方案。
