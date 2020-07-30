
## 1 MortgageCalculation.py
### 1.1 usage
```
python MortgageCalculation.py -h
usage: MortgageCalculation.py [-h] -l LOAN_AMOUNT -i ANNUAL_INTEREST_RATE -s
                              STAGE_NUM

计算房贷利率算法，等额本金/等额本息

optional arguments:
  -h, --help            show this help message and exit
  -l LOAN_AMOUNT, --loan_amount LOAN_AMOUNT
                        输入贷款总额，如120000
  -i ANNUAL_INTEREST_RATE, --annual_interest_rate ANNUAL_INTEREST_RATE
                        输入年化利率，如0.06
  -s STAGE_NUM, --stage_num STAGE_NUM
                        输入贷款期数
```

### 1.2 example
```
python MortgageCalculation.py -l 120000 -i 0.06 -s 12
```

### 1.3 output

```
**************************************************
**************************************************
贷款额：120000 ； 年化利率：0.06；还款期数：12
##################################################
等额本金贷款：
第1期，本金：10000.00, 利息：600.00
第2期，本金：10000.00, 利息：550.00
第3期，本金：10000.00, 利息：500.00
第4期，本金：10000.00, 利息：450.00
第5期，本金：10000.00, 利息：400.00
第6期，本金：10000.00, 利息：350.00
第7期，本金：10000.00, 利息：300.00
第8期，本金：10000.00, 利息：250.00
第9期，本金：10000.00, 利息：200.00
第10期，本金：10000.00, 利息：150.00
第11期，本金：10000.00, 利息：100.00
第12期，本金：10000.00, 利息：50.00
#总利息3900.00
##################################################
等额本息贷款：
每月还款额：10327.97
第1期，本金：9727.97, 利息：600.00
第2期，本金：9776.61, 利息：551.36
第3期，本金：9825.49, 利息：502.48
第4期，本金：9874.62, 利息：453.35
第5期，本金：9924.00, 利息：403.98
第6期，本金：9973.62, 利息：354.36
第7期，本金：10023.48, 利息：304.49
第8期，本金：10073.60, 利息：254.37
第9期，本金：10123.97, 利息：204.00
第10期，本金：10174.59, 利息：153.38
第11期，本金：10225.46, 利息：102.51
第12期，本金：10276.59, 利息：51.38
#总利息3935.66
##################################################
```

## 2. clean_reads_extract.pl

### 2.1 Usage
```
Usage: script for extracting paired sequences for MGS/MLG according to bowtie mapping results

```

### 2.2 Examples 
```
#building index
bwa index total.bin.fa

#extracting reads of  bins/mgs/mlg
bwa mem -t 10 total.bin.fa sample_1.fastq.gz sample_2.fastq.gz 2> bwa.log | perl clean_reads_extract.pl -1 sample_1.fastq.gz -2 sample_2.fastq.gz -mgs bin.cluster  -bam -  --outdir ./
```

## 3 coverage_calculate_from_genomeCoverageBed.pl

## 4 search_paper_info.py ***Searching paper infomation according to keywords and title***

### 4.1 Usage
```
python search_paper_info.py  -h

usage: search_paper_info.py [-h] -l LIST [-m MAXITERM]

Rerieve published paper infomation from pubmed
(https://pubmed.ncbi.nlm.nih.gov/) according to article title or keywords.

optional arguments:
  -h, --help            show this help message and exit
  -l LIST, --list LIST  input list include article title or keywords
  -m MAXITERM, --maxiterm MAXITERM
                        Max iterms when using esearch function, default is 20
```

```
python search_paper_info.py -l plist -m 10
```

```
#All pubmed ids: 32345342, 32698005

Pubmed ID: 32345342

Title: Influenza infection elicits an expansion of gut population of endogenous Bifidobacterium animalis which protects mice against infection.

Journal: Genome biology (Genome Biol.)

Date: 2020-04

Abstrct_EN: Influenza is a severe respiratory illness that continually threatens global health. It has been widely known that gut microbiota modulates the host response to protect against influenza infection, but mechanistic details remain largely unknown. Here, we took advantage of the phenomenon of lethal dose 50 (LD<sub>50</sub>) and metagenomic sequencing analysis to identify specific anti-influenza gut microbes and analyze the underlying mechanism. Transferring fecal microbes from mice that survive virulent influenza H7N9 infection into antibiotic-treated mice confers resistance to infection. Some gut microbes exhibit differential features to lethal influenza infection depending on the infection outcome. Bifidobacterium pseudolongum and Bifidobacterium animalis levels are significantly elevated in surviving mice when compared to dead or mock-infected mice. Oral administration of B. animalis alone or the combination of both significantly reduces the severity of H7N9 infection in both antibiotic-treated and germ-free mice. Functional metagenomic analysis suggests that B. animalis mediates the anti-influenza effect via several specific metabolic molecules. In vivo tests confirm valine and coenzyme A produce an anti-influenza effect. These findings show that the severity of influenza infection is closely related to the heterogeneous responses of the gut microbiota. We demonstrate the anti-influenza effect of B. animalis, and also find that the gut population of endogenous B. animalis can expand to enhance host influenza resistance when lethal influenza infection occurs, representing a novel interaction between host and gut microbiota. Further, our data suggest the potential utility of Bifidobacterium in the prevention and as a prognostic predictor of influenza.

Abstrct_CN: 流感是一种严重的呼吸道疾病，不断威胁着全球健康。众所周知，肠道菌群可调节宿主反应以预防流感感染，但机制细节仍然未知。在这里，我们利用了致命剂量50（LD <sub> 50 </ sub>）的现象和宏基因组测序分析来识别特定的抗流感肠道微生物并分析其潜在机制。将粪便微生物从在H7N9流感强毒中幸存下来的小鼠转移到抗生素治疗的小鼠中，使其具有抗药性。一些肠道微生物表现出与致命流感感染不同的特征，具体取决于感染的结果。与死亡或模拟感染的小鼠相比，存活小鼠中的假双歧杆菌和动物双歧杆菌水平显着升高。单独口服动物双歧杆菌或两者结合可显着降低抗生素治疗和无菌小鼠中H7N9感染的严重程度。功能宏基因组学分析表明，动物双歧杆菌通过几种特定的代谢分子介导抗流感作用。体内试验证实缬氨酸和辅酶A产生抗流感作用。这些发现表明，流感感染的严重程度与肠道菌群的异质性反应密切相关。我们证明了动物双歧杆菌的抗流感作用，并且还发现当致死性流感感染发生时，内源性动物双歧杆菌的肠道种群可以扩展以增强宿主对流感的抵抗力，代表宿主与肠道菌群之间的新型相互作用。此外，我们的数据表明双歧杆菌在预防流感和作为流感的预后指标方面具有潜在的实用性。

####################################################################################################

Pubmed ID: 32698005

Title: Local Necrotic Cells Trigger Systemic Immune Activation via Gut Microbiome Dysbiosis in Drosophila.

Journal: Cell reports (Cell Rep)

Date: 2020-Jul

Abstrct_EN: Necrotic cells elicit an inflammatory response through their endogenous factors with damage-associated molecular patterns. Blocking apoptosis in Drosophila wings leads to the necrosis-driven systemic immune response by unknown mechanisms. Here, we demonstrate that immune activation in response to necrotic cells is mediated by commensal gut microbiota. Removing the microbiome attenuates hyperactivation of the innate immune signaling IMD pathway in necrosis-induced flies. Necrotic cells in wings trigger Gluconobacter expansion in the gut. An isolated Gluconobacter sp. strain is sufficient for pathological IMD activation in necrosis-induced flies, while it is not inflammatory for control animals. In addition, bacterial colonization shifts the host metabolome and shortens the lifespan of necrosis-induced flies. This study shows that local necrosis triggers a pathological systemic inflammatory response through interaction between the host and the dysbiotic gut microbiome.

Abstrct_CN: 坏死细胞通过其内源性因子与损伤相关的分子模式引起炎症反应。果蝇翅膀的凋亡被未知机制导致坏死驱动的全身免疫反应。在这里，我们证明了响应性坏死细胞的免疫激活是由共生肠道菌群介导的。去除微生物组减弱了坏死诱导的果蝇中先天免疫信号IMD途径的过度激活。翅膀中的坏死细胞触发肠道内的葡糖杆菌膨胀。分离的葡糖杆菌该菌株足以使坏死诱导的果蝇发生病理性IMD活化，而对对照动物却没有炎症。此外，细菌定植会改变宿主的代谢组并缩短坏死诱导的果  1 Influenza infection elicits an expansion of gut population of endogenous Bifidobacterium animalis which protects mice against infection
蝇的寿命。这项研究表明，局部坏死通过宿主与营养不良肠道微生物组之间的相互作用触发了病理性全身炎症反应。

####################################################################################################
```
