
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

