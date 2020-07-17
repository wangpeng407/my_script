
## 1 MortgageCalculation.py
### 1.1 usage
```
python MortgageCalculation.py -h
usage: MortgageCalculation.py [-h] -l LOAN_AMOUNT -i ANNUAL_INTEREST_RATE -s
                              STAGE_NUM

计算房贷利率算法，等额本金/等息本金

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
##################################################
贷款额：120000 ； 年化利率：0.06；还款期数：12
##################################################
等额本金贷款：
第1期，本金：10600.00, 利息：600.00
第2期，本金：10550.00, 利息：550.00
第3期，本金：10500.00, 利息：500.00
第4期，本金：10450.00, 利息：450.00
第5期，本金：10400.00, 利息：400.00
第6期，本金：10350.00, 利息：350.00
第7期，本金：10300.00, 利息：300.00
第8期，本金：10250.00, 利息：250.00
第9期，本金：10200.00, 利息：200.00
第10期，本金：10150.00, 利息：150.00
第11期，本金：10100.00, 利息：100.00
第12期，本金：10050.00, 利息：50.00
#总利息3900.00
##################################################
等额本息贷款：
每月还款额：10327.97
#总利息3935.66
##################################################
```
