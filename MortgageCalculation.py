#!/usr/bin/python
# -*- coding: utf-8 -*-

from sympy import *
import argparse
#贷款额；年利率；贷款期数

#等额本金/等额本息

def parse_args():
    parser = argparse.ArgumentParser(description="计算房贷利率算法，等额本金/等息本金")
    parser.add_argument("-l", "--loan_amount", type=int, required=True, help="输入贷款总额，如120000")
    parser.add_argument("-i", "--annual_interest_rate", type=float, required=True, help="输入年化利率，如0.06")
    parser.add_argument("-s", "--stage_num", type=int, required=True,help="输入贷款期数")
    args = parser.parse_args()
    return args

def equal_principal_payment(loan_amount, annual_interest_rate, stage_num):
    '''
    #贷款额：120000； 年利率：6%；偿还期数：12
    #本金：10000/M
    #1. 120000 * 0.5% = 600
    #2. 110000 * 0.5% = 550
    #......
    #12. 10000 * 0.5% = 50
    '''
    fixed_loan_per_month = float(loan_amount / stage_num)
    IR_perM = annual_interest_rate / 12.0
    total_interest = sum([(loan_amount - fixed_loan_per_month * i) * IR_perM for i in range(0, stage_num)])
    return [fixed_loan_per_month + (loan_amount - fixed_loan_per_month * i) * IR_perM for i in range(0, stage_num)]

def equal_loan_payment(loan_amount, annual_interest_rate, stage_num):
    '''
    #贷款额：120000； 年利率：6%；偿还期数：12
    #每月欠款分别为，a0, a1, a2, , , , , a12；每月还款为A（包括本金和利息）
    #a0 = 120000
    #a1 = a0 * (1+0.5%) - A
    #a2 = a1 * (1+0.5%) - A
    #......
    #a12 = a11 * (1+0.5%) - A = 0
    #求A即可
    '''
    A = symbols('A')
    IR_perM = annual_interest_rate / 12.0
    #a0 = loan_amount
    res_expr = loan_amount * (1 + IR_perM) - A
    for _ in range(1, stage_num):
        res_expr = res_expr * (1 + IR_perM) - A
    res = solve(res_expr, A)
    return res[0]

def main():

    args = parse_args()

    EPP = equal_principal_payment(loan_amount=args.loan_amount, annual_interest_rate=args.annual_interest_rate, stage_num=args.stage_num)
    print("#" * 50)
    print("贷款额：{} ； 年化利率：{}；还款期数：{}".format(args.loan_amount, args.annual_interest_rate, args.stage_num))
    print("#" * 50)
    print('等额本金贷款：')

    each_M = args.loan_amount/args.stage_num
    for i in range(0, args.stage_num):
        print('第{}期，本金：{:.2f}, 利息：{:.2f}'.format(i+1, EPP[i], (EPP[i] - each_M)))

    total_inter = sum(EPP) - args.loan_amount

    print('#总利息{:.2f}'.format(total_inter))

    print("#"*50)
    ELP = equal_loan_payment(loan_amount=args.loan_amount, annual_interest_rate=args.annual_interest_rate, stage_num=args.stage_num)
    print('等额本息贷款：')
    print('每月还款额：{:.2f}'.format(ELP))
    print('#总利息{:.2f}'.format(ELP*args.stage_num - args.loan_amount))
    print("#" * 50)
if __name__ == '__main__':
    main()

