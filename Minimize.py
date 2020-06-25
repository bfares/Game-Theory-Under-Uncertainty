# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 20:30:08 2019

@author: Bernard Fares
"""
from scipy import optimize

class Minimize(object):
    A_eq=""
    b_eq=""
    A_ineq=""
    b_ineq=""
    c=""
    ax=""
    colour=""
    
    def __init__(self,c, A_ineq, b_ineq, A_eq, b_eq):
        self.A_eq= A_eq
        self.b_eq= b_eq
        self.A_ineq= A_ineq
        self.b_ineq= b_ineq
        self.c= c
        
    def compute(self):
        A_eq=self.A_eq 
        b_eq=self.b_eq 
        A_ineq = self.A_ineq 
        b_ineq=self.b_ineq 
        c=self.c
        
        res=optimize.linprog(c=c, A_ub=A_ineq, b_ub=b_ineq, A_eq=A_eq, b_eq=b_eq, method='interior-point', bounds=(0, 1) )
        if not res.success:
            print("LP failed:", res.message)
        #print('Optimal value:', res.fun, '\nX:', res.x)
        return res.x