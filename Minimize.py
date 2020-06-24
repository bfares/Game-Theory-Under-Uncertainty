# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 20:30:08 2019

@author: User
"""


from scipy.optimize import minimize
import numpy as np

def minimize2():
        c= [-5, 4]
        c0=[2]
        x0_bounds = (0, None)
        x1_bounds = (0, None)
        res=optimize.linprog(c=c, method="interior-point", bounds=(x0_bounds, x1_bounds) )
        #res=optimize.linprog(c,A_ineq,b_ineq, method="interior-point" )
        if res.success:
            print("LP Solved")
        else:
            print("LP failed:", res.message)
        #print('Optimal value:', res.fun, '\nX:', res.x)
        return res.x



def eq2( x ):
    
    f1 = 2 - 5 * x[0] + 4 * x[1]
    f2 = 3 + 0 * x[0] - x[1]
    f3 = 4 + 5 * x[0] - 5 * x[1]
    return (f1, f2, f3) 

bnds = ( (0, None), (0, None) )
cons = ( { 'type' : 'ineq', 'fun': lambda x: x[0]+[p1]+[p2] - 1} )

x0=np.array([0.3,0.3])




def eq( p ):
    s1,s2,s3 = p 
    f1 = 1.1**3 / s1*1.1**1+s2*1.1**2+s3*1.1**3
    f2 = 0.9**1 / s1*0.9**1+s2*0.9**2+s3*0.9**3
    return (f1, f2) 

bnds = ( (0, None), (0, None), (0, None) )
cons = ( { 'type' : 'ineq', 'fun': lambda p: p[0]+p[1]+p[2] - 1} )


res=minimize( eq, [0.3,0.3,0.3],  bounds=bnds,method='SLSQP', constraints=cons )
print (res)
