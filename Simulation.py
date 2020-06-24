# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 22:49:34 2019

@author: Bernard Fares
"""
from Polytope import Polytope
from Game import Game
import numpy as np
import matplotlib.pyplot as plt

#Game examples
chicken_p1 = np.array( [[4,1], [5,0]])
chicken_p2 = np.array( [[4,5], [1,0]])
battleOfTheSexes_p1 = np.array( [[2,0], [0,1]])
battleOfTheSexes_p2= np.array( [[1,0], [0,2]])
shapley_p1 = np.array( [[1,0,0], [0,1,0], [0,0,1]])
shapley_p2 = np.array( [[0,1,0], [0,0,1], [1,0,0]])
matchingPennies_p1_RiskAverse = np.array( [[0.293,-0.414], [-0.414,0.293]])
matchingPennies_p2_RiskAverse = np.array( [[-0.414,0.293], [0.293,-0.414]])
matchingPennies_p1_MarginalUtilities = np.array( [[0.245,0.49], [0.49,0.245]])
matchingPennies_p2_MarginalUtilities = np.array( [[0.49,0.245], [0.245,0.49]])

Thesis_Test_p1 = np.array( [['f',0], [0,1]])
Thesis_Test_p2 = np.array( [[1,0], [0,2]])

lP={'f':1.33}
uP={'f':1.83}

#cGame=Game(chicken_p1, chicken_p2)
#cGame=Game(shapley_p1, shapley_p2)
#cGame=Game(battleOfTheSexes_p1, battleOfTheSexes_p2)
cGame=Game(Thesis_Test_p1, Thesis_Test_p2, lP, uP)


#ineq_A=cGame.A_ineq_RiskAverse(matchingPennies_p1_MarginalUtilities,matchingPennies_p2_MarginalUtilities)
A_ineq= cGame.A_ineq()
c=cGame.c()
b_ineq= cGame.b_ineq()
A_eq= cGame.A_eq()
b_eq= cGame.b_eq()
print("Game CE constraints: \n", A_ineq)
print ("Game objective function: \n",c)
cPolytope=Polytope(c, A_ineq, b_ineq , A_eq,b_eq)

CEqui=cPolytope.minimize()
print ("Game Correlated Equilibria: \n", CEqui)
vertices=cPolytope.vertices()
print ("Polytope Vertices: \n", vertices)
colourPolytope=(0.6,1,0.6, 0.1)
cPolytope.plotConvexHull(colourPolytope)
#cPolytope.plotConvexHull()
#print(np.hstack((cGame.A_ineq() , cGame.b_ineq())))
#print(cGame.c())
 #cGame.c()
identity_A_ineq= cGame.identity_A_ineq()
identity_b_ineq= cGame.identity_b_ineq()
cSimplex=Polytope(c, identity_A_ineq, identity_b_ineq, cGame.A_eq(),cGame.b_eq(), cPolytope.getFigure())
colourSimplex=(0.878,0.878,0.878, 0.05)
cSimplex.plotConvexHull(colourSimplex)


ax = cPolytope.getFigure()
ax.dist=10
ax.azim=30
ax.elev=10
ax.set_xlim([0,1])
ax.set_ylim([0,1])
ax.set_zlim([0,1])
zAxis=0.0

for itr in range (20):
    zAxis=zAxis+0.05
    
    x=[0.0, zAxis]
    y=[0.0, 1-zAxis]
    z=[zAxis, 0.0 ]
        
    ax.plot(x, y, z, color="black")
    
#x = [0,0,0,1]
#y = [0,0,1,0]
#z = [0,1,0,0]
#
#
#ax.plot(x, y, z, 'o-')
colourinter=(1,0,0, 0.6)

 
#ax.scatter(0.2222,0.44444,0.111111,color=colourinter) Red Dot for battle of the sexes
ax.scatter(1,0,0,color=colourinter)
ax.scatter(0,0,0,color=colourinter) 

ax.text(-0.1,-0.1,0,  '%s' % (str('BR')), size=10, zorder=1,  color='k')


ax.text(0,0,1.02,  '%s' % (str('BL')), size=10, zorder=1,  color='k')


ax.text(0,1.02,0,  '%s' % (str('TR')), size=10, zorder=1,  color='k')


ax.text(1,-0.08,-0.05,  '%s' % (str('TL')), size=10, zorder=1,  color='k')
plt.show()