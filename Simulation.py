# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 22:49:34 2019

@author: Bernard Fares
"""
from Polytope import Polytope
from Game import Game
import numpy as np
import matplotlib.pyplot as plt
from Minimize import Minimize
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

Thesis_Test_p1 = np.array( [[2,0], [0,1]])
Thesis_Test_p2 = np.array( [[1,0], [0,2]])

lP={'f':1.33}
uP={'f':1.83}

#cGame=Game(chicken_p1, chicken_p2)
#cGame=Game(shapley_p1, shapley_p2)
#cGame=Game(battleOfTheSexes_p1, battleOfTheSexes_p2)
cGame=Game(Thesis_Test_p1, Thesis_Test_p2, lP, uP)


#ineq_A=cGame.A_ineq_RiskAverse(matchingPennies_p1_MarginalUtilities,matchingPennies_p2_MarginalUtilities)

c=cGame.c()
A_ineq= cGame.A_ineq()
b_ineq= cGame.b_ineq()
A_eq= cGame.A_eq()
b_eq= cGame.b_eq()
identity_A_ineq= cGame.identity_A_ineq()
identity_b_ineq= cGame.identity_b_ineq()
print("Game CE constraints: \n", cGame.CE_A_ineq())
print ("\n Game objective function: \n",c)


minimizeLP=Minimize(-1*c, -1*A_ineq, -1*b_ineq , A_eq,b_eq)
CEequi=minimizeLP.compute()
print ("\n Game Correlated Equilibria: \n", CEequi)

cPolytope=Polytope(c, A_ineq, b_ineq , A_eq,b_eq)
vertices=cPolytope.vertices()
print ("\n Polytope Vertices: \n", vertices)

polytopeColour=(0.6,1,0.6, 0.1)
cPolytope.plotConvexHull(polytopeColour)

cSimplex=Polytope(c, identity_A_ineq, identity_b_ineq, cGame.A_eq(),cGame.b_eq(), cPolytope.getFigure())
simplexColour=(0.878,0.878,0.878, 0.05)
cSimplex.plotConvexHull(simplexColour)

#Plot labels and dots
ax = cPolytope.getFigure()

zAxis=0.0
for itr in range (20):
    zAxis=zAxis+0.05
    
    x=[0.0, zAxis]
    y=[0.0, 1-zAxis]
    z=[zAxis, 0.0 ]
        
    ax.plot(x, y, z, color="black")
    
colourinter=(1,0,0, 0.6)
#ax.scatter(0.2222,0.44444,0.111111,color=colourinter) Red Dot for battle of the sexes
ax.scatter(1,0,0,color=colourinter)
ax.scatter(0,0,0,color=colourinter) 

ax.text(-0.1,-0.1,0,  '%s' % (str('BR')), size=10, zorder=1,  color='k')
ax.text(0,0,1.02,  '%s' % (str('BL')), size=10, zorder=1,  color='k')
ax.text(0,1.02,0,  '%s' % (str('TR')), size=10, zorder=1,  color='k')
ax.text(1,-0.08,-0.05,  '%s' % (str('TL')), size=10, zorder=1,  color='k')
plt.show()