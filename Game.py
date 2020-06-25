# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 07:50:49 2019

@author: Bernard Fares
"""
import numpy as np
import math
from collections import OrderedDict
class Game(object ):
    player1PayoffMatrix=""
    player2PayoffMatrix=""
    LowerPrevisionsDict=""
    UpperPrevisionsDict=""
    rows=0
    cols=0
    def __init__(self, player1PayoffMatrix, player2PayoffMatrix, LowerPrevisionsDict=None, UpperPrevisionsDict=None):
        self.player1PayoffMatrix=player1PayoffMatrix
        self.player2PayoffMatrix=player2PayoffMatrix
        self.LowerPrevisionsDict=LowerPrevisionsDict
        self.UpperPrevisionsDict=UpperPrevisionsDict
        self.rows= int(len(player1PayoffMatrix))
        self.cols= int(len(player1PayoffMatrix[0])) 
        
    
    def displayPayoff(self):
        player1PayoffMatrix= self.player1PayoffMatrix
        player2PayoffMatrix= self.player2PayoffMatrix
        m=self.rows
        n=self.cols
        
        for row in range (m):
            for column in range(n):
                print(player1PayoffMatrix[row][column],",",player2PayoffMatrix[row][column],"    ", end='')
            print(" ")
             
    def b_ineq(self):   
        return np.concatenate((self.CE_b_ineq(), self.identity_b_ineq()), axis=None)
    
    def identity_b_ineq(self):
        m=self.rows
        n=self.cols
        return np.zeros(m*n)
    
    def CE_b_ineq(self):
        m=self.rows
        n=self.cols
        return np.zeros(math.factorial(m) + math.factorial(n))
            
    def A_eq(self):
        m=self.rows
        n=self.cols
        return np.ones((1,m*n),dtype=int)
    
    def b_eq(self):
        return np.array([1])
    
    def initConstraints(self):
        m=self.rows
        n=self.cols
        constraints=OrderedDict()
        for row in range (m):
            for column in range(n):
              constraints[str(row)+str(column)]=0
        return constraints

    def A_ineq_RiskAverse(self, P1_MarginalUtilitiesMatrix, P2_MarginalUtilitiesMatrix):
        player1PayoffMatrix= self.player1PayoffMatrix
        player2PayoffMatrix= self.player2PayoffMatrix
        m=self.rows
        n=self.cols
        A=[]
        
        #Player 1 constraints
        for row in range (m): 
            for q in range (m):
                if(row==q):
                     continue
                constraints=self.initConstraints()
                for column in range(n):
                    constraints[str(row)+str(column)]= (player1PayoffMatrix[row][column] - player1PayoffMatrix[q][column])/P1_MarginalUtilitiesMatrix[row][column]            
                
                A+=[list(constraints.values())]
        
        #Player 2 constraints
        for column in range (m): 
            for k in range (m):
                if(column==k):
                     continue
                constraints=self.initConstraints()
                for row in range(n):
                    constraints[str(row)+str(column)]= (player2PayoffMatrix[row][column] - player2PayoffMatrix[row][k] )/P2_MarginalUtilitiesMatrix[row][column]            
                A+=[list(constraints.values())] 
        #x>=0 constraints
        for row in range (m):    
            for column in range(n):
                constraints=self.initConstraints()
                constraints[str(row)+str(column)]= 1
                A+=[list(constraints.values())]
        return np.array(A)
    
    def A_ineq(self):     
        return np.concatenate((self.CE_A_ineq(), self.identity_A_ineq()), axis=0)
    
    def identity_A_ineq(self):
        m=self.rows
        n=self.cols
        return np.identity(m*n)
    
    def CE_A_ineq(self):
        player1PayoffMatrix= self.player1PayoffMatrix
        player2PayoffMatrix= self.player2PayoffMatrix
        lP=self.LowerPrevisionsDict
        uP=self.UpperPrevisionsDict
        m=self.rows
        n=self.cols
        A=[]
        
        #Player 1 constraints
        for row in range (m): 
            for q in range (m):
                if(row==q):
                     continue
                constraints=self.initConstraints()
                for column in range(n):
                    if not self.isNumber(player1PayoffMatrix[row][column]):                        
                        if player1PayoffMatrix[row][column].startswith('-'):
                            player1Payoff= -uP[player1PayoffMatrix[row][column]]
                        else:
                            player1Payoff=lP[player1PayoffMatrix[row][column]]
                    else:
                        player1Payoff=player1PayoffMatrix[row][column]
                    
                    if not self.isNumber(player1PayoffMatrix[q][column]):
                         if player1PayoffMatrix[q][column].startswith('-'):
                            player1Payoffq= -lP[player1PayoffMatrix[q][column]]
                         else:
                            player1Payoffq=uP[player1PayoffMatrix[q][column]]
                            
                    else:
                        player1Payoffq=player1PayoffMatrix[q][column]
                   
                    constraints[str(row)+str(column)]= float(player1Payoff) - float(player1Payoffq)             
                    
                A+=[list(constraints.values())]
        
        #Player 2 constraints
        for column in range (m): 
            for k in range (m):
                if(column==k):
                     continue
                constraints=self.initConstraints()
                for row in range(n):
                     if not self.isNumber(player2PayoffMatrix[row][column]):  
                        
                        if player2PayoffMatrix[row][column].startswith('-'):
                            player2Payoff= -uP[player2PayoffMatrix[row][column]]
                        else:
                            player2Payoff=lP[player2PayoffMatrix[row][column]]
                     else:
                        player2Payoff=player2PayoffMatrix[row][column]
                    
                     if not self.isNumber(player2PayoffMatrix[row][k]):
                         if player2PayoffMatrix[row][k].startswith('-'):
                            player2Payoffk= -lP[player2PayoffMatrix[row][k]]
                         else:
                            player2Payoffk=uP[player2PayoffMatrix[row][k]]
                     else:
                        player2Payoffk=player2PayoffMatrix[row][k]
                        
                     constraints[str(row)+str(column)]= float(player2Payoff) - float(player2Payoffk)             
                A+=[list(constraints.values())] 
        return np.array(A)
    
    def isNumber(self, input):
        isNumeric=True
             
        try:
             float(input)
        except ValueError:
            isNumeric=False
                
        return isNumeric
    
    def c(self):
        #commonKnowledgeMatrix= self.A_ineq_gamble( lP, uP)     
        #print (commonKnowledgeMatrix)
        #return np.sum(commonKnowledgeMatrix, axis=0)
        lP=self.LowerPrevisionsDict
        uP=self.UpperPrevisionsDict
        player1PayoffMatrix= self.player1PayoffMatrix
        player2PayoffMatrix= self.player2PayoffMatrix
        m=self.rows
        n=self.cols
        c=[]
        # maximize matrix c
        c=self.initConstraints()
        for row in range (m): 
            for column in range (n):
                if not self.isNumber(player1PayoffMatrix[row][column]):                        
                        if player1PayoffMatrix[row][column].startswith('-'):
                            player1Payoff= -uP[player1PayoffMatrix[row][column]]
                        else:
                            player1Payoff=lP[player1PayoffMatrix[row][column]]
                else:
                        player1Payoff=player1PayoffMatrix[row][column]
                        
                if not self.isNumber(player2PayoffMatrix[row][column]):                        
                        if player2PayoffMatrix[row][column].startswith('-'):
                            player2Payoff= -uP[player2PayoffMatrix[row][column]]
                        else:
                            player2Payoff=lP[player2PayoffMatrix[row][column]]
                else:
                        player2Payoff=player2PayoffMatrix[row][column]
                        
                c[str(row)+str(column)] = float(player1Payoff) + float(player2Payoff) # sum of payoffs for both players
        
        return np.array(list(c.values()))
    
    
    
    
    
   