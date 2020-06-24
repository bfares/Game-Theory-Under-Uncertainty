import numpy as np

import pypoman
import matplotlib.pyplot as plt
from Polytope import Polytope


def vertices(A,b):
       
        vertices = pypoman.compute_polytope_vertices(A, b)

        return np.array(vertices)
    
def b_eq():
        B=[1]
        return np.array(B)

def A_eq():
        
        #sigma x =1
        constraints = [1 for i in range(3)]
        A=[constraints]
        return np.array(A)  
    
def c():
        
        #sigma x =1
        constraints = [1 for i in range(3)]
        A=[constraints]
        return np.array(A) 
    
    
ineq_A = np.array([
    [1,  0,  0],
    [0, -1,  0],
    [1,  0, -1],
    [0, 0,  1],
    
    [1, 0,  0],
    [0, 1,  0],
    [0, 0,  1],
    ])
b = np.array([0.5, -0.5, 0, 0.333,0,0,0])   



cPolytope=Polytope(c(), ineq_A, b, A_eq(),b_eq())

vertices1=cPolytope.vertices()

print ("Polytope Vertices: \n", vertices1)



