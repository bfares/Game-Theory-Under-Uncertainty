# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 19:55:30 2019

@author: Bernard Fares
"""

#from scipy.spatial import HalfspaceIntersection
from scipy.spatial import ConvexHull
#from mpl_toolkits.mplot3d import Axes3D
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy import optimize
#from numpy.linalg import solve
#import pylab
import scipy as sp
import numpy as np
import pypoman
import pylab
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
from sklearn.decomposition import PCA
#import cdd
class Polytope(object):
    A_eq=""
    b_eq=""
    A_ineq=""
    b_ineq=""
    c=""
    ax=""
    colour=""
    def __init__(self,c, A_ineq, b_ineq, A_eq, b_eq, ax=""):
        self.A_eq= A_eq
        self.b_eq= b_eq
        self.A_ineq= A_ineq
        self.b_ineq= b_ineq
        self.c= c
        
        if ax=="":
            self.ax= a3.Axes3D(plt.figure())
        else:
            self.ax=ax
     
    def setFigure(self, ax):
        self.ax=ax
        
    def getFigure(self):
        return self.ax
    
    def minimize(self):
        A_eq=self.A_eq 
        b_eq=self.b_eq 
        A_ineq =-1* self.A_ineq 
        b_ineq=-1*self.b_ineq 
        c=-1*self.c
        
        res=optimize.linprog(c=c, A_ub=A_ineq, b_ub=b_ineq, A_eq=A_eq, b_eq=b_eq, method='interior-point', bounds=(0, 1) )
        #res=optimize.linprog(c,A_ineq,b_ineq, method="interior-point" )
        if res.success:
            print("LP Solved")
        else:
            print("LP failed:", res.message)
        #print('Optimal value:', res.fun, '\nX:', res.x)
        return res.x
        
    def vertices(self):
        A_ineq = -1 * self.A_ineq 
        b_ineq = -1 * self.b_ineq 
        A_eq= self.A_eq
        b_eq= self.b_eq
        
        A= np.vstack((A_ineq, A_eq,-1*A_eq))
        b= np.concatenate((b_ineq, b_eq,-1*b_eq))
        #At=np.array([[-2,1,0,0], [0,0,2,-1],[-1,0,2,0],[0,1,0,-2],[-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1], [1,1,1,1], [-1,-1,-1,-1]])
        #bt = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1,-1])
        
        vertices = pypoman.compute_polytope_vertices(A, b)

        return np.array(vertices)
        
    def plotConvexHull(self, colour=(sp.rand(),sp.rand(),sp.rand(), 0.1)):
        A = -1 * self.A_ineq
        b = -1 * self.b_ineq
        C = self.A_eq
        d = self.b_eq
        ineq = (A, b)  # A * x <= b
        eq = (C, d)    # C * x == d
        n = len(A[0])  # dimension of the original polytope
        p = 3   # dimension of the projected polytope
         
       
        # Projection is proj(x) = [x_0 x_1]
        E = np.zeros((p, n))
        E[0, 0] = 1.
        E[1, 1] = 1.
        E[2, 2] = 1.
       
        f = np.zeros(p)
        proj = (E, f)  # proj(x) = E * x + f
        #
        vertices = pypoman.project_polytope(proj, ineq, eq, method='cdd')
        verts= np.around(np.array(vertices),4)
        
                
        #verts = np.around(self.vertices(), 4) #np.matrix([[ 0.4, 0.0, 0.2, 0.4], [1.0, 0.0, 0.0, 0.0] , [ 0.25, 0.5, 0.0, 0.25], [ 0.222, 0.444, 0.111, 0.222],[ 0.0, 0.0, 0.0, 1.0]])
        hull = ConvexHull(verts, qhull_options="QJ")
        faces = hull.simplices
        
        ax = self.ax
        ax.dist=10
        ax.azim=30
        ax.elev=10
        ax.set_xlim([0,1])
        ax.set_ylim([0,1])
        ax.set_zlim([0,1])
        for s in faces:
            sq = [
                [verts[s[0], 0], verts[s[0], 1], verts[s[0], 2]],
                [verts[s[1], 0], verts[s[1], 1], verts[s[1], 2]],
                [verts[s[2], 0], verts[s[2], 1], verts[s[2], 2]]
            ]
            f = a3.art3d.Poly3DCollection([sq])#, alpha=0.1, linewidths=1)
            f.set_color(colors.rgb2hex(sp.rand(3)))
            f.set_edgecolor('k')
            f.set_facecolor(colour)
            
            ax.add_collection3d(f)
            #ax.set_axis_off()
        plt.show()


    def projected3DVertices(self):
            A = -1 * self.A_ineq
            b = self.b_ineq
            C = self.A_eq
            d = self.b_eq
            ineq = (A, b)  # A * x <= b
            eq = (C, d)    # C * x == d
            n = len(A[0])  # dimension of the original polytope
            p = 3   # dimension of the projected polytope
             
           
            # Projection is proj(x) = [x_0 x_1]
            E = np.zeros((p, n))
            E[0, 0] = 1.
            E[1, 1] = 1.
            E[2, 2] = 1.
           
            f = np.zeros(p)
            proj = (E, f)  # proj(x) = E * x + f
            #
            vertices = pypoman.project_polytope(proj, ineq, eq, method='cdd')
            verts= np.around(np.array(vertices),4)
                    
            return verts
    
    def plotConvexHull_NoProj(self):
        
        verts= np.around(np.array(self.vertices()),4)
        hull = ConvexHull(verts, qhull_options="Qt")
        faces = hull.simplices
        
        ax = self.ax
        ax.dist=10
        ax.azim=30
        ax.elev=10
        ax.set_xlim([0,1])
        ax.set_ylim([0,1])
        ax.set_zlim([0,1])
        for s in faces:
            sq = [
                [verts[s[0], 0], verts[s[0], 1], verts[s[0], 2]],
                [verts[s[1], 0], verts[s[1], 1], verts[s[1], 2]],
                [verts[s[2], 0], verts[s[2], 1], verts[s[2], 2]]
            ]
            f = a3.art3d.Poly3DCollection([sq])#, alpha=0.1, linewidths=1)
            f.set_color(colors.rgb2hex(sp.rand(3)))
            f.set_edgecolor('k')
            
            f.set_facecolor(self.colour)
           
            ax.add_collection3d(f)
            #ax.set_axis_off()
        plt.show()

    def plotConvexHull_CS_test(self):
            
            verts= np.around(np.array(self.vertices()),4)

            hull = ConvexHull(verts)
            faces = hull.simplices
            
            ax = self.ax
            ax.dist=10
            ax.azim=30
            ax.elev=10
            ax.set_xlim([0,1])
            ax.set_ylim([0,1])
            ax.set_zlim([0,1])
            for s in faces:
                sq = [
                    [verts[s[0], 0], verts[s[0], 1], verts[s[0], 2]],
                    [verts[s[1], 0], verts[s[1], 1], verts[s[1], 2]],
                    [verts[s[2], 0], verts[s[2], 1], verts[s[2], 2]]
                ]
                f = a3.art3d.Poly3DCollection([sq])#, alpha=0.1, linewidths=1)
                f.set_color(colors.rgb2hex(sp.rand(3)))
                f.set_edgecolor('k')
                
                f.set_facecolor(self.colour)
               
                ax.add_collection3d(f)
                #ax.set_axis_off()
            plt.show()      
            
    def projectedNVertices(self,p):
            # dimension of the projected polytope
            A = -1 * self.A_ineq
            b = -1 * self.b_ineq
            C = self.A_eq
            d = self.b_eq
            ineq = (A, b)  # A * x <= b
            eq = (C, d)    # C * x == d
            n = len(A[0])  # dimension of the original polytope
             # dimension of the projected polytope
             
           
            # Projection is proj(x) = [x_0 x_1]
            E = np.zeros((p, n))
            E[0, 0] = 1.
            E[1, 1] = 1.
            f = np.zeros(p)
            proj = (E, f)  # proj(x) = E * x + f
            #
            vertices = pypoman.project_polytope(proj, ineq, eq, method='bretl')
            verts= np.asmatrix(vertices)
            hull = ConvexHull(vertices, qhull_options="QJ")
            faces = hull.simplices
            
            print ("verts:") 
            print (verts)
            
            print ("faces:") 
            print (faces)
            
            ax = self.ax
            ax.dist=10
            ax.azim=30
            ax.elev=10
            ax.set_xlim([0,1])
            ax.set_ylim([0,1])
            ax.set_zlim([0,1])
            for s in faces:
                sq = [
                    [verts[s[0], 0], verts[s[0], 1]],
                    [verts[s[1], 0], verts[s[1], 1]]
                    
                ]
                plt.plot([sq])
                #f = a3.art3d.Poly3DCollection([sq])#, alpha=0.1, linewidths=1)
                #f.set_color(colors.rgb2hex(sp.rand(3)))
                #f.set_edgecolor('k')
                #f.set_facecolor('Red')
                
                #ax.add_collection3d(f)
                #ax.set_axis_off()
            plt.show()
                    
            return verts
        
    def minimize_test(self):
        
        A_eq=self.A_eq 
        
        b_eq=self.b_eq 
        A_ineq=  [[2,0,1,0], [-1,0,0,-1], [0,-2,-2,0], [0,1,0,2], [-1,0,0,0], [0,-1,0,0], [0,0,-1,0], [0,0,0,-1]] 
        b_ineq=  [[-1,-1,-1,-1,0,0,0,0]] 
        c= [[0,0,0,0]]
        
        res=optimize.linprog(c=c, A_ub=A_ineq, b_ub=b_ineq, A_eq=A_eq, b_eq=b_eq, method='interior-point', bounds=(0, 1) )
        #res=optimize.linprog(c,A_ineq,b_ineq, method="interior-point" )
        if res.success:
            print("LP Solved")
        else:
            print("LP failed:", res.message)
        #print('Optimal value:', res.fun, '\nX:', res.x)
        return res.x