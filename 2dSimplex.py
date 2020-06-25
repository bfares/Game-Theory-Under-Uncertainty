import ternary
import numpy as np

scale=1
def plotSimplex( data, labels):
    
    scale = 1
    figure, tax = ternary.figure(scale=scale)
    
    # Draw Boundary and Gridlines
    tax.boundary(linewidth=2.0)
    tax.gridlines(color="blue", multiple=0.1)
    
    # Set Axis labels and Title
    fontsize = 12
    #offset = 0.14
    
    tax.right_corner_label("G", fontsize=fontsize)
    tax.top_corner_label("B", fontsize=fontsize)
    tax.left_corner_label("N", fontsize=fontsize)
    
    # Draw an arbitrary line, ternary will project the points for you
    
    tax.left_parallel_line(0.5, linewidth=1., color='darkblue', linestyle="-")
    tax.right_parallel_line(0.333, linewidth=1., color='darkblue', linestyle="-")
    tax.horizontal_line(0.5, linewidth=1., color='darkblue', linestyle="-")
    tax.line([0.5,0,0.5],[0,1,0], linewidth=1., color='darkblue', linestyle="-")
    
    data= data[np.lexsort((data[:,0], data[:,1]))][::-1]        #reverse
    
    for p in range (len(data)):
       if p== len(data)-1:
        p1= [data[p][0],data[p][1],data[p][2]]
        p2=[data[0][0],data[0][1],data[0][2]]
        
       else:  
        p1= [data[p][0],data[p][1],data[p][2]]
        p2=[data[p+1][0],data[p+1][1],data[p+1][2]] 
       tax.line(p1, p2, linewidth=1.,  color='darkblue', linestyle="-")
       
    labelColor= 'darkblue'
    
    
    ternary.plt.annotate("",
            xy=(0.670, 0.180), xycoords='data',
            xytext=(0.622, 0.22), textcoords='data',
            arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", color=labelColor),
            )
    ternary.plt.text(0.67, 0.22, '$d_1$', size=9, ha='center', va='center', color=labelColor)
    
   
    ternary.plt.annotate("",
            xy=(0.3, 0.37), xycoords='data',
            xytext=(0.3, 0.44), textcoords='data',
            arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", color=labelColor),
            )
    ternary.plt.text(0.33, 0.4, '$d_2$', size=9, ha='center', va='center', color=labelColor)
    
    ternary.plt.annotate("",
            xy=(0.55, 0.5), xycoords='data',
            xytext=(0.5, 0.5), textcoords='data',
            arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", color=labelColor),
            )
    
    ternary.plt.text(0.55, 0.52, '$d_3$', size=9, ha='center', va='center', color=labelColor)
    
    ternary.plt.annotate("",
            xy=(0.43, 0.31), xycoords='data',
            xytext=(0.47, 0.35), textcoords='data',
            arrowprops=dict(arrowstyle="-|>", connectionstyle="arc3", color=labelColor),
            )
    ternary.plt.text(0.43, 0.36, '$d_4$', size=9, ha='center', va='center', color=labelColor)
#    p1 = (1, 0, 0)
#    p2 = (0.8, 0, 0.2)
#    tax.line(p1, p2, linewidth=1.,  color='blue', linestyle=":")
#    
#    p1 = (0.8, 0, 0.2)
#    p2= (0, 0.44, 0.55)
#    
#    tax.line(p1, p2, linewidth=1., marker='s', color='blue', linestyle=":")
#    
#    p1 = (0, 0.44, 0.55)
#    p2 = (0,0.5,0.5)
#    tax.line(p1, p2, linewidth=1., marker='s', color='blue', linestyle=":")
#    
#    p1 = (0,0.5,0.5)
#    p2 = (1,0,0)
#    
#    tax.line(p1, p2, linewidth=1., marker='s', color='blue', linestyle=":")
    
    #tax.ticks(axis='lbr', multiple=5, linewidth=1, offset=0.025)
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()
    tax.show()

    
data = np.array(

[[0.667, 0,    0.333],
 [0.5,   0.167, 0.333],
 [0.5,   0,    0.5  ]])

labels=['X', 'Y', 'Z']

plotSimplex(data, labels)

    