# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 09:01:24 2014

@author: David Moorhouse
"""

import numpy as np
import matplotlib.pyplot as plt
import FDM as fdm
from scipy.constants import epsilon_0 as E0
from scipy.constants import mu_0 as u0
from scipy.constants import c

resolutions = np.r_[1:12]

#geometry='stripline'
geometry = 'edge coupled stripline'
trace_width = 200e-6
trace_spacing = 150e-6
L1_thickness = 35e-6
L1_to_L2_thickness = 200e-6
trace_thickness = 35e-6
L2_to_L3_thickness = 200e-6
L3_thickness = 35e-6

#E0=(1/(36*np.pi))*1e-9  #Vacuum Permittivity
Er=4.8                   #Relative Permittivity
E=E0*Er                 #Permittivity
#u0=4*np.pi*1e-7         #Vacuum Permeability
ur=1                    #Relative Permeability
u=u0*ur                 #Permeability
vp=1/(np.sqrt(u*E))     #Phase Velocity
ev=5                   #Electrode Voltage

cell_size = min(L1_thickness, L1_to_L2_thickness, trace_thickness,
                   L2_to_L3_thickness, L3_thickness)
Grow_min = np.round((L1_thickness + L1_to_L2_thickness + trace_thickness 
    + L2_to_L3_thickness+L3_thickness)/cell_size)
i=0
for resolution in resolutions:
    Grow = Grow_min*resolution             #Resolution of geometry matrix in Rows
    trace_width_cells = np.round(resolution*trace_width/cell_size)
    trace_thickness_cells = np.round(resolution*trace_thickness/cell_size)
    trace_spacing_cells = np.round(resolution*trace_width/cell_size)
    Gcol = np.round((trace_width * 5)/cell_size)*resolution
    if np.mod(Gcol,trace_width_cells)!=0:
        Gcol+=1
    if np.mod(Grow,2)==0:
        Grow+=1
    
    G1,traces = fdm.geometry_matrix(geometry,Grow,Gcol,trace_width_cells,
                                    trace_spacing_cells,trace_thickness_cells)
    Q, nodes, A,v = fdm.process_matrix(G1,ev)
    q,C,Z0,L = fdm.trace_values(Q,ev,E,u,vp,Er,c,traces)
    
    
    #Print results of calculations
#    print('Resolution of matrix = %i Rows')%(Grow);
#    print('Quantity of Nodes = %i Nodes')%(nodes);
#    print('Capacitance per unit length = %.2f pF')%(C*1e12);
#    print('Inductance per unit length = %.2f nH')%(L*1e9);
    print('Characteristic Impedance = %.2f Ohms')%(Z0);
    #Display contour plot
    plt.figure(1)
    plt.ion()
    a=np.r_[-ev:ev+0.1:0.2]
    plt.contourf(Q,a)
    title= geometry+' '+(str(int(nodes)))+' Nodes'
    plt.title(title)
    plt.xlabel('Matrix Column');
    plt.ylabel('Matrix Row')
    plt.xlim((0,Gcol-1))
    plt.ylim((0,Grow-1))
    plt.pause(0.001)
    plt.show()
    
    if i==0:
        impedances=[Z0]
    else:
        impedances=np.append(impedances,Z0)
        if np.abs(impedances[i-1]-impedances[i])<impedances[i]*0.01:
            break
    i+=1