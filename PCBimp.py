# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 09:01:24 2014

@author: David Moorhouse
"""

import numpy as np
import matplotlib.pyplot as plt
import FDM as fdm

resolution = 4

#geometry='stripline'
geometry = 'edge coupled stripline'
trace_width = 200e-6
trace_spacing = 200e-6
L1_thickness = 35e-6
L1_to_L2_thickness = 500e-6
L2_thickness = 35e-6
L2_to_L3_thickness = 500e-6
L3_thickness = 35e-6

E0=(1/(36*np.pi))*1e-9  #Vacuum Permittivity
Er=4.8                   #Relative Permittivity
E=E0*Er                 #Permittivity
u0=4*np.pi*1e-7         #Vacuum Permeability
ur=1                    #Relative Permeability
u=u0*ur                 #Permeability
vp=1/(np.sqrt(u*E))     #Phase Velocity
ev=10                   #Electrode Voltage

cell_size = min(L1_thickness, L1_to_L2_thickness, L2_thickness,
                   L2_to_L3_thickness, L3_thickness)
Grow_min = np.round((L1_thickness + L1_to_L2_thickness + L2_thickness 
    + L2_to_L3_thickness+L3_thickness)/cell_size)
Grow = Grow_min*resolution             #Resolution of geometry matrix in Rows
trace_width_cells = np.round(resolution*trace_width/cell_size)
trace_spacing_cells = np.round(resolution*trace_width/cell_size)
Gcol = np.round((trace_width * 5)/cell_size)*resolution
if np.mod(Gcol,trace_width_cells)!=0:
    Gcol+=1
if np.mod(Grow,2)==0:
    Grow+=1

G1 = fdm.geometry_matrix(geometry,Grow,Gcol,trace_width_cells,trace_spacing_cells)
Q, nodes, A,v = fdm.process_matrix(G1,ev)
q,C,Z0,L = fdm.trace_values(Q,ev,E,u)


#Print results of calculations
print('Resolution of matrix = %i Rows')%(Grow);
print('Quantity of Nodes = %i Nodes')%(nodes);
print('Capacitance per unit length = %.2f pF')%(C*1e12);
print('Inductance per unit length = %.2f nH')%(L*1e9);
print('Characteristic Impedance = %.2f Ohms')%(Z0);
#Display contour plot
plt.figure(1);
a=np.r_[0:10.1:0.1]
plt.contourf(Q,a)
title= geometry+' '+(str(int(nodes)))+' Nodes'
plt.title(title)
plt.xlabel('Matrix Column');
plt.ylabel('Matrix Row')
plt.xlim((0,Gcol-1))
plt.ylim((0,Grow-1))
plt.show()

#    #Display capacitance plot
#    plt.figure(t+1);
#    plt.semilogx(nodes,CpF,'--rs','LineWidth',2,...
#     'MarkerEdgeColor','k',...
#     'MarkerFaceColor','g',...
#     'MarkerSize',10);
#    plt.grid('on');
#    plt.title('Capacitance Change with Calculation Precision');
#    plt.xlabel('Quadrant Nodes');
#    plt.ylabel('Capacitance pf per unit length');
#    #Display inductance plot
#    plt.figure(t+2);
#    plt.semilogx(nodes,LnH,'--rs','LineWidth',2,...
#     'MarkerEdgeColor','k',...
#     'MarkerFaceColor','g',...
#     'MarkerSize',10);
#    plt.grid;
#    plt.title('Inductance Change with Calculation Precision');
#    plt.xlabel('Quadrant Nodes');
#    plt.ylabel('Inductance nH per unit length');
#    #Display characteristic impedance plot
#    plt.figure(t+3);
#    plt.semilogx(nodes,Z0,'--rs','LineWidth',2,...
#     'MarkerEdgeColor','k',...
#     'MarkerFaceColor','g',...
#     'MarkerSize',10);
#    plt.grid;
#    plt.title('Characteristic Impedance Change with Calculation Precision');9
#    plt.xlabel('Quadrant Nodes');
#    plt.ylabel('Characteristic Impedance');