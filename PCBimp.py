# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 09:01:24 2014

@author: david.moorhouse
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
import FDM as fdm

nodes=20
height=np.sqrt(nodes)
width=np.sqrt(nodes)
cross_sec=np.zeros((width,height))

cross_sec[np.round(width/2),np.round(height/2)]=1
cross_sec[:,0]=-1
cross_sec[:,width-1]=-1
cross_sec[0,:]=-1
cross_sec[width-1,:]=-1

#Geometry Matrix
#Full Height of Geometry
H=2;
#Full Width of Geometry
W=6;
#Full Conductor Height
Ch=1;
#Vacuum Permittivity
E0=(1/(36*np.pi))*1e-9;
#Relative Permittivity
Er=1;
#Permittivity
E=E0*Er;
#Vacuum Permeability
u0=4*np.pi*1e-7;
#Relative Permeability
ur=1;
#Permeability
u=u0*ur;
#Phase Velocity2
vp=1/(np.sqrt(u*E));
#Electrode Voltage
ev=10;
#Resolution of Quadrant in Rows
#GrowMat=[10, 26, 50, 100, 250, 500, 1000];
GrowMat=[20]
C=np.zeros((1,np.size(GrowMat)))
CpF=np.zeros((1,np.size(GrowMat)))
L=np.zeros((1,np.size(GrowMat)))
LnH=np.zeros((1,np.size(GrowMat)))
L2=np.zeros((1,np.size(GrowMat)))
Z0=np.zeros((1,np.size(GrowMat)))
nodes=np.zeros((1,np.size(GrowMat)))

for t in range(0,np.size(GrowMat)):
    #Resolution of quadrant in rows
    Grow=GrowMat[t];
    # Set column np.size from aspect ratio
    Gcol=(W/2)*Grow;
    #Initialise geometry matrix
    G1=np.zeros((Grow,Gcol))
    # Conductor height in quadrant in terms of rows
    Cq=2*Grow*Ch/H;
    #Matrix G2
#    G2=np.zeros((Grow,Gcol))

#    #Set position of conductor in complete matrix (4 quadrants)
#    Qrow=(Grow*2)-1;
#    Qcol=(Gcol*2)-1;
    #Populate geometry matrix for conductor and ground
    G1[0,:] = -1
    G1[Grow-1,:] = -1
    G1[:,0] = -1  
    G1[(Grow/2)-(Cq/2):(Grow/2)+(Cq/2),(Gcol-1)/2] = -2
    #Populate geometry matrix for nodes
    
    Q, node =fdm.process_matrix(G1,ev)
    
#    loop_iter=1;
    #Copy top left quadrant voltages across horizontal mirror
#    Q[Grow:,:Gcol]=(Q[:Grow-1,:Gcol])[::-1,:]
    #Copy left half voltages across vertical mirror
#    Q[:,Gcol-1:]=(Q[:,:Gcol])[:,::-1]
#    for r in range(Grow,Qrow):
#        for c in range(0,Gcol):
#            Q[r,c]=Q[Grow-1-loop_iter,c];
#    
#        loop_iter=loop_iter+1;
#    
#    loop_iter=1;
#    #Copy left half voltages across vertical mirror
#    for r in range(0,Qrow):
#        for c in range(Gcol+1,Qcol):
#            Q[r,c]=Q[r,Gcol-loop_iter];
#            loop_iter=loop_iter+1;
#        loop_iter=1;
    
    #used to verify Von and Vin summations
    #D=np.zeros(Qrow,Qcol);
    Von=0;
    #Perform column summation of outer "Von" node voltages neglecting corners
    for c in range(3,Gcol-2+1):
        for r in range(2,Grow-1+1,Grow-3):
     #Set values adjacent to corners to 0.5
            if c==3 or c==Gcol-2:
                Von=Von+(Q[r,c]/2);
                #D[r,c]=0.5;
            else:
                Von=Von+Q[r,c];
                #D[r,c]=1;
    
    #Perform row summation of outer "Von" node voltages neglecting corners
    for r in range(3,Grow-2+1):
        for c in range(2,Gcol-1+1,Gcol-3):
    #Set values adjacent to corners to 0.5
            if r==3 or r==Grow-2:
                Von=Von+(Q[r,c]/2);
                #D[r,c]=0.5;
            else:
                Von=Von+Q[r,c];
                #D[r,c]=1;
    Vin=0;
    #Perform column summation of inner "Vin" node voltages neglecting corners
    for c in range(4,Gcol-3+1):
        for r in range(3,Grow-2+1,Grow-5):
    #Set values adjacent to corners to 0.5
            if c==4 or c==Gcol-3:
                Vin=Vin+(Q[r,c]/2);
                #D[r,c]=-0.5;
            else:
                Vin=Vin+Q[r,c];
                #D[r,c]=-1;
    #Perform row summation of inner "Vin" node voltages neglecting corners
    for r in range(4,Grow-3+1):
        for c in range(3,Gcol-2+1,Gcol-5):
    #Set values adjacent to corners to 0.5
            if r==4 or r==Grow-3:
                Vin=Vin+(Q[r,c]/2);
                #D[r,c]=-0.5;
            else:
                Vin=Vin+Q[r,c];
                #D[r,c]=-1;
    #Calculate charge
    q=E*(Vin-Von);
    #Calculate capacitance
    C[t]=q/ev;
    #Scale capacitance to pico farads
    CpF[t]=C[t]*1e12;
    #Calculate characteristic impedance
    Z0[t]=(np.sqrt(u*E))/C[t];
    #Calculate inductance
#    L[t]=(1/(vp*np.sqrt(C[t])))^2
    #Calculate inductance using a different method to check above value
    L2[t]=C[t]*Z0[t]**2;
    #Scale inductance to nano henrys
    LnH[t]=L[t]*1e9;
    #Add current node value to node matrix
    nodes[t]=node;
    #Print results of calculations
    print('Resolution of 1 Quadrant\n= %i Rows\n\n')%(Grow);
    print('Quantity of Quadrant Nodes\n= %i Nodes\n\n')%(nodes[t]);
    print('Capacitance per unit length\n= %.2f pF\n\n')%(CpF[t]);
    print('Inductance per unit length\n= %.2f nH\n\n')%(LnH[t]);
    print('Characteristic Impedance\n= %.2f Ohms\n\n\n\n')%(Z0[t]);
    #Display contour plot
    plt.figure(t);
    a=[1,2,3,4,5,6,7,8,9,10];
    plt.contourf(Q,a);
    plt.title(['Equipotentials at 1 Volt Interval ', str(nodes[t]),' Quadrant Nodes']);
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