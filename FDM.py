# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 20:50:53 2014

@author: David Moorhouse
"""

import numpy as np
from scipy import sparse as sp
from scipy.sparse.linalg import spsolve


def process_matrix(G1,ev):
    #initialise node numbering
    node=-1; 
    Grow=G1.shape[0]
    Gcol=G1.shape[1]
    for r in range(0,Grow):
        for c in range(0,Gcol):
            if G1[r,c]==0:  
                node+=1
                G1[r,c]=node
                
    
    #Define matrix sizes for node voltage equation matricies
    h=(node*4)-(Grow+Gcol)-1;
    Ar=np.zeros((h))
    Ac=np.zeros((h))
    Ad=np.zeros((h))
    Ainc=1;
    b=np.zeros((node+1))
    #Process geometry matrix cell by cell
    for r in range(0,Grow):
        for c in range(0,Gcol):
            y=c;
    #Process cells 1 row below and 1 row above current cell
            for x in [r-1,r+1]:
                if (x>0 and y>0 and G1[r,c]>0 and x<Grow and y<=Gcol):
                    e=G1[r,c];
                    f=G1[x,y];
    #Process cells that fall on the conductor
                    if G1[x,y] == -2:
                        b[(G1[r,c])]=b[(G1[r,c])]+ev
                    elif G1[x,y] == -3:
                        b[(G1[r,c])]=b[(G1[r,c])]-ev
    #Do nothing with ground cells
                    elif G1[x,y] == -1:
                        pass
    #Process nodes that are in the final row only
                    else:                      
                        if r==Grow:
                            Ar[Ainc]=e;
                            Ac[Ainc]=f;
                            Ad[Ainc]=Ad[Ainc]-1;

                        #Process nodes
                        Ar[Ainc]=e;
                        Ac[Ainc]=f;
                        Ad[Ainc]=Ad[Ainc]-1;
                        Ainc=Ainc+1;
            x=r;
    #Process cells 1 column left and 1 column right of the current cell
            for y in [c-1,c+1]:
                if (x>0 and y>0 and G1[r,c]>0 and x<Grow and y<Gcol):
                    e=G1[r,c];
                    f=G1[x,y];
        #Process cells that fall on the conductor            
                    if G1[x,y]== -2:
                        b[(G1[r,c])]=b[(G1[r,c])]+ev
                    elif G1[x,y]== -3:
                        b[(G1[r,c])]=b[(G1[r,c])]-ev
        #Do nothing with ground cells
                    elif G1[x,y]== -1:
                        pass
        #Process nodes that are in the final column only
                    else:
                        if c==Gcol:
                            Ar[Ainc]=e
                            Ac[Ainc]=f
                            Ad[Ainc]=Ad[Ainc]-1
                        #Process nodes
                        Ar[Ainc]=e
                        Ac[Ainc]=f
                        Ad[Ainc]=Ad[Ainc]-1
                        Ainc=Ainc+1
    
    #Load node equation data into sparse matrix
    n=np.size(np.nonzero(Ac)[0])
    A=sp.csr_matrix((Ad,(Ar,Ac)),shape=(max(Ar)+1,max(Ac)+1))
#    A=A.todense()
    #Initialise matricies for diagonal line
    #Populate diagonal
    RC=np.r_[0:node]
    #for r in range(0,node):
    #    RC[r]=r;
    #    DiagD[r]=4;
    #Add diagonal values
    #D=coo_matrix(DiagD,(RC,RC))
    D=sp.csr_matrix(np.identity(node+1)*4)
    A=A+D
    #Perform node voltage equation
    v=spsolve(A,b)
#    v=np.linalg.solve(A,b)
#    v=A/b
    Q=np.zeros((Grow,Gcol))
    #Process node voltage matrix cell by cell
    for r in range(0,Grow):
        for c in range(0,Gcol):
    #Set conductor voltage
            if G1[r,c]== -2:
                Q[r,c]=ev
            elif G1[r,c]== -3:
                Q[r,c]=-ev
    #Set ground voltage
            elif G1[r,c]== -1:
                Q[r,c]=0
     #Set node voltage
            else:
                Q[r,c]=v[(G1[r,c])]
    return Q, node, A, v
    
def trace_values(Q,ev,E,u,vp,Er,c,traces):
    Von=0
    Grow=Q.shape[0]
    Gcol=Q.shape[1]
    

    if traces == 1:
       pass 
    elif traces ==2:
        Von = np.sum(Q[1,1:(Gcol/2)]) + np.sum(Q[Grow-2,1:(Gcol/2)]) + np.sum(Q[2:Grow-2,1])
        + np.sum(Q[2:Grow-2,(Gcol/2)-1]) 
        Vin = np.sum(Q[2,2:(Gcol/2)-1]) + np.sum(Q[Grow-3,2:(Gcol/2)-1]) + np.sum(Q[3:Grow-3,2]) 
        + np.sum(Q[3:Grow-3,(Gcol/2)-2]) 
    #Calculate charge
    q=E*(Vin-Von)
    #Calculate capacitance
    C=q/(2*ev)    
    #Calculate inductance
    L=1/(vp*np.sqrt(C))**2
    #Calculate characteristic impedance
    Z0=(np.sqrt(u*E))/C
    
#    Z0 = (np.sqrt(Er))/(c*C)    
    return q,C,Z0,L
    
def geometry_matrix(geometry,Grow,Gcol,trace_width_cells,trace_spacing_cells,
                    trace_thickness_cells):
    #Initialise geometry matrix
    G1=np.zeros((Grow,Gcol))
    if geometry == 'stripline':
        #Populate geometry matrix for conductor and ground
        G1[0,:] = -1
        G1[Grow-1,:] = -1
        G1[np.round((Grow/2))+1,np.round(((Gcol/2))-
        (trace_width_cells/2)):np.round(((Gcol/2))+(trace_width_cells/2))] = -2
        traces=1
    elif geometry == 'edge coupled stripline':
        #Populate geometry matrix for conductor and ground
        if np.round(trace_thickness_cells/2)>1:
            gr1=((Grow-1)/2)-np.round(trace_thickness_cells/2)
            gr2=((Grow-1)/2)+np.round(trace_thickness_cells/2)
        else:
            gr1=((Grow-1)/2)
            gr2=gr1+1
        G1[0,:] = -1
        G1[Grow-1,:] = -1
        G1[gr1:gr2,np.round(((Gcol/2))-trace_width_cells-(trace_spacing_cells/2)):np.round(((Gcol/2))-(trace_spacing_cells/2))] = -2
        G1[gr1:gr2,np.round(((Gcol/2))+(trace_spacing_cells/2)):
            np.round(((Gcol/2))+trace_width_cells+(trace_spacing_cells/2))] = -3
        traces=2
    elif geometry == 'edge coupled microstripline':
        #Populate geometry matrix for conductor and ground
        G1[0,:] = -1
        G1[Grow-1,np.round(((Gcol/2))-trace_width_cells-
        (trace_spacing_cells/2)):np.round(((Gcol/2))-(trace_spacing_cells/2))] = -2
        G1[Grow-1,np.round(((Gcol/2))+(trace_spacing_cells/2)):
            np.round(((Gcol/2))+trace_width_cells+(trace_spacing_cells/2))] = -3
        traces=2
    return G1,traces