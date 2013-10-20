# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 20:18:08 2012

@author: Sebastien Gros

Assistant Professor
 
Department of Signals and Systems
Chalmers University of Technology
SE-412 96 Gšteborg, SWEDEN
grosse@chalmers.se

Python/casADi Module:
A Fast Algorithm for Power Smoothing of Wind Farms based on Distributed Optimal Control

Requires the installation of the open-source Python module casADi together with the NLP solver ipopt

Required version of CasADi: v1.7.x
 
"""

from casadi import *
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from matplotlib import interactive
interactive(True)

import time as timer

import TrafficFlow
reload(TrafficFlow)
from TrafficFlow import *

Horizon = 5
T = FreeWay(Path = 'incident data_correct VMS', Meas = ['rho','v'], Slacks = ['Salpha', 'Sbeta', 'SVe'], ExtParam = ['Ve_max'])

SimTime = 800

TrafficParameters = { 'a'     :   4.04, 
                      'tau'   :   16.1,
                      'kappa' :   6.84,
                      'rho_cr':  28.56,
                      'eta'   :  32.13,
                      'T'     :  10/3600.,
                      'delta' :     0.,
                      'vfree' : 107.11}


vfree  = TrafficParameters[ 'vfree'  ]
a      = TrafficParameters[ 'a'      ]
rho_cr = TrafficParameters[ 'rho_cr' ]

# vv_3lane   -> speed (v)
# rr_3lane/3 -> density (rho)
# VMS_d      -> variable speed limit 0> Ve = min

Q = {'rho': 10., 'v': 10., 'Ve': 10.}

Costs = {'Stage': 0, 'Terminal': 0, 'Arrival': 0}
Const = []
for key in Costs.keys():
    for i in range(T.NumSegment): 
        Costs[key] += (T.VSpace['Slacks']['Salpha',i])
        Costs[key] += (T.VSpace['Slacks']['Sbeta', i])
        Costs[key] += 1e-3*(T.VSpace['Slacks']['SVe',   i])**2
        
        Costs[key] += Q['rho']*(T.VSpace['States']['rho',i] - T.Meas['rho',i])**2  
        Costs[key] += Q[  'v']*(T.VSpace['States']['v',  i] - T.Meas['v',  i])**2 
    
        if key == 'Stage':
            dalpha = T.VSpace['Inputs']['alpha',i] - T.VSpacePrev['Inputs']['alpha',i]
            dbeta  = T.VSpace['Inputs']['beta', i] - T.VSpacePrev['Inputs']['beta', i]
             
            Const.append(-dbeta  - T.VSpace['Slacks']['Sbeta', i]) # <= 0 
            Const.append( dbeta  - T.VSpace['Slacks']['Sbeta', i]) # <= 0
            Const.append(-dalpha - T.VSpace['Slacks']['Salpha',i]) # <= 0 
            Const.append( dalpha - T.VSpace['Slacks']['Salpha',i]) # <= 0
            
            Ve_arg   = (1 + T.VSpace['Inputs']['alpha',i])*T.VSpace['States']['rho',i]
            Ve_i_ref = vfree*exp(-((Ve_arg/rho_cr)**a)/a)
            Costs[key] += Q['Ve']*(T.VSpace['Inputs']['Ve',i] - Ve_i_ref)**2
            
                
            
            #Const.append( T.VSpace['Inputs']['Ve',i]  - T.ExtParam['Ve_max',i] - T.VSpace['Slacks']['SVe', i])
            
        if key == 'Arrival':
            Costs[key] += T.VSpace['States']['rho',i] - T.VSpacePrev['States']['rho',i]
            Costs[key] += T.VSpace['States']['v',  i] - T.VSpacePrev['States']['v',  i]

        
            
            
    T.setCost(Costs[key], Type = key)

T.setIneqConst(Const)

T.BuildMHE(Horizon = Horizon, SimTime = SimTime, Tol = 1e-3)


#Assign initial Data and Parameters
Data, EPData = T.GenerateData(VStruct = T.VSim(), EPStruct = T.EPSim(), ExtParam = TrafficParameters, Simulated = True)#, AddNoise = True)


print "HIT RETURN TO LAUNCH MHE"

raw_input()
start = timer.time()
#################################

init, EP, lbV, ubV = T.PassMHEData(Data = Data, EP = EPData, ExtParam = TrafficParameters, time = 0)

MHETraj = T.VSim()
StatusLog = []
for time in range(SimTime):
    print time

    ubV['Inputs',:,'Ve',:] = 200
    EP['ExtParam',:,'Ve_max',:] = 200
    for i in range(T.NumSegment):
        for k in range(T.Horizon-1):
            VeBound = np.mean(T.Data['VMS_d'][time+k,3*i:3*i+3])
            if np.isnan(VeBound):
                VeBound = 200.
            EP['ExtParam',k,'Ve_max',i] = VeBound
    
    
    #Cost,g = T.Check(init = init, EP = EP)
    #print float(Cost)
    #print g['EqConst',veccat]
    #assert(0==1)
    X, Mu, Status = T.SolveMHE(EP = EP, lbV = lbV, ubV = ubV, init = init)
    StatusLog.append(Status)
    
    ##Display
    #timeDisp = {'Sim'      : [k      for k in range(SimTime)    ],
    #            'MHEState' : [k+time for k in range(Horizon)    ],
    #            'MHEInput' : [k+time for k in range(Horizon-1)  ],
    #            'Data'     : [k      for k in range(Horizon+SimTime) ]}
    #
    #if (time > -1):
    #    plt.close('all')
    #    for i in range(T.NumSegment):
    #        plt.figure(i+1)
    #        plt.hold('on')
    #        plt.subplot(2,2,1)
    #        plt.plot(timeDisp['MHEState'], X['States',:,'v',i],color = 'r',linewidth = 2)
    #        plt.plot(timeDisp['Data'],     Data['States',:,'v',i],color = 'g')
    #        plt.title('v')
    #        
    #        plt.subplot(2,2,2)
    #        plt.plot(timeDisp['MHEState'], X['States',:,'rho',i],color = 'r',linewidth = 2)
    #        plt.plot(timeDisp['Data'],     Data['States',:,'rho',i],color = 'g')
    #        plt.title('rho')
    #        
    #        plt.subplot(2,2,3)
    #        plt.plot(timeDisp['MHEInput'], X['Inputs',:,'Ve',i],color = 'r',linewidth = 2)
    #        plt.plot(timeDisp['Data'],     Data['Inputs',:,'Ve',i],color = 'g')
    #        plt.ylim([-10,120])
    #        plt.grid()
    #        plt.title('Ve')
    #
    #        
    #        plt.subplot(2,2,4)
    #        plt.plot(timeDisp['MHEInput'], X['Inputs',:,'alpha',i],color = 'k', label = 'alpha', linestyle = '--',linewidth = 2)
    #        plt.plot(timeDisp['MHEInput'], X['Inputs',:,'beta',i], color = 'k', label = 'beta',linewidth = 2)
    #        plt.ylim([-0.1,1.1])
    #        plt.grid()
    #        plt.legend(loc = 0)
    #    raw_input()
    #    #assert(0==1)


    MHETraj[...,time] = X[...,0]
    #QP = T.PrepareQP(T.Solver, Primal = X, Adjoint = Mu, EP = EP, lbV = lbV, ubV = ubV)
    
    _, EP, lbV, ubV = T.PassMHEData(Data = Data, EP = EPData, ExtParam = TrafficParameters, time = time+1)
    init  = T.Shift(X, EP)

elapsed = (timer.time() - start)

print "Total time: ", elapsed

timeDisp = {'Sim'      : [k      for k in range(SimTime)    ],
            'MHEState' : [k+time for k in range(Horizon)    ],
            'MHEInput' : [k+time for k in range(Horizon-1)  ],
            'Data'     : [k      for k in range(T.DataTime) ]}
    
plt.close('all')
for i in range(T.NumSegment):
    plt.figure(i+1)
    plt.hold('on')
    plt.subplot(2,2,1)  
    plt.plot(timeDisp['Sim'],      MHETraj['States',:SimTime,'v',i],color = 'r',linewidth = 2)
    plt.plot(timeDisp['Sim'],      Data['States',:SimTime,'v',i],color = 'b')
    plt.grid()
    plt.title('v')
    
    plt.subplot(2,2,2)
    plt.plot(timeDisp['Sim'],      MHETraj['States',:SimTime,'rho',i],color = 'r',linewidth = 2)
    plt.plot(timeDisp['Sim'],      Data['States',:SimTime,'rho',i],color = 'b')
    plt.grid()
    plt.title('rho')
    
    plt.subplot(2,2,3)
    plt.plot(timeDisp['Sim'],     np.mean(T.Data['VMS_d'][:SimTime,3*i:3*(i+1)],axis=1), linestyle = 'none', marker = '.',color = 'k')
    plt.plot(timeDisp['Sim'],     MHETraj['Inputs',:SimTime,'Ve',i],color = 'r',linewidth = 2)
    plt.plot(timeDisp['Sim'],     Data['Inputs',:SimTime,'Ve',i],color = 'g')
    plt.ylim([-10,120])
    plt.grid()
    plt.title('Ve')
    
    plt.subplot(2,2,4)
    plt.plot(timeDisp['Sim'], MHETraj['Inputs',:SimTime,'alpha',i],color = 'k', label = 'alpha', linestyle = '--',linewidth = 2)
    plt.plot(timeDisp['Sim'], MHETraj['Inputs',:SimTime,'beta',i], color = 'k', label = 'beta',linewidth = 2)
    plt.ylim([-0.1,1.1])
    plt.grid()
    plt.legend(loc = 0)
    
#plt.figure(42)
#plt.hold('on')
#for i in range(1,T.NumSegment-1):
#    plt.subplot(T.NumSegment-2,1,i)
#    plt.plot(timeDisp['Sim'], MHETraj['Inputs',:SimTime,'alpha',i],color = 'k', label = 'alpha', linestyle = '--',linewidth = 2)
#    plt.plot(timeDisp['Sim'], MHETraj['Inputs',:SimTime,'beta',i], color = 'k', label = 'beta',linewidth = 2)
#    plt.ylim([-0.1,1.1])
#    plt.grid()
#plt.legend()

plt.figure(99)
plt.plot(StatusLog, linestyle = 'none', marker = '.',color = 'k')
plt.ylim([-0.2,1.2])
raw_input()
    
#
