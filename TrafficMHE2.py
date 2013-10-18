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


import TrafficFlow2
reload(TrafficFlow2)
from TrafficFlow2 import *

Horizon = 10
T = FreeWay(Path = 'incident data_correct VMS', Meas = ['rho','v'], Slacks = ['Salpha','Sbeta'])

SimTime = 1600

TrafficParameters = { 'a'     :   4.04, 
                      'tau'   :   16.1,
                      'kappa' :   6.84,
                      'rho_cr':  28.56,
                      'eta'   :  32.13,
                      'T'     :  10/3600.,
                      'delta' :   0.,
                      'vfree' : 107.11}


vfree  = TrafficParameters[ 'vfree'  ]
a      = TrafficParameters[ 'a'      ]
rho_cr = TrafficParameters[ 'rho_cr' ]

# vv_3lane   -> speed (v)
# rr_3lane/3 -> density (rho)
# VMS_d      -> variable speed limit 0> Ve = min

Q = {'rho': 1., 'v': 1.}

Costs = {'Stage': 0, 'Terminal': 0, 'Arrival': 0}
Const = []
for key in Costs.keys():
    for i in range(T.NumSegment): 
        Costs[key] += (T.VSpace['Slacks']['Segment',i,'Salpha'])
        Costs[key] += (T.VSpace['Slacks']['Segment',i, 'Sbeta'])
        Costs[key] += Q['rho']*(T.VSpace['States']['Segment',i,'rho'] - T.Meas['Segment',i,'rho'])**2  
        Costs[key] += Q[  'v']*(T.VSpace['States']['Segment',i,  'v'] - T.Meas['Segment',i,  'v'])**2 
    
        if key == 'Stage':
            dalpha = T.VSpace['Inputs']['Segment',i,'alpha'] - T.VSpacePrev['Inputs']['Segment',i,'alpha']
            dbeta  = T.VSpace['Inputs']['Segment',i, 'beta'] - T.VSpacePrev['Inputs']['Segment',i, 'beta']
             
            Const.append(-dbeta  - T.VSpace['Slacks']['Segment',i, 'Sbeta']) # <= 0 
            Const.append( dbeta  - T.VSpace['Slacks']['Segment',i, 'Sbeta']) # <= 0
            Const.append(-dalpha - T.VSpace['Slacks']['Segment',i,'Salpha']) # <= 0 
            Const.append( dalpha - T.VSpace['Slacks']['Segment',i,'Salpha']) # <= 0
            
            Ve_arg   = (1 + T.VSpace['Inputs']['Segment',i,'alpha'])*T.VSpace['States']['Segment',i,'rho']
            Ve_i_ref = vfree*exp(-((Ve_arg/rho_cr)**a)/a)
            Costs[key] += 1e1*(T.VSpace['Inputs']['Segment',i,'Ve'] - Ve_i_ref)**2
            
        #if key == 'Arrival':
        #    Costs[key] += T.VSpace['States']['Segment',i,'rho'] - T.VSpacePrev['States']['Segment',i,'rho']
        #    Costs[key] += T.VSpace['States']['Segment',i,  'v'] - T.VSpacePrev['States']['Segment',i,  'v']

        
            
            
    T.setCost(Costs[key], Type = key)

T.setIneqConst(Const)

T.BuildMHE(Horizon = Horizon, SimTime = SimTime)




#Assign initial Data and Parameters
Data, EPData = T.GenerateData(VStruct = T.VSim(), EPStruct = T.EPSim(), ExtParam = TrafficParameters, Simulated = False)#, AddNoise = True)


print "HIT RETURN TO LAUNCH MHE"

raw_input()
#################################

init, EP, lbV, ubV = T.PassMHEData(Data = Data, EP = EPData, ExtParam = TrafficParameters, time = 0)

MHETraj = T.VSim()
StatusLog = []
for time in range(SimTime):
    print time

    #Cost, g = T.Check(init = init, EP = EP)
    #print "Check initial guess:"
    #print "Cost0", Cost, " | Const", list(g['EqConst',veccat])
    
    #Error = []
    #for i in range(T.NumSegment):
    #    Error.append(veccat(EP['Meas',:,'Segment',i,'rho']) - veccat(init['States',:,'Segment',i,'rho']))
    #
    #assert(time == 0)
    X, Mu, Status = T.SolveMHE(EP = EP, lbV = lbV, ubV = ubV, init = init)
    StatusLog.append(Status)
    
    #Display
    timeDisp = {'Sim'      : [k      for k in range(SimTime)    ],
                'MHEState' : [k+time for k in range(Horizon)    ],
                'MHEInput' : [k+time for k in range(Horizon-1)  ],
                'Data'     : [k      for k in range(Horizon+SimTime) ]}
    
    #if (time > -1):
    #    plt.close('all')
    #    for i in range(T.NumSegment):
    #        plt.figure(i+1)
    #        plt.hold('on')
    #        plt.subplot(2,1,1)
    #        plt.plot(timeDisp['MHEState'], X['States',:,'Segment',i,'v'],color = 'r',linewidth = 2)
    #        plt.plot(timeDisp['Data'],     Data['States',:,'Segment',i,'v'],color = 'g')
    #        plt.title('v')
    #        
    #        plt.subplot(2,1,2)
    #        plt.plot(timeDisp['MHEState'], X['States',:,'Segment',i,'rho'],color = 'r',linewidth = 2)
    #        plt.plot(timeDisp['Data'],     Data['States',:,'Segment',i,'rho'],color = 'g')
    #        plt.title('rho')
    #        
    #    plt.figure(42)
    #    plt.hold('on')
    #    for i in range(1,T.NumSegment-1):
    #        plt.subplot(T.NumSegment-2,1,i)
    #        plt.plot(timeDisp['MHEInput'],X['Inputs',:,'Segment',i,'alpha'],color = 'k', label = 'alpha', linestyle = '--')
    #        plt.plot(timeDisp['MHEInput'],X['Inputs',:,'Segment',i,'beta'], color = 'k', label = 'beta')
    #        plt.ylim([-0.1,1.1])
    #        plt.grid()
    #    plt.legend()
    #    raw_input()
    #    assert(0==1)


    MHETraj[...,time] = X[...,0]
    #QP = T.PrepareQP(T.Solver, Primal = X, Adjoint = Mu, EP = EP, lbV = lbV, ubV = ubV)
    
    _, EP, lbV, ubV = T.PassMHEData(Data = Data, EP = EPData, ExtParam = TrafficParameters, time = time+1)
    init  = T.Shift(X, EP)

    

timeDisp = {'Sim'      : [k      for k in range(SimTime)    ],
            'MHEState' : [k+time for k in range(Horizon)    ],
            'MHEInput' : [k+time for k in range(Horizon-1)  ],
            'Data'     : [k      for k in range(T.DataTime) ]}
    
plt.close('all')
for i in range(T.NumSegment):
    plt.figure(i+1)
    plt.hold('on')
    plt.subplot(3,1,1)  
    plt.plot(timeDisp['Sim'],      MHETraj['States',:SimTime,'Segment',i,'v'],color = 'r',linewidth = 2)
    plt.plot(timeDisp['Sim'],      Data['States',:SimTime,'Segment',i,'v'],color = 'b')
    plt.grid()
    plt.title('v')
    
    plt.subplot(3,1,2)
    plt.plot(timeDisp['Sim'],      MHETraj['States',:SimTime,'Segment',i,'rho'],color = 'r',linewidth = 2)
    plt.plot(timeDisp['Sim'],      Data['States',:SimTime,'Segment',i,'rho'],color = 'b')
    plt.grid()
    plt.title('rho')
    
    plt.subplot(3,1,3)
    plt.plot(timeDisp['Sim'],     T.Data['VMS_d'][:SimTime,3*i], linestyle = 'none', marker = '.',color = 'k')
    plt.plot(timeDisp['Sim'],     MHETraj['Inputs',:SimTime,'Segment',i,'Ve'],color = 'r',linewidth = 2)
    plt.plot(timeDisp['Sim'],     Data['Inputs',:SimTime,'Segment',i,'Ve'],color = 'g')
    plt.ylim([-10,120])
    plt.grid()
    plt.title('Ve')
    
plt.figure(42)
plt.hold('on')
for i in range(1,T.NumSegment-1):
    plt.subplot(T.NumSegment-2,1,i)
    plt.plot(timeDisp['Sim'], MHETraj['Inputs',:SimTime,'Segment',i,'alpha'],color = 'k', label = 'alpha', linestyle = '--',linewidth = 2)
    plt.plot(timeDisp['Sim'], MHETraj['Inputs',:SimTime,'Segment',i,'beta'], color = 'k', label = 'beta',linewidth = 2)
    plt.ylim([-0.1,1.1])
    plt.grid()
plt.legend()

plt.figure(99)
plt.plot(StatusLog, linestyle = 'none', marker = '.',color = 'k')
plt.ylim([-0.2,1.2])
raw_input()
    
#
