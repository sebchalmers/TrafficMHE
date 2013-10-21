# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16, 2014

@author:

Sebastien Gros
Assistant Professor
 
Department of Signals and Systems
Chalmers University of Technology
SE-412 96 Gšteborg, SWEDEN
grosse@chalmers.se

Python/casADi Module:
An MHE module for Freeway Traffic Incident Detection

Requires the installation of the open-source Python module casADi together with the NLP solver ipopt

Required version of CasADi: v1.7.x
 
"""

from casadi import *
from casadi.tools import *
#import math as math

import numpy as np
import scipy.special as sc_spec
from scipy import linalg
import scipy.io
import random as rand

import matplotlib.pyplot as plt
from pylab import matshow
from matplotlib import interactive

def assertList(var):
    if not(isinstance(var,list)):
        var = [var]
    return var

#Solver Constructor (used in both Turbine and WindFarm class)
def _setSolver(self, V,Cost,g,P, Tol):

    lbg = g()
    ubg = g()
    if ('IneqConst' in g.keys()):
        lbg['IneqConst'] = -inf
        ubg['IneqConst'] =  0
    EqConst          = list(veccat(g.i['EqConst'  ]))
    
    if 'IneqConst' in g.keys():   
        IneqConst = list(veccat(g.i['IneqConst']))
    else:
        IneqConst = []

    nl = MXFunction(nlpIn(x=V,p=P),nlpOut(f=Cost,g=g))
    nl.init()
    
    # set-up solver
    solver = IpoptSolver(nl)  
    solver.setOption("expand",True)
    solver.setOption("print_level",0)  
    solver.setOption("hessian_approximation","exact")
    solver.setOption("max_iter",1500)
    solver.setOption("tol",Tol)
    solver.setOption("linear_solver","ma27")
    
    solver.init()
    
    print "Solver constructed, Extract Hessian & gradient"
    
    H = solver.hessLag()
    H.init()
    dg = solver.jacG()
    dg.init()
    
    gfunc = MXFunction([V,P],[g])
    gfunc.init()
    
    f = solver.gradF() 
    f.init()
    
    Solver = {'solver': solver, 'H': H, 'f': f, 'dg': dg, 'g': gfunc, 'lbg': lbg, 'ubg': ubg, 'V': V, 'EqConst':EqConst,  'IneqConst':IneqConst}
    return Solver

def _CreateQP(solver,V):

    solver['H'].evaluate()
    H = solver['H'].output()# + np.eye(solver['H'].output().shape[0])
    
    solver['dg'].evaluate()
    dg = solver['dg'].output()
    
    QPsolver = NLPQPSolver(qpStruct(h=H.sparsity(),a=dg.sparsity()))
    QPsolver.setOption({"nlp_solver":IpoptSolver, "nlp_solver_options": {"tol": 1e-8,"verbose": 0} })
   
    QPsolver.init()

    return QPsolver  

def _CreateFunc(Input,Output):
    Func = MXFunction(Input,Output)        
    Func.init()
    return Func


class FreeWay:
    def _StructBuilder(self, Vars):
            EntryList = []
            for key in Vars:
                Struct = struct_ssym([key])
                EntryList.append(entry(key, struct = Struct, repeat = self.NumSegment))      
            Struct = struct_msym(EntryList)
            return Struct
        
    def __init__(self, Path = [], Slacks = [], Meas = [], ExtParam = []):
        print "------------------"
        print "Initialize FreeWay"
        print "------------------"
        print "Load Data"
        try:
            
            self.Data = scipy.io.loadmat(Path)
        except:
            print "Unkown Data file, give path/name", Path
                           
        
        print "Build Symbolic variables"
        
        self.NumSegment = self.Data['L'].shape[1]+1
        self.DataTime   = self.Data['vv_3lane'].shape[0]-1
        
        VarList = {
                    'States': ['rho', 'v'],
                    'Inputs': ['r','f','alpha','beta', 'Ve']
                  }
        
        GlobalParams  = struct_ssym(['vfree', 'kappa', 'tau', 'a', 'rho_cr', 'eta', 'T', 'delta'])
        Li            = struct_ssym(['L'])
        
        
        #Build structure to hold the data (will merge into EP)
        if not(Meas == []):
            self.Meas = self._StructBuilder(Meas)
        else:
            print "WARNING: NO DATA PROVIDED"
            

        if not(Slacks == []):
            # Slacks declared by the user (for constraints relaxation)
            print "Slack variables detected:"
            print Slacks
            VarList['Slacks'] = Slacks
                                 
        #Spatial Structure
        self.VSpace     = {}
        self.VSpacePrev = {}
        for var in VarList.keys():  
            self.VSpace[var]     = self._StructBuilder(VarList[var])
            self.VSpacePrev[var] = self._StructBuilder(VarList[var])
        
        #Create Structure for External parameters
        ParamList =                 [
                                        entry('Global',  struct = GlobalParams                                      ),
                                        entry('States0', struct = self.VSpace['States']                             ),
                                        entry('Inputs0', struct = self.VSpace['Inputs']                             ),
                                        entry('Dist',    struct = Li,                      repeat = self.NumSegment )
                                    ]
        self.Param  = struct_msym(ParamList)
        
        if not(ExtParam == []):
            ExtParamList = []
            for key in ExtParam:
                Struct = struct_ssym([key])
                ExtParamList.append(entry(key,    struct = Struct, repeat = self.NumSegment ))
            self.ExtParam = struct_msym(ExtParamList)
        
        
        
        T      = self.Param['Global',     'T']
        a      = self.Param['Global',     'a']
        delta  = self.Param['Global', 'delta']
        kappa  = self.Param['Global', 'kappa']
        tau    = self.Param['Global',   'tau']
        vfree  = self.Param['Global', 'vfree']
        rho_cr = self.Param['Global','rho_cr']
        eta    = self.Param['Global',   'eta']

        
        print "Build Stage-wise Dynamic Function"
        self._Shoot = ['Boundary, dynamics n/a']
        for i in range(1,self.NumSegment-1):
            ##Attribute spatial dynamic elements
            rho_i   = self.VSpace['States']['rho',  i] 
            v_i     = self.VSpace['States']['v',    i] 
            
            rho_im  = self.VSpace['States']['rho',i-1]  
            v_im    = self.VSpace['States']['v',  i-1]              
                    
            rho_ip  = self.VSpace['States']['rho',i+1]  
            v_ip    = self.VSpace['States']['v',  i+1] 
                            
            q_im    = rho_im*v_im
            q_i     = rho_i *v_i
            
            r_i     = self.VSpace['Inputs']['r',    i]  
            f_i     = self.VSpace['Inputs']['f',    i]  
            
            beta_i  = self.VSpace['Inputs']['beta', i]  
            alpha_i = self.VSpace['Inputs']['alpha',i]  
            
            Li = self.Param['Dist',i,'L']
            
            #Ve function
            #Ve_arg  = (1+alpha_i)*rho_i
            #Ve_i    = vfree*exp(-((Ve_arg/rho_cr)**a)/a)
            Ve_i    = self.VSpace['Inputs']['Ve',  i]  
            
            #Construct Dynamics        
            Dyn = []
                    
            Dyn.append(   rho_i + T*(q_im - q_i + r_i - f_i)/Li )
            
            Dyn.append(     v_i + T*(beta_i*Ve_i - v_i)/tau    \
                                + T*v_i*(v_im - v_i)/Li  \
                                - beta_i*(1 - alpha_i)*(eta*T)*(rho_ip - rho_i)/(rho_i + kappa)/(tau*Li)         \
                                - delta*T*r_i*v_i/Li/(rho_i + kappa) )
            

            Shoot = _CreateFunc([self.VSpace['Inputs'],self.VSpace['States'],self.Param],[veccat(Dyn)])
            self._Shoot.append(Shoot)
            
        self._Shoot.append('Boundary, dynamics n/a')

    def _BuildFunc(self, Expr, Type):

        X           = self.VSpace['States']
        U           = self.VSpace['Inputs']
            
        Xprev       = self.VSpacePrev['States']
        Uprev       = self.VSpacePrev['Inputs']
                        
        #CAREFUL, THE INPUT SCHEME MUST MATCH THE FUNCTION CALLS !!
        if Type == 'Terminal':
            listFuncInput = [X, Xprev]             #, Slacks; if applicable
        else:
            listFuncInput = [X, Xprev, U, Uprev]   #, Slacks; if applicable         
        
        if ('Slacks' in self.VSpace.keys()):
            listFuncInput.append(self.VSpace['Slacks'])#['Slacks'])
            
        listFuncInput.append(self.Param)
        listFuncInput.append(self.Meas)
        
        if hasattr(self,'ExtParam'):
            listFuncInput.append(self.ExtParam)
        
        self.listFuncInput = listFuncInput
        Func = MXFunction(listFuncInput,[Expr])
        Func.init()
        
        return Func
    
    def setIneqConst(self, Const, Type = 'Stage'):
        """
        Sets the Cost of your MHE problem. Send in a symbolic expression constructed from .VSpace symbols as argument.
        Specify Terminal = True if this is your Terminal constraint.
        """
        if not(isinstance(Const,list)):
            Const = [Const]
        
        #ConstFunc = self._BuildFunc(veccat(Const), Terminal)
        if (Type == 'Stage'):
            self._StageConst    = self._BuildFunc(veccat(Const), Type)
        if (Type == 'Arrival'):
            self._ArrivalConst  = self._BuildFunc(veccat(Const), Type)
        if (Type == 'Terminal'):
            self._TerminalConst = self._BuildFunc(veccat(Const), Type)
                    
    def setCost(self, Expr, Type = 'Stage'):
        """
        Sets the Cost of your MHE problem. Send in a symbolic expression constructed from .VSpace symbols as argument.
        Specify Arrival = True if this is your arrival cost.
        """
        if (Type == 'Stage'):
            self._StageCost    = self._BuildFunc(Expr, Type)
        if (Type == 'Arrival'):
            self._ArrivalCost  = self._BuildFunc(Expr, Type)
        if (Type == 'Terminal'):
            self._TerminalCost = self._BuildFunc(Expr, Type)            
    
    def _AddCostAndConst(self, Cost, IneqConst, k):
        
        StageInputList = []
        AddCost        = 0
                
        if (k==0):
            StageInputList.append(self.V['States', k    ])
            StageInputList.append(self.EP['Param','States0'     ])
            StageInputList.append(self.V[ 'Inputs', k   ])
            StageInputList.append(self.EP['Param','Inputs0'     ])
            
        if (k > 0) and (k < self.Horizon-1):
            StageInputList.append(self.V['States', k    ])
            StageInputList.append(self.V['States', k-1  ])
            StageInputList.append(self.V['Inputs', k    ])
            StageInputList.append(self.V['Inputs', k-1  ])
            
        if (k == self.Horizon-1):
            StageInputList.append(self.V['States', k    ])
            StageInputList.append(self.V['States', k-1  ])
            
        
        if ('Slacks' in self.VSpace.keys()):
            StageInputList.append(self.V['Slacks',k])

        StageInputList.append(self.EP['Param'])
        StageInputList.append(self.EP['Meas', k])
        
        if ('ExtParam' in self.EP.keys()):
            StageInputList.append(self.EP['ExtParam',k])
    
        if k == 0:
            if hasattr(self, '_ArrivalCost'):
                print "Arrival Cost Detected, build."
                [AddCost ]    = self._ArrivalCost.call( StageInputList)
                Cost += AddCost
                
            if hasattr(self, '_ArrivalConst'):
                print "Arrival Constraint Detected, build."
                [AppendConst] = self._ArrivalConst.call(StageInputList)
                IneqConst.append( AppendConst )
        
        if k == self.Horizon-1:
            if hasattr(self, '_TerminalCost'):
                print "Terminal Cost Detected, build."
                [AddCost ]    = self._TerminalCost.call( StageInputList)
                Cost += AddCost
                
            if hasattr(self, '_TerminalConst'):
                print "Terminal Constraint Detected, build."
                [AppendConst] = self._TerminalConst.call(StageInputList)
                IneqConst.append( AppendConst )
                
        if (k < self.Horizon-1):
            if hasattr(self, '_StageCost'):
                print "Stage Cost Detected, build."
                [AddCost ] = self._StageCost.call(           StageInputList)
                Cost += AddCost
                
            if hasattr(self, '_StageConst'):
                print "Stage Constraint Detected, build."
                [AppendConst] = self._StageConst.call(   StageInputList)
                IneqConst.append( AppendConst )
        

        
        
        return Cost, IneqConst
    
    def _CreateEPStruct(self, TimeHorizon):
        EPList = [
                    entry('Param', struct = self.Param),
                    entry('Meas',  struct = self.Meas,      repeat = TimeHorizon)
                  ]
        
        if hasattr(self,'ExtParam'):
            EPList.append(entry('ExtParam',  struct = self.ExtParam,      repeat = TimeHorizon))
            
        EP = struct_msym(EPList)
        return EP
    
    def BuildMHE(self, Horizon = 1, SimTime = 1, Tol = 1e-3):
        """
        Construct the MHE problem. Specify your horizon lengths as Horizon = ...
        """
        if (Horizon > 0):
            self.Horizon = Horizon
        else:
            print "Warning, no horizon specified. Default value is 1 !!"
        
        self.SimTime = SimTime + Horizon
        
        NumElement = {'States': Horizon, 'Inputs' : Horizon-1}
        if ('Slacks' in self.VSpace.keys()):
            NumElement['Slacks'] = Horizon
        
        print "Build Freeway Variable Structures"
        VarStruct = []
        for var in NumElement.keys():
            VarStruct.append(  entry(var, struct = self.VSpace[var], repeat = NumElement[var]  ))
        
        
        SimStruct = []
        for var in NumElement.keys():
            SimStruct.append(  entry(var, struct = self.VSpace[var], repeat = self.SimTime  ))
        
        self.V    = struct_msym(VarStruct)
        self.VSim  = struct_msym(SimStruct)
        
        self.EP    = self._CreateEPStruct(self.Horizon)        
        self.EPSim = self._CreateEPStruct(self.SimTime)

        #Construct Dynamic Constraints
        print "Building Dynamic Constraints"
        EquConst  = []
        for k in range(0,Horizon-1):
            #print "Time", k
            for i in range(1,self.NumSegment-1):
                #print "Segment", i
                [Xplus]     = self._Shoot[i].call([self.V['Inputs',k],self.V['States',k],self.EP['Param']])
                self.Xplus = Xplus
                EquConst.append(  veccat(self.V['States',k+1,...,i]) - Xplus )
                
        #Construct Stage Cost & Inequality Constraints
        print "Building Cost & Inequality Constraints"
        IneqConst = []
        Cost      = 0 
        for k in range(Horizon):
            Cost, IneqConst = self._AddCostAndConst(Cost, IneqConst, k)

        print "\n"
        print "------------------"
        print "Construct Solvers"
        print "------------------"
        print "\n"
        print "NLP"
        
        gList = [entry('EqConst',   expr = EquConst)]
        if not(IneqConst == []):
            gList.append(entry('IneqConst', expr = IneqConst))
        self.g = struct_MX(gList)
        
        ##Setup Central Solver
        self.Solver = _setSolver(self, self.V, Cost, self.g, self.EP, Tol = Tol)
        
        self.gfunc = MXFunction([self.V, self.EP],[self.g])
        self.gfunc.init()
        self.Costfunc = MXFunction([self.V, self.EP],[Cost])
        self.Costfunc.init()
        
        ##Setup Central Solver
        print "QP"
        self.QPSolver = _CreateQP(self.Solver, self.V)
        print "\n"
        print "Construction Terminated."
        
        
    def PrepareQP(self, solver = [], Primal = [], Adjoint = [], EP = [], lbV = [], ubV = []):
            
        solver['H'].setInput(Primal, 0)
        solver['H'].setInput(EP     ,1)
        solver['H'].setInput(1.     ,2)
        solver['H'].setInput(Adjoint,3)
        solver['H'].evaluate()
    
        solver['f'].setInput(Primal, 0)
        solver['f'].setInput(EP,     1)
        solver['f'].evaluate()
        
        solver['dg'].setInput(Primal, 0)
        solver['dg'].setInput(EP,     1)
        solver['dg'].evaluate()
        
        solver['g'].setInput(Primal, 0)
        solver['g'].setInput(EP,     1)
        solver['g'].evaluate()
    
        QP =       {
                    'H'  : DMatrix(solver[ 'H'].output()),
                    'f'  : DMatrix(solver[ 'f'].output()),
                    'dg' : DMatrix(solver['dg'].output()),
                    'g'  : DMatrix(solver[ 'g'].output()),
                    'lbX': DMatrix(  lbV.cat - Primal.cat),
                    'ubX': DMatrix(  ubV.cat - Primal.cat),
                    'lbg': DMatrix(solver['lbg'].cat - solver['g'].output()),
                    'ubg': DMatrix(solver['ubg'].cat - solver['g'].output())
                    }
        
        return QP
    
    def PlotSparsity(self, solver):
        """
        Call in with your T.Solver as argument.
        """
        solver['H'].evaluate()
        H = solver['H'].output()

        solver['dg'].evaluate()
        dg = solver['dg'].output()
        
        plt.figure(1)
        plt.spy(np.concatenate([H,dg.T],axis=1))
        plt.show()
    
    def _setParameters(self, EP = [], ExtParam = []):
        """
        Sets the FreeWay parameters. Call in with a Dictionary of your external parameters.
        Segment lengths are automatically assigned from your specified data file.
        """

        for key in ExtParam.keys():
            EP['Param','Global',key] = ExtParam[key]
            
        for i in range(self.Data['L'].shape[1]):
            EP['Param','Dist',i,'L'] = self.Data['L'][0,i]
        
        return EP
    
    def _AssignData(self, Simulated, Data = 0, VStruct = [], EPStruct = [], ExtParam = []):
        
        TimeRange = len(VStruct['States',:,...,veccat])
        
        EPStruct = self._setParameters(EP = EPStruct, ExtParam = ExtParam)
        
        SegmentList = range(self.NumSegment)            
            
        #Assign Measurements
        DataDic = {'v': ['vv_3lane', 1.], 'rho': ['rr_3lane', 3.]}
        for k in range(TimeRange):
            for i in SegmentList:
                for key in DataDic.keys():
                    EPStruct['Meas', k,key,i]         = Data[DataDic[key][0]][k, i]/DataDic[key][1]
                    VStruct['States',k,key,i]         = Data[DataDic[key][0]][k, i]/DataDic[key][1]
                    EPStruct['Param','States0',key,i] = Data[DataDic[key][0]][0, i]/DataDic[key][1]
        
        #Assign accident-free values
        VStruct['Inputs',:,'alpha',:] = 0.
        VStruct['Inputs',:,'beta',:] = 1.
        EPStruct['Param','Inputs0','alpha',:] = 0.
        EPStruct['Param','Inputs0','beta',:]  = 1.
         
        for key in ['rho','v']:
            EPStruct['Meas',:,key,:] = VStruct['States',:,key,:]
 
        #Assign Ve
        vfree  = EPStruct['Param','Global', 'vfree'  ]
        a      = EPStruct['Param','Global', 'a'      ]
        rho_cr = EPStruct['Param','Global', 'rho_cr' ]
        for k in range(TimeRange):
            for i in SegmentList:    
                Ve_arg   = (1 + VStruct['Inputs',k,'alpha',i])*VStruct['States',k,'rho',i]
                VStruct['Inputs',k,'Ve',i] = vfree*exp(-((Ve_arg/rho_cr)**a)/a)
                
        return VStruct, EPStruct
    
    def GenerateData(self, VStruct = [], EPStruct = [], ExtParam = [], TimeRange = 0, Simulated = False, AddNoise = False):
        """
        Generate data for testing purposes
        """
        print "Assign Data"
        StateNoise = 0

        Data, EP = self._AssignData(Simulated, Data = self.Data, VStruct = VStruct, EPStruct = EPStruct, ExtParam = ExtParam)
        
        if AddNoise:
            rhoMean = np.mean(Data['States',:,'rho',:])
            vMean   = np.mean(Data['States',:,'v',:])
        
        if TimeRange == 0:
            TimeRange = len(Data['States',:,...,veccat])
        
        if Simulated:
            print "Construct Simulated Data"
            for k in range(TimeRange-1):
                for i in range(1,self.NumSegment-1):
                         
                    for index, setInput in enumerate([Data['Inputs',k],Data['States',k],EP['Param']]):
                        self._Shoot[i].setInput(setInput,index)        
                    self._Shoot[i].evaluate()
                    
                    if AddNoise:
                        StateNoise = np.array([
                                                rand.normalvariate(0,1e-3*rhoMean),
                                                rand.normalvariate(0,1e-3*vMean)
                                              ])
                    Data['States',k+1,...,i] = list(self._Shoot[i].output() + StateNoise)

            for key in ['rho','v']:
                EP['Meas',:,key,:] = Data['States',:,key,:]
            
        print "Done"
        
        return Data, EP
    
    def PassMHEData(self,Data = [], EP = [], time = 0, ExtParam = []):
        """
        Assign data to the MHE solver
        """
        
        initMHE = self.V()
        EPMHE   = self._setParameters(EP = self.EP(), ExtParam = ExtParam)
        
        for key in ['States0','Inputs0']:
            EPMHE['Param',key] = EP['Param',key]
        
        for k in range(self.Horizon):
            initMHE['States',k,...]   =  Data['States',time+k,...]
            EPMHE['Meas',k,...,:] =  EP['Meas',time+k,...,:]
        
        for k in range(self.Horizon-1):
            initMHE['Inputs',k,...]   =  Data['Inputs',time+k,...]
            
        initMHE['Inputs',:,'alpha',:] = 0.
        initMHE['Inputs',:,'beta',:]  = 1.
        
        
        #Set variable bounds & initial guess
        lbV  = self.V(-inf)
        ubV  = self.V( inf)
        
        #Collapsed variables
        for key in ['f','r']:
            lbV['Inputs',:,key,:] = 0.
            ubV['Inputs',:,key,:] = 0.
        
        #Bounds on alpha, beta
        for key in ['alpha', 'beta']:
            lbV['Inputs',:,key,:] = 0
            ubV['Inputs',:,key,:] = 1
        
        lbV['Slacks'] = 0.
        lbV['States'] = 0.
        lbV['Inputs',:,'Ve',:] = 0.
        
        for i in range(self.NumSegment):
            for k in range(self.Horizon-1):
                VeBound = np.mean(self.Data['VMS_d'][time+k,3*i:3*i+3])
                if np.isnan(VeBound):
                    VeBound = 200.
                ubV['Inputs',k,'Ve',i] = VeBound
        
        return initMHE, EPMHE, lbV, ubV
    
    def SolveMHE(self, EP = [], lbV = [], ubV = [], init = []):
        print "Solve Instance of MHE"
        
        self.Solver['solver'].setInput(init,                   "x0")
        self.Solver['solver'].setInput(self.Solver['lbg'],    "lbg")
        self.Solver['solver'].setInput(self.Solver['ubg'],    "ubg")
        self.Solver['solver'].setInput(lbV,                   "lbx")
        self.Solver['solver'].setInput(ubV,                   "ubx")
        self.Solver['solver'].setInput(EP,                      "p")
            
        self.Solver['solver'].solve()
        
        self._lbg = self.Solver['lbg']
        self._ubg = self.Solver['ubg']
        
        Adjoints = self.g(np.array(self.Solver['solver'].output('lam_g')))
        Primal   = self.V(np.array(self.Solver['solver'].output('x')))
        Status   = int(self.Solver['solver'].getStats()['return_status'] == 'Solve_Succeeded')
        
        return Primal, Adjoints, Status
    
    def Shift(self, Primal, EP):
        
        #Assign the previous 1st state and input to EP
        EPShifted = EP
        EPShifted['Param','States0'] = Primal['States',0]
        EPShifted['Param','Inputs0'] = Primal['Inputs',0]
        
        #Primal Shift
        PrimalShifted = self.V(Primal)
        PrimalShifted[...,:-1] = Primal[...,1:]
        PrimalShifted[...,-1]  = Primal[...,-1]
        
        #Last step completion
        for i in range(1,self.NumSegment-1):
            for index, setInput in enumerate([Primal['Inputs',-1],Primal['States',-1],EP['Param']]):
                self._Shoot[i].setInput(setInput,index)        
            self._Shoot[i].evaluate()
            PrimalShifted['States',-1,...,i] = list(self._Shoot[i].output())
        
        #Assign measurement to the boundary conditions
        for i in [0,self.NumSegment-1]:
            PrimalShifted['States',-1,...,i] = EP['Meas',-1,...,i]
            
        return PrimalShifted
    

    def Check(self, init = 0, EP = 0):
        #Check initial guess
        self.Costfunc.setInput(init, 0)
        self.Costfunc.setInput(EP, 1)
        self.Costfunc.evaluate()
        Cost0 = self.Costfunc.output()
        
        self.gfunc.setInput(init, 0)
        self.gfunc.setInput(EP, 1)
        self.gfunc.evaluate()
        g0 = self.gfunc.output()
        
        return Cost0, self.g(g0)
        
    