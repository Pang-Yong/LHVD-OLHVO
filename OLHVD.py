# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 20:25:00 2021

@author: pangyong: pangy@mail.dlut.edu.cn
"""
import numpy as np
from sampling import mm_sampling,lhvd_sampling
from sa import csa

class olhvd:
    """
    lhvd and olhvd design methods for DoE in constrained desgin space.
    lhvd can generate only feasible design.
    olhvd can generate feasible design with optimal space-filling.
      
    parameters
    n_points: int
        number of the design points
    n_dim:int
        number of the dimenisons
    lb: array
        low bound of the design space
    ub: array
        up bound of the design space
    Imax: int
        itrations under each temperature, default: n_dim*n_points*(n_points-1)/2
    opt: bool
        True:olhvd
        False:lhd
    con_opt: list of functions
        The constraints of the design space，the function in the list is
        paramters
        x: array
            a design points from the design
        return       
        0: if design point satisfies the constraint
        1: if not
        
    con_MM: list of functions
        the constraints of MC sampling, equal to con_opt if it is not given.the type of 
        function is the same as con_opt
    rate_intial: float
        intial accept rate set in simulated annealing,default: 0.9
    rate_stopping: float
        stopping accept rate set in simulated annealing ,default: 0.01
    Alpha：int
        samples size of MC for representing the feasible design space,default: 1E4
    theta: int
        cooling coefficient of simulated annealing,default:0.9
    
    
    """
    
    
    def __init__(self, n_points, n_dim, lb, ub,con_opt,con_MM=[],Imax=0,opt=True,rate_intial=0.9,rate_stopping=0.01, Alpha=10000,theta=0.9 ):
        self.n_points=n_points
        self.n_dim=n_dim
        self.lb=lb
        self.ub=ub
        self.rate_intial=rate_intial
        self.rate_stopping=rate_stopping
        self.Alpha=Alpha
        self.theta=theta
        self.opt=opt
        if Imax==0:
            self.Imax=int(n_dim*n_points*(n_points-1)/2)
        else:
            self.Imax=Imax        
        self.con_opt=con_opt
        if len(con_MM) > 0:
            self.con_MM=con_MM
        else:
            self.con_MM=con_opt
    
                    
        
    def run(self):
        """
        

        Returns
        -------
        result : array
            desgin from LHVD or OLHVD

        """
        samp=mm_sampling(self.n_dim,10000,self.lb,self.ub,self.con_MM)
        self.select_samples=samp.run()
        if self.select_samples.shape[0]<self.Alpha:
            print(self.select_samples.shape[0])
            Nm=int(10000*self.Alpha/self.select_samples.shape[0])
            samp2=mm_sampling(self.n_dim,Nm,self.lb,self.ub,self.con_MM)
            self.select_samples=samp2.run()        
          
        if self.select_samples.shape[0]==0:
            return("no feasible region due to too many constraints ")
        else:
               
            lhdsamp=lhvd_sampling(self.select_samples)
            lhd=lhdsamp.run(self.n_points)
            sa1=csa(lhd, self.lb,self.ub,self.con_opt,self.opt,self.Imax,self.rate_intial,self.rate_stopping,self.theta)
            result_n=sa1.optimize()
            
            result=result_n*(self.ub-self.lb)+self.lb
            
            return result
        
        
        
        
        
        
        
        
        
        
        
        

  

