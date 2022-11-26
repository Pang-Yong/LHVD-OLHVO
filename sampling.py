# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 20:18:22 2021

@author: pangyong
"""
import numpy as np
import matplotlib.pyplot as plt
import copy



class mm_sampling:
    
    """
        Conducting Monte Carlo sampling in the design space and dividing the samples into 
    feasible subset and infeasible subset
    """
    
    
    def __init__(self,n_dim, m,lb,ub,constraint_ueq=tuple()):
        self.lb=lb
        self.ub=ub
        self.has_constraint = len(constraint_ueq) > 0
        if self.has_constraint:
            self.constraint_ueq=list(constraint_ueq)
        self.n_dim=n_dim
        self.m=m
        
        
    def mmsampling(self ):
        """
            MM sampling

        Returns
        -------
        samples : array
            MM samples

        """
        sampling_num=self.m
        samples=np.random.rand(sampling_num,self.n_dim) 
        #plt.scatter(samples[:,0],samples[:,1])
        return samples
    
    
    def select(self,samples):
        for i in range(samples.shape[0]-1,-1,-1):
            #print(samples.shape)
            index=[self.constraint_ueq[j](samples[i,:]*(self.ub-self.lb)+self.lb) for j in range(len(self.constraint_ueq))]
            #print(index)
            label= np.any(np.array(index)==1)
            if label:
               samples=np.delete(samples,i,axis = 0)
                
            
        return samples
        
        
    def run(self):
        mmsamples=self.mmsampling()
        select_samples=self.select(copy.deepcopy(mmsamples))
        return select_samples
    
    
    
    
    
            
class  lhvd_sampling:
    """
        generate design points based on the MM samples
    
    
    
    """
    def __init__(self,mm_samples):
        self.mm_samples=mm_samples
        self.n_mmsamples=self.mm_samples.shape[0]
        self.n_dim=self.mm_samples.shape[1]
        
        
        
        
    def run(self,n_lhsamples):
        
        
        
        
        
        n_eachcube=int(self.n_mmsamples/n_lhsamples)
        remainder=self.n_mmsamples%n_lhsamples
        lhsamples=np.zeros((n_lhsamples,self.n_dim))
        
        for i in range(self.n_dim):#对每个维度
            temp_mmsamples=copy.deepcopy(self.mm_samples)
            temp_mmsamples=temp_mmsamples[temp_mmsamples[:,i].argsort()]
            for j in range(n_lhsamples):
                 if j <remainder:
                     temp_x=j*(n_eachcube+1)+int(np.ceil((n_eachcube+1)/2))
                     lhsamples[j,i]=temp_mmsamples[temp_x-1,i]
                 else:
                     temp_x=j*(n_eachcube)+int(np.ceil((n_eachcube)/2))+remainder
                     lhsamples[j,i]=temp_mmsamples[temp_x-1,i]
        
            shuffleindex=np.random.permutation(np.arange(n_lhsamples))   
            lhsamples[:,i]=lhsamples[:,i][shuffleindex]

            
        return  lhsamples
             
              
                
            
        
        
        
        
          
            
            
            
            
            
            
            
            
if __name__ == '__main__':
    #建立约束，满足约束返回1，不满足返回0
    def cons1(x):
        #print(x)
        y=x[1]-x[0]
        if y<0:
            return 1
        else:
            return 0
        
    
    lb=np.array([0,0])
    ub=np.array([1,2])
    samp=mm_sampling(2,10000,lb,ub,constraint_ueq={cons1})
    select_samples=samp.run()
    plt.scatter(select_samples[:,0],select_samples[:,1])
    
    
    lhdsamp=clhd_sampling(select_samples)
    lhd=lhdsamp.run(10)
    plt.scatter(lhd[:,0],lhd[:,1])
    
    
    