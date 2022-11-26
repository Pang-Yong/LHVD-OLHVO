# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 20:26:23 2021

@author: pangyong
"""

import numpy as np
import random


class csa:
    """
    constrained simulated annealing
    
    parameter
    lhd:array
        inital design 
    
    
    """
    def __init__(self,lhd, lb,ub,con_opt,opt,Imax,rate_intial,rate_stopping,theta):
        self.H=lhd
        self.n_sam=lhd.shape[0]
        self.n_dim=lhd.shape[1]
        self.lb=lb
        self.ub=ub
        self.con_opt=con_opt
        self.opt=opt
        self.rate_intial=rate_intial
        self.rate_stopping=rate_stopping
        self.Imax=Imax
        self.theta=theta        
            
            
            
            
    


    def randomexchange(self,samples):
        """
        completely random exchange for new design
        
        """
        
        samples2=samples.copy()
        
        column_index=random.sample(range(self.n_dim),1)       
        row_index_couple=random.sample(range(self.n_sam),2)
        #print(column_index,row_index_couple)
        samples2[row_index_couple[1],column_index]=samples[row_index_couple[0],column_index]
        samples2[row_index_couple[0],column_index]=samples[row_index_couple[1],column_index]
        #print("samples",samples,"samples2",samples2)
        return samples2,row_index_couple
    

    
    def cal_con(self,H,r,row):
        """
        calculate S given the index of the elements for exchanging

        """
        R=r.copy()
        for i in row:
            R[i]=0
            for j in range(len(self.con_opt)):
                R[i]=R[i]+self.con_opt[j](H[i,:]*(self.ub-self.lb)+self.lb)

        return R
        
    
    def cal_con1(self,H): 
        """
        calculate S

        """
        R=np.zeros(H.shape[0])
        for i in range(H.shape[0]):
            for j in range(len(self.con_opt)):
                R[i]=R[i]+self.con_opt[j](H[i,:]*(self.ub-self.lb)+self.lb)
        return R
        
        
        
        
        
        
    def cal_obj(self,H):
        
        """
        Calculates the space-filling criterion 
        from an open toolbox:pykring http://www.pykriging.com/
        

        """
        p=2
        q=2
        #calculate the distances between all pairs of
        #points (using the p-norm) and build multiplicity array J
        J,d = self._jd(H,p)
        #print(H)
        #the sampling plan quality criterion
        
        Phiq = (np.sum(J*(d**(-q))))**(1.0/q)
        return Phiq
    
        
    
    def _jd(self, H,p=1):
        ######p

        #number of points in the sampling plan
        n = H.shape[0]
        #computes the distances between all pairs of points
        d = np.zeros((n*(n-1)//2))
    
    
    
    #    for i in xrange(n-1):
    #        for j in xrange(i+1,n):
    #            if i == 0:
    #                d[i+j-1] = np.linalg.norm((rld[0,:]-rld[j,:]),2)
    #            else:
    #                d[((i-1)*n - (i-1)*i/2 + j - i  )] = np.linalg.norm((H[i,:] - H[j,:]),2)
    
        #an alternative way of the above loop
        list = [(i,j) for i in range(n-1) for j in range(i+1,n)]
        for k,l in enumerate(list):
            d[k] = np.linalg.norm((H[l[0],:]-H[l[1],:]),p)
    
    
        #remove multiple occurences
        distinct_d, J = np.unique(d, return_counts=True)
    
        return J, distinct_d    



    
    def initialization(self):
        """
        initialization of the parameters in constrained simulated annealing
        """

        self.max_iter=self.Imax
        #设置初始温度    
        H = self.H.copy()
        Phi_H=self.cal_obj(H)
        temp=0
        count=0
        for i in range(self.max_iter):
            H_try,row_index=self.randomexchange(H)
            Phi_try=self.cal_obj(H_try)
            if Phi_try>Phi_H:
                count=count+1
                temp=temp+Phi_try-Phi_H
            H=H_try
            Phi_H=Phi_try         
        ave_cost=temp/count
        self.T=(-ave_cost)/(np.log(self.rate_intial))
        print("T=",self.T)
        
        
        
    
    
    
    
    def compare(self,H,H_try,R,R_try):
        """
        update the current design
        """
        
        n_R=np.sum(R)
        n_try=np.sum(R_try)
        if n_try<n_R:
            self.count=self.count+1
            return H_try,R_try
        elif n_try==n_R:
            Phi_H=self.cal_obj(H)
            Phi_try=self.cal_obj(H_try)
            if Phi_try<=Phi_H:
                self.count=self.count+1
                return H_try,R_try
            else:                              
                if random.random()<np.exp(-(Phi_try-Phi_H)/self.T):
                    self.count=self.count+1
                    return H_try,R_try
                else:
                    return H,R     
        else:
            return H,R
        
    
    def comparebest(self,H,H_try,R,R_try):
        """
        update the best design
        """        
        n_R=np.sum(R)
        n_try=np.sum(R_try)
        if n_try<n_R:
            return H_try,R_try
        elif n_try==n_R:
            Phi_H=self.cal_obj(H)
            Phi_try=self.cal_obj(H_try)
            if Phi_try<Phi_H:
                return H_try,R_try
            else:                              
                return H,R     
        else:
            return H,R    


         

    def con_exchange1(self,samples,c_index):
        """ 
        selective random exchange
        """
        
        samples2=samples.copy()
        column_index=random.sample(range(self.n_dim),1)       
        row_index_couple=random.sample(list(c_index),1)
        remain=np.delete(np.array(range(self.n_sam)),row_index_couple)
        row_index_couple.extend(random.sample(list(remain),1))
        
        
        
        #print(column_index,row_index_couple)
        samples2[row_index_couple[1],column_index]=samples[row_index_couple[0],column_index]
        samples2[row_index_couple[0],column_index]=samples[row_index_couple[1],column_index]
        #print("samples",samples,"samples2",samples2)
        
        return samples2,row_index_couple   






    def accel(self,H,R,H_best,R_best):
        """
        feasibility acceleration
        
        
        """
        for I in range(self.max_iter):
            c_index=np.array(np.where(R>0)).flatten()
            if c_index.shape[0]<1:
                break
            else:
                H_try,row_index=self.con_exchange1(H,c_index)
                R_try=self.cal_con(H_try,R,row_index)
                H,R=self.compare(H,H_try,R,R_try)
                                              
                H_best,R_best=self.comparebest(H_best,H_try,R_best,R_try) 
                print(np.sum(R_best))
                    
        return H,R,H_best,R_best






    
    def optimize(self):
        """
            conducting constrained simulated annealing
    
        """       
              
        H = self.H.copy()
        R=self.cal_con1(H) 
        H_best=H.copy()
        R_best=R.copy()
        self.initialization()
        self.count=0
        H,R,H_best,R_best=self.accel(H,R,H_best,R_best)
        if np.sum(R_best)>0:
            print("warning:acceleration module cannot find feasible design,please make sure the correctness of the setting")
        
        
        if self.opt:
            while 1:
                self.count=0
                for i in range(self.Imax):
                    H_try,row_index=self.randomexchange(H)
                    R_try=self.cal_con(H_try,R,row_index)
                    H,R=self.compare(H,H_try,R,R_try)
                    
                                  
                    H_best,R_best=self.comparebest(H_best,H_try,R_best,R_try) 
                    
        
                accept_rate=self.count/self.max_iter    
                print("Accept rate",accept_rate)
                    
                if ( accept_rate<self.rate_stopping):
                    if np.sum(R_best)==0:
                        print("successfully ")
                        break
                    else:
                        print("fail to find feasible design")
                        break
                
                self.T=self.T*self.theta
                print("T=",self.T)
            
        return H_best
        
        
    
    
    
    