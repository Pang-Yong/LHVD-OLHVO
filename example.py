# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 20:16:59 2022

@author: pang
"""

import numpy as np
from OLHVD import olhvd
import matplotlib.pyplot as plt







def setting(name):

    if name=='example1':
        def example1_1(x):
            #print(x)
            y=-(x[1]-1)**2-(x[0]-1)**2+1
            if y<0:
                return 0
            else:
                return 1
            
        lowb=np.array([0,0])
        upb=np.array([1,1])    
        return lowb,upb, [example1_1]
        
        
        
    if name=='example2': 
        def example2_1(x):
            #print(x)
            y=min([max([0.5-x[0],0.5-x[1]]),(x[0]**2+x[1]**2-0.25)])
            if y<0:
                return 0
            else:
                return 1
              
        lowb=np.array([0,0])
        upb=np.array([1,1])    
        return lowb,upb, [example2_1]


if __name__ == '__main__':
    


    e='example1'
    lowb,upb,const=setting(e)
    n_points=20
    n_dim=2
    
    
    lhd=olhvd(n_points,n_dim,lowb,upb,const,opt=True)  
    result=lhd.run()
    mm_samples=lhd.select_samples
       
    plt.figure(2)
    plt.scatter(mm_samples[:,0],mm_samples[:,1])
    plt.scatter(result[:,0],result[:,1])
