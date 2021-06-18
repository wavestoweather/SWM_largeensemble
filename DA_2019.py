#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:35:31 2019

@author: Yvonne.Ruckstuhl
"""
import numpy as np
import numpy.linalg as lin
import pdb
import math
import scipy.linalg
from constants import *


def EnKF(obs,model,seed):
    obs_perturb = np.zeros((3*n,nens))
    np.random.seed(seed)
    for i in range(3):
        if i == 2:
            obs_perturb[i*n:(i+1)*n,:] = (np.random.lognormal(-8,sig[i],(n,nens)).T + obs[i*n:(i+1)*n]).T
        else:
            obs_perturb[i*n:(i+1)*n,:] = (np.random.normal(0,sig[i],(n,nens)).T + obs[i*n:(i+1)*n]).T
    P = np.cov(model)*LocMat(n,loc_radius)
    PH = P[:,H]
   # pdb.set_trace()
    S = np.dot(PH,lin.inv(PH[H,:]+np.diag(rdiag)))
    model = model + np.dot(S,obs_perturb[H,:]-model[H,:])
    model[2*n:3*n,:] = np.where(model[2*n:3*n,:]< 0, 0.0,model[2*n:3*n,:])
    return model

def etkf(obs,model,seed):
    
    obs = obs[H]
    synobs = model[H,:]
    synobsmean = np.average(synobs,axis = 1)
    synobs = synobs.T - synobsmean
    modelmean = np.average(model,axis = 1)
    model =  model.T - modelmean
   # pdb.set_trace()
    C = np.divide(synobs,rdiag) 
    pinverse = ((nens-1.)/1.)*np.identity(nens) + np.dot(C, synobs.transpose())
    p = lin.inv(pinverse)
    (U,S,VT) = lin.svd(p)
    D = np.diag(np.sqrt(S))
    w = np.sqrt((nens-1.))*np.dot(np.dot(U,D),VT)
    wmean = np.dot(p, np.dot(C, (obs-synobsmean)))
    w = w + wmean
    model = (np.dot(w,model) + modelmean).T
    model[2*n:3*n,:] = np.where(model[2*n:3*n,:]< 0, 0.0,model[2*n:3*n,:])
    return model

def letkf(obs,model,seed):           # seed redundant here but need for msw main script

      obs = obs[H]
      synobs = model[H,:].T
      model = model.T
      nobsnew=len(H)

      
      maxrad=int(2*loc_radius)
      wlen=np.minimum(2*maxrad+1,n)  # number of local observations used
      
      xx=np.arange(n)  # the gridpoints
      synobsmean = np.average(synobs,0)
      synobs = synobs - synobsmean
      # 2) forecast ensemble - average and perturbation     
      modelmean = np.average(model,0).reshape(-1)
      model =  np.reshape(model,(nens,-1)) - modelmean
    
      oldmodel = np.reshape(model,(nens,3,n))
      modelmean=np.reshape(modelmean,(3,n))
      model=np.reshape(model,(nens,3,n))
          
     

      for j in range(0,n):  # Loop over the analysis grid points
          #pdb.set_trace()
          window = np.zeros(3*wlen)
          window_temp=np.roll(xx,maxrad-j)
          window[0:wlen]=window_temp[0:wlen]
          window[wlen:2*wlen]=n+window_temp[0:wlen]
          window[2*wlen:3*wlen]=2*n+window_temp[0:wlen]

          lmodel=np.reshape(oldmodel[:,:,j],(nens,3))
          lmodelmean=np.reshape(modelmean[:,j],3)
          # Get the local variables. 

          [lsynobs,lrdiag,lsynobsmean,lobs] = winobs(j,n,synobs,synobsmean,obs,rdiag,oldmodel,H,window,gasparicohn(n,loc_radius),nobsnew)
         
          C = np.divide(lsynobs,lrdiag) # only works if R is diagonal, C has dimension k x nobs          
          pinverse = ((nens-1.)/1.0)*np.identity(nens) + np.dot(C, lsynobs.transpose())
          p = lin.inv(pinverse)

          try:                          # solving problem of non-convergence with numpy svd
            (U,S,VT) = scipy.linalg.svd(p)
          except:                       # here using the relation between a SVD and an eigenvalue decomposition
            (S,VT) = lin.eig(np.dot(p.T,p))
            (S,U) = lin.eig(np.dot(p,p.T))
            S = S**0.5

          D = np.diag(np.sqrt(S))
          w = np.sqrt((nens-1.))*np.dot(np.dot(U,D),VT)

          wmean = np.dot(p, np.dot(C, (lobs-lsynobsmean)))
          w = w + wmean
   
          model[:,:,j] = (np.dot(w, lmodel.reshape(nens,3)) + lmodelmean.reshape(-1)).reshape(nens,3)
              

      model = model.reshape(nens,-1).T      
      model[2*n:3*n,:] = np.where(model[2*n:3*n,:]< 0, 0.0,model[2*n:3*n,:])
      return model
      
      

def get_obs(truth,seed):
    n = np.int(len(truth)/3)
    obs = np.zeros((np.int(3*n)))            
    np.random.seed(seed)
    for i in range(3):
        if i == 2:
            obs[i*n:(i+1)*n] = truth[i*n:(i+1)*n] + (np.random.lognormal(-8,sig[i],n))
        else:
            obs[i*n:(i+1)*n] = truth[i*n:(i+1)*n] + np.random.normal(0,sig[i],n)
    return obs

def winobs(anpos,n,synobs,synobsmean,obs,rdiag,oldmodel,obspos,window,gcf,nobsnew):
    #determine the subset of those observations that lie in the current local window.


    locwindow=np.intersect1d(window,obspos) # chose the position of local obs
                                            # depending on radius and obspos
    maxlocw=len(locwindow)                                        

    # Now get the indices of all local observations in the big obs array.
    locobsposind=np.zeros(1).astype(int)

    for i in range(0,nobsnew):
        for j in range(0,maxlocw):
            if obspos[i]==locwindow[j]:
                locobsposind=np.append(locobsposind,i)
    #pdb.set_trace()
    locobsposind=locobsposind[1:len(locobsposind)]


    lsynobs=synobs[:,locobsposind]
    lrdiag=rdiag[locobsposind]

    lsynobsmean=synobsmean[locobsposind]
    lobs=obs[locobsposind]  # original: lobs=obs[locobsposind]
    
    # Calculate the distance and apply gaspari-cohn weights to the R-Matrix
    for i in range(0,len(lobs)):
        
        dis=np.remainder(abs(obspos[locobsposind[i]]-anpos),n)

        if dis>n/2:
            dis=abs(n-dis)
            
        lrdiag[i]=lrdiag[i]/gcf[dis]

    return lsynobs,lrdiag,lsynobsmean,lobs


def LocMat(nx,c):
    LocMat = np.zeros((nx,nx))

    for l in range(2*c):
        b = float(c)
        z = l/b
        for j in range(0,nx):
            if l < c:
                if l+j < nx:
                    #pdb.set_trace()
                    LocMat[j+l,j] = -0.25*math.pow(z,5)+0.5*math.pow(z,4)+(5.0/8.0)*math.pow(z,3)-(5.0/3.0)*math.pow(z,2)+1.0
                    LocMat[j,j+l] = LocMat[j+l,j]
                else:
                    #pdb.set_trace()
                    LocMat[l-nx+j,j] = -0.25*math.pow(z,5)+0.5*math.pow(z,4)+(5.0/8.0)*math.pow(z,3)-(5.0/3.0)*math.pow(z,2)+1.0
                    LocMat[j,l-nx+j] = LocMat[l-nx+j,j]
            else:
                if l+j < nx:
                    #pdb.set_trace()
                    LocMat[j+l,j] = (1.0/12.0)*math.pow(z,5)-0.5*math.pow(z,4)+(5.0/8.0)*math.pow(z,3)+(5.0/3.0)*math.pow(z,2) -5.0*z+4.0-(2.0/3.0)*(b/l)
                    LocMat[j,j+l] = LocMat[j+l,j]
                else:
                    #pdb.set_trace()
                    LocMat[l-nx+j,j] = (1.0/12.0)*math.pow(z,5)-0.5*math.pow(z,4)+(5.0/8.0)*math.pow(z,3)+(5.0/3.0)*math.pow(z,2) -5.0*z+4.0-(2.0/3.0)*(b/l)
                    LocMat[j,l-nx+j] = LocMat[l-nx+j,j]
    LocMat=np.tile(LocMat,(3,3))
    return LocMat
def gasparicohn(n,irad):
    gcf=np.zeros(n)
    c=float(irad)
    for i in range(n):
        z=float(i)
#        print i,z
        if i <= irad:
            gcf[i]=-0.25*(z/c)**5+0.5*(z/c)**4+(5.0/8.0)*(z/c)**3-(5.0/3.0)*(z/c)**2+1.0
        elif i <= 2*irad:
            gcf[i]=(1.0/12.0)*(z/c)**5-0.5*(z/c)**4+(5.0/8.0)*(z/c)**3+(5.0/3.0)*(z/c)**2-5*(z/c)+4-(2.0/3.0)*(c/z)
        else:
            gcf[i]=0.0;
    
#    print gcf[0:10]
    return gcf
