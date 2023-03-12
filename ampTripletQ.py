#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 17:47:38 2018

@author: xin liu at Stanford
liuxin@stanford.edu
liuxine@gmail.com


Liu, X., Ben-Zion, Y., & Zigone, D. (2015). Extracting seismic attenuation coefficients from cross-correlations of ambient noise at linear triplets of stations, Geophys. J. Int., 203, 1149–1163, 
https://doi.org/10.1093/gji/ggv357
All rights reserved
"""

import numpy as np
import scipy
import cmath,math
import scipy.signal
import matplotlib.pyplot as plt

import multiprocessing as mp
from scipy import interpolate


class ampInvTriplet(object):
    def __init__(self,freqvals,dictAmp=None,hasAmpError=True):
        """
        self.freqvals: contains frequency values used in triplet inversion
        dictAmp: a python dictionary that contains c12 and c13  keys
        """
        self.freqvals = freqvals # A SIMPLE NUMPY 1D ARRAY
        self.nfreqs = len(freqvals)
        self.hasAmpError =hasAmpError
        if dictAmp is not None:
            self.setAmpDicts(c12=dictAmp['c12'],c13=dictAmp['c13'],
                         c23=dictAmp['c23'])
    
    def setAmpDicts(self,c12=None,c13=None,c23=None):
        """
        The following DICTIONARIES:
        c12, c13 or c23, contains amparr, stdamparr, phaseslowarr (s/m) keys
        each of the keys has the same length as self.freqvals;
        c12 (and c13, c23) also contains 'distance' value in meter
        #
            amparr: amplitude versus frequency (1D array)
            stdamparr: amplitude error versus frequency (1D array) - if no error data, set values to 1
            dtphasearr: phase delay time of current station pair        
        """
        self.c12=c12
        self.c13=c13
        self.c23=c23
        
    def measureAmpDecay(self,ixcorr=2,ifgetsiteR=True,weightedLSQ=True,
                        applyGeometricCorrection=True,distratio_13_12=None):
        """
        ixcorr denotes the pair we derive for Q inversion:
            0==c12; 
            1==c13; 
            2==c23 (default)
        """
        aa=5
        if ixcorr==2:
            aa=5
            # for phase slowness/Q between C23
            #"""
            # G matrix for least-square inversion:
#            if ifgetsiteR:
#                G=np.concatenate(- (self.c23['distance']*np.pi*self.freqvals*self.c23['phaseslowarr'] ).T, np.ones(self.nfreqs,1) )
#            else: # BELOW: ASSUME log(site_ratio) TO BE ZERO
#                G = - (self.c23['distance']*np.pi*self.freqvals*self.c23['phaseslowarr'] ).T             
            # when dist_23 << dist_12 or dist_13
            diff23 = -np.pi*self.freqvals*(self.c13['dtphasearr'] - self.c12['dtphasearr']) 
            if ifgetsiteR:
                G=np.concatenate( ( diff23[:,None], np.ones( (self.nfreqs,1)) ),axis=1 )
            else: # BELOW: ASSUME log(site_ratio) TO BE ZERO
                G =  diff23[:,None]             
            # yerr IS THE STD ERROR FOR THE \delta C = ln(C13) - ln(C12):
            yerr=np.sqrt( (self.c12['stdamparr'])**2 / (self.c12['amparr'])**2 + 
                            (self.c13['stdamparr'])**2 / (self.c13['amparr'])**2)
            Gw= G/yerr[:, np.newaxis]  #bsxfun(@times,1./yerr{ixcorr},G);
            # yw: weighted data yvec
            # yvec: unweighted data vector
            if applyGeometricCorrection:
                if distratio_13_12 is not None:
                    yw=(0.5*np.log(distratio_13_12)+np.log(self.c13['amparr'])-np.log(self.c12['amparr'] ))/yerr;
                else:
                    yw=(np.log(self.c13['amparr']/np.sqrt(1.0/self.c13['dtphasearr']))-np.log(self.c12['amparr']/np.sqrt(1.0/self.c12['dtphasearr']) ))/yerr;
            else:
                yw=(np.log(self.c13['amparr'])-np.log(self.c12['amparr']))/yerr;
            yvec=yw*yerr
            indnan=np.isnan(yvec)
            if sum(~indnan) < 1:
                return None,None,None
            if weightedLSQ:
                # weighted least-square solution:
#                print yw
                beta2,resid,rank,s = np.linalg.lstsq(Gw[~indnan,:],yw[~indnan])
            else:
                # unweighted least-square solution:
                beta2,resid,rank,s = np.linalg.lstsq(G[~indnan,:],yvec[~indnan])
#            beta2=Gw\yw;
#            beta2=G\yvec{ixcorr};
            logampapprox2=np.dot(G,beta2)
#            covmat=np.linalg.inv(np.dot( Gw.T,Gw)  )
            covmat=np.linalg.inv(np.dot( Gw[~indnan,:].T,Gw[~indnan,:])  )
            Qinverr=covmat[0,0]**.5
            Qinv=beta2[0]
            if ifgetsiteR:
                logsiteR=beta2[-1]
            else:
                logsiteR= 0 #
            # NOTE: logsiteR can be zero if we force it to be that
            return Qinv,Qinverr,logsiteR
        if ixcorr==1:
            aa=5
            # for phase slowness/Q between C13
            #"""
            # G matrix for least-square inversion:
#            if ifgetsiteR:
#                G=np.concatenate(- (self.c23['distance']*np.pi*self.freqvals*self.c23['phaseslowarr'] ).T, np.ones(self.nfreqs,1) )
#            else: # BELOW: ASSUME log(site_ratio) TO BE ZERO
#                G = - (self.c23['distance']*np.pi*self.freqvals*self.c23['phaseslowarr'] ).T             
            # when dist_23 << dist_12 or dist_13
            diff13 = -np.pi*self.freqvals*(self.c23['dtphasearr'] + self.c12['dtphasearr']) 
            if ifgetsiteR:
                G=np.concatenate( ( diff13[:,None], np.ones( (self.nfreqs,1)) ),axis=1 )
            else: # BELOW: ASSUME log(site_ratio) TO BE ZERO
                G =  diff13[:,None]             
            # yerr IS THE STD ERROR FOR THE \delta C = ln(C13) - ln(C12):
            yerr=np.sqrt( (self.c12['stdamparr'])**2 / (self.c12['amparr'])**2 + 
                            (self.c23['stdamparr'])**2 / (self.c23['amparr'])**2)
            Gw= G/yerr[:, np.newaxis]  #bsxfun(@times,1./yerr{ixcorr},G);
            # yw: weighted data yvec
            # yvec: unweighted data vector
            if applyGeometricCorrection:
                if distratio_13_12 is not None:
                    yw=(0.5*np.log(distratio_13_12)+np.log(self.c23['amparr'])-np.log(self.c12['amparr'] ))/yerr;
                else:
                    yw=(np.log(self.c23['amparr']/np.sqrt(1.0/self.c23['dtphasearr']))-np.log(self.c12['amparr']/np.sqrt(1.0/self.c12['dtphasearr']) ))/yerr;
            else:
                yw=(np.log(self.c23['amparr'])-np.log(self.c12['amparr']))/yerr;
            yvec=yw*yerr
            indnan=np.isnan(yvec)
            if sum(~indnan) < 1:
                return None,None,None
            if weightedLSQ:
                # weighted least-square solution:
#                print yw
                beta2,resid,rank,s = np.linalg.lstsq(Gw[~indnan,:],yw[~indnan])
            else:
                # unweighted least-square solution:
                beta2,resid,rank,s = np.linalg.lstsq(G[~indnan,:],yvec[~indnan])
#            beta2=Gw\yw;
#            beta2=G\yvec{ixcorr};
            logampapprox2=np.dot(G,beta2)
#            covmat=np.linalg.inv(np.dot( Gw.T,Gw)  )
            covmat=np.linalg.inv(np.dot( Gw[~indnan,:].T,Gw[~indnan,:])  )
            Qinverr=covmat[0,0]**.5
            Qinv=beta2[0]
            if ifgetsiteR:
                logsiteR=beta2[-1]
            else:
                logsiteR= 0 #
            # NOTE: logsiteR can be zero if we force it to be that
            return Qinv,Qinverr,logsiteR
        
            #"""
#            np.cosh()

        

if __name__ == "__main__":
    """
    HERE I ONLY WANT TO TEST THE ampInvTriplet CLASS
    """  
    readsegy=True
    basedir='/Users/xinliu/Dropbox/Programing/python'
    fname='hahah.txt'
    ##
    # freqvals=np.arange([0.5, 0.6,0.7]) # frequencies
    nfreqs = 10
    freqvals=np.linspace(0.15,0.25, num=nfreqs)
    ampinvobj=ampInvTriplet(freqvals)
    ##
    """
        c12 (and c13, c23) also contains .distance value in meter
        #
        amparr: amplitude versus frequency (1D array)
        stdamparr: amplitude error versus frequency (1D array) - if no error data, set values to 1
        dtphasearr: phase delay time of current station pair
    """
    Qm = 25 #35 # modelled Q
    phasevel = np.ones(nfreqs)*2.0*1e3 # phase velocity in m/s
    phaseslow = 1./phasevel  #phase  slowness (s/m) versus frequency (1D array)
    stdamparr= np.ones(nfreqs)
    # interstation distances
    dist12 = 30*1e3 # in meter
    dist13 = 35*1e3 #
    dist23 = dist13-dist12
    amparrc12 = np.ones(nfreqs)/np.sqrt(dist12)
    amparrc13 = np.exp(-2*np.pi*freqvals*dist23*phaseslow/(2*Qm))/np.sqrt(dist13)
    ab=5
    #
    c12 = dict(distance=dist12, amparr=amparrc12, stdamparr=stdamparr,dtphasearr=dist12*phaseslow)
    c13 = dict(distance=dist13, amparr=amparrc13, stdamparr=stdamparr,dtphasearr=dist13*phaseslow)
    #
    ampinvobj.setAmpDicts(c12=c12,c13=c13)
    Qinv,Qinverr,logsiteR = ampinvobj.measureAmpDecay(ixcorr=2,ifgetsiteR=False)
    print('inverted Q value %.5f'% (1.0/Qinv))
    ab=6