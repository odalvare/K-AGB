# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 17:35:17 2023

@author: oscar
"""

import sys
import platform
import gc

import rasterio
import time

#Detect operaring system and apply folder delimiter
ost=platform.system()
delim='/'
if(ost=='Windows'):
    delim='//'

import numpy as np
import pandas as pd

from scipy import interpolate
from scipy import ndimage

class MonteCarloBiomass:
        
    #Begining of constructor --------------------------------------------------
    def __init__(self,project,params,foldername):
        """
        Construction of an instance for reading K-AGB basic execution parameters

        Parameters:
        ----------
            project : string
                The folder location of the project
            parameters : dictionary
                Basic parameters of the project
        
        Returns:
        -------
            null
            
        """
        
        #Constructor begins here
        
        self.__workspace=project
        self.__params=params
        self.__folder=foldername
        
        # Open biomass probability distribution        
        rtcpf=self.__workspace+'data'+delim+'montecarlo'+delim+'Biomass_cpf_cond.csv'
        self.__dfcpf=pd.read_csv(rtcpf,sep=',')
        
        self.__categories=self.__params['categories']
        self.__ncat=self.__params['n_categories']
        self.__years=self.__params['years']
        self.__nyears=self.__params['n_years']
        self.__names=['landcover_'+str(year) for year in self.__years] 
        self.__nreas=self.__params['mc_reas']
        self.__coverreg=[c+1 for c in range(self.__ncat)]
        
        self.__metacover={}
        
    #End of constructor -------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def run_mcSimulations(self):
        """
        Performs the biomass MonteCarlo simulations

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            names : list
                array of names
        """
            
        #code begins here
        
        # configure the interpolators of the biomass probability distributions
        back_funcs=[]
        for c in range(self.__ncat):            
            p=self.__dfcpf['Prob.'].values
            perc=self.__dfcpf[self.__categories[c].decode("utf-8")].values
            back_funcs.append(interpolate.interp1d( p/100.0,
                                                    perc,
                                                    kind='linear',
                                                    fill_value="extrapolate"))
        
        # Run all target years ------------------------------------------------
        
        for k in range(self.__nyears):
            
            t = time.time()
            t1 = time.localtime()
            print('Simulating year '+str(self.__years[k])+' -----------------')
            print('Starting time '+str(time.strftime("%H:%M:%S", t1)))
        
            loc=self.__workspace+'data'+delim+'montecarlo'+delim+'landcover'+delim+self.__folder+delim
            rt=loc+self.__names[k]+'.tif'
            
            cover=rasterio.open(rt)
            self.__metacover=cover.meta.copy()
            
            nxcover=cover.width                     # Number of columns
            nycover=cover.height                    # Number of rows
            rxcover,rycover=cover.res
            
            arrcover=cover.read(1)
            arrcover=arrcover.astype(np.float32)
            arrcover[arrcover[:]<=0]=np.nan
            
            # Create an array with random numbers uniformly distributed
            
            coverdata=np.reshape(arrcover,(nxcover*nycover))
            coverdata=np.tile(coverdata, (self.__ncat, 1)).T
            covercheck=np.tile(self.__coverreg, (nxcover*nycover, 1))
            
            biomassen=np.zeros((self.__nreas,nycover,nxcover),dtype=np.float32)
            
            mask=coverdata==covercheck
            
            probsim=np.random.rand(self.__nreas,nycover*nxcover,self.__ncat)
            simbio=np.zeros_like(probsim)
            
            for rea in range(self.__nreas):
                
                print('Simulation '+str(rea+1))
                
                probsim=np.random.rand(nycover*nxcover,self.__ncat)
                simbio=np.zeros_like(probsim)
                
                for c in range(self.__ncat):
                    bf=back_funcs[c]
                    simbio[:,c]=bf(probsim[:,c])
                
                mask=np.where(mask,1,0)
                simbio1=simbio*mask
                simbio1=simbio1.sum(axis=1)
                mask1=np.where(np.isnan(arrcover),np.nan,1.0)
                simbio2=np.reshape(simbio1,(nycover,nxcover))*mask1
                simbio2=ndimage.percentile_filter(simbio2, percentile=50, size=5)
                biomassen[rea,:,:]=simbio2
            
            # ...
            t2 = time.localtime()
            print('Ending time '+str(time.strftime("%H:%M:%S", t2)))
            print('Biomass Montecarlo simulation executed')
            print(np.round_(time.time() - t, 3), 'sec elapsed')
            
            loc=self.__workspace+'output'+delim
            rt=loc+'Biomass_sim_'+str(self.__years[k])+'.tif'
                
            # Write to TIFF
            kwargs = cover.meta
            kwargs.update(
                dtype=rasterio.float32,
                count=self.__nreas,
                compress='lzw')
            
            with rasterio.open(rt, 'w', **kwargs) as dst:
                dst.write(biomassen)
                
            cover.close()
        
        del simbio,simbio1,simbio2,mask,mask1,coverdata,covercheck,probsim,bf,arrcover        
        del back_funcs
        
        gc.collect()
        
        return biomassen
    
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def getNames(self):
        """
        Returns the landcover names

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            names : list
                array of names
        """
            
        #code begins here
        return self.__names
    
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def getMetaCover(self):
        """
        Returns the landcover names

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            names : list
                array of names
        """
            
        #code begins here
        return self.__metacover
    
    #End of method ------------------------------------------------------------
            

            
            
    