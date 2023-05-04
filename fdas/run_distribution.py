# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 09:06:46 2023

@author: oscar

Calculates the biomass cumulative probabily distribution for the secondary
information for land cover categories of interest

"""

import sys
import platform
import gc

import rasterio

#Detect operaring system and apply folder delimiter
ost=platform.system()
delim='/'
if(ost=='Windows'):
    delim='//'

import numpy as np
import pandas as pd

from downscaling.downscaling_help import resample_ras_int

class ProbabilityDistributions:
    
    #Begining of constructor --------------------------------------------------
    def __init__(self,project,params,refyr):
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
        self.__refyr=refyr
        
        # Opend simulated biomass (250 m) -------------------------------------
        
        rt=self.__workspace+'data'+delim+'fdas'+delim+'biosim_down.tif'
        
        biomass=rasterio.open(rt)        
        
        rxbio,rybio=biomass.res
        nreas=biomass.count
        
        biomass_bnds = biomass.bounds
        minxbio=biomass_bnds[0]
        minybio=biomass_bnds[1]
        maxxbio=biomass_bnds[2]
        maxybio=biomass_bnds[3]    
        
        extentbio = [minxbio,maxxbio,minybio,maxybio]
        
        nreas=200
        bands=[i+1 for i in range(nreas)]
        
        self.__arrbio=biomass.read(bands).astype(np.float32)        
        self.__names=['landcover_'+str(year) for year in params['years']]
        
        # Open land cover and resample-----------------------------------------
        
        rt=self.__workspace+'data'+delim+'fdas'+delim+'history_landcover'+delim+self.__names[self.__refyr]+'.tif'
        coveroriginal=rasterio.open(rt)
        print(coveroriginal.meta)        
        
        rt=self.__workspace+'data'+delim+'fdas'+delim+self.__names[self.__refyr]+'_res.tif'
        resample_ras_int(coveroriginal,                                        # origin cover (fine) - target biomass (coarse)
                         biomass,
                         rt,extentbio,
                         'EPSG:3116') 
        coveroriginal.close()
        self.__cover = rasterio.open(rt)
        
    #End of constructor -------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def estimate_cpd(self):
        """
        Estimates the prior biomass probability distribution for each land
        cover category

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            arrbiodown : numpy float array (nsim,nrow,ncols)
                Downscaled biomass realizations (from SGSIM)
            
        """
        
        nxcover=self.__cover.width
        nycover=self.__cover.height
        rxcover,rycover=self.__cover.res       
        
        arrcover=self.__cover.read(1)
        arrcover=arrcover.astype(np.float32)
        arrcover[arrcover[:]<=0]=np.nan
        
        arrcoveren=np.zeros((self.__params['sgsim_reas'],nycover,nxcover),dtype=np.float32)
        
        for rea in range(self.__params['sgsim_reas']):
            arrcoveren[rea,:,:]=arrcover[:,:]
        
        # Calculate probability functions -------------------------------------
        
        coverreg=[i+1 for i in range(self.__params['n_categories'])]
        print(coverreg)
        
        covernames=[self.__params['categories'][i].decode("utf-8") for i in range(self.__params['n_categories'])]
        print(covernames)
        
        ncoverreg=len(coverreg)
        
        cont=0
        
        # Parameters for the probability distribution estimation
        p = np.linspace(0, 100, 101)
        expdist=np.zeros((101,ncoverreg),dtype=np.float32)
        
        # Output
        index=[i for i in range(101)]
        columns=covernames
        covernames.append('Prob.')        
        dfcpf=pd.DataFrame(index=index,columns=columns)
        
        # Masks for computing efficiency
        mask1=np.logical_not(np.isnan(arrcoveren))
        mask2=self.__arrbio>0.0
        mask12=mask1*mask2
        
        for c in coverreg:
            
            cont=c-1
            
            biomassdata=self.__arrbio*mask12
            biomassdata[biomassdata[:,:,:]<=0.0]=np.nan
            coverdata=arrcoveren*mask12
            
            print('Fitting model for land cover '+covernames[cont])
            
            mask=coverdata[:,:,:]==c            
            
            biomassdataar=biomassdata[mask]
            logbiomassdata=np.log10(biomassdataar)
            
            expdist[:,cont]=np.percentile(logbiomassdata,p)
            print('Log Average='+str(np.mean(logbiomassdata)))
            print('Average='+str(np.power(10.0,np.mean(logbiomassdata))))
            print('St. Dev.='+str(np.std(logbiomassdata)))
            
            dfcpf.loc[:,columns[cont]]=np.power(10.0,expdist[:,cont])
            
        dfcpf.loc[:,'Prob.']=p
        
        rtcpf=self.__workspace+'data'+delim+'fdas'+delim+'Biomass_cpf_unc.csv'
        dfcpf.to_csv(rtcpf,sep=',')
        
        # del arrcoveren,mask,biomassdata,coverdata,logbiomassdata,biomassdataar
        # del mask1,mask2,mask12,arrcover
        
        # gc.collect()
        
        # return expdist
        
        return self.__arrbio,arrcoveren,mask,biomassdata,coverdata,expdist
        
        
    #End of method ------------------------------------------------------------