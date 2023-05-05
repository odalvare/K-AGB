# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 14:09:14 2023

@author: oscar

Performs the biomass downscaling procedure using Linear, Exponential or
Power regression

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

from sklearn.linear_model import LinearRegression

from downscaling.downscaling_help import resample_ras

class DownScaling:
    
    #Begining of constructor --------------------------------------------------
    def __init__(self,project,params):
        """
        Construction of an instance for performing linear multivariate bias
        corrected downscaling of biomass

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
        
        # Open simulated biomass
        rt=self.__workspace+'data'+delim+'downscaling'+delim+'sgsim_biosim.tif'
        self.__biomassrea=rasterio.open(rt)
        self.__metaSRC=self.__biomassrea.meta
        
        nxbio=self.__biomassrea.width                     # Number of columns
        nybio=self.__biomassrea.height                    # Number of rows
        rxbio,rybio=self.__biomassrea.res
        
        biomass_bnds = self.__biomassrea.bounds
        minxbio=biomass_bnds[0]
        minybio=biomass_bnds[1]
        maxxbio=biomass_bnds[2]
        maxybio=biomass_bnds[3]    
        
        self.__extentbio = [minxbio,maxxbio,minybio,maxybio]
        
        nsim=self.__params['sgsim_reas']
        self.__arrbiosim=np.zeros((nsim,nybio,nxbio),dtype=np.float32)
        for i in range(nsim):
            self.__arrbiosim[i,:,:]=self.__biomassrea.read(i+1)
        self.__arrbiosim[self.__arrbiosim[:,:,:]>150.0]=150.0
        
    #End of constructor -------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def run(self):
        """
        Runs the downscaling process

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            arrbiodown : numpy float array (nsim,nrow,ncols)
                Downscaled biomass realizations (from SGSIM)
            
        """
        
        # ******************* RESAMPLIG TO COARSE RESOUTION *******************
        
        # Open NDVI (250 m) target resolution ---------------------------------
        
        rt=self.__workspace+'data'+delim+'downscaling'+delim+'ndvi_original.tif'
        ndvior = rasterio.open(rt)
        
        nxndvior=ndvior.width                     # Number of columns
        nyndvior=ndvior.height                    # Number of rows
        rxndvior,ryndvior=ndvior.res
        
        ndvi_bndsor = ndvior.bounds
        minxndvior=ndvi_bndsor[0]
        minyndvior=ndvi_bndsor[1]
        maxxndvior=ndvi_bndsor[2]
        maxyndvior=ndvi_bndsor[3]
        extentndvior = [minxndvior,maxxndvior,minyndvior,maxyndvior]
        
        arrndvior=ndvior.read(1)
        arrndvior=arrndvior.astype(np.float32)
        arrndvior[arrndvior[:]<=0]=np.nan
        
        # Resample NDVI to 1 km spatial resolution (biomass resolution)
        
        rtres=self.__workspace+'data'+delim+'downscaling'+delim+'ndvi_res.tif'
        resample_ras(ndvior,                  # origin ndvi (fine) - target biomass (coarse)
                     self.__biomassrea,
                     rtres,
                     self.__extentbio,
                     self.__params['EPSG'])                                    
        
        ndvi = rasterio.open(rtres)
        arrndvi=ndvi.read(1)
        arrndvi=arrndvi.astype(np.float32)
        arrndvi[arrndvi[:,:]<=0]=np.nan
        
        # Open EVI (250 m) target resolution ---------------------------------
        
        rt=self.__workspace+'data'+delim+'downscaling'+delim+'evi_original.tif'
        evior = rasterio.open(rt)
        
        arrevior=evior.read(1)
        arrevior=arrevior.astype(np.float32)
        arrevior[arrevior[:]<=0]=np.nan
        
        # Resample EVI to 1 km spatial resolution (biomass resolution)
        
        rtres=self.__workspace+'data'+delim+'downscaling'+delim+'evi_res.tif'
        resample_ras(evior,                  # origin evi (fine) - target biomass (coarse)
                     self.__biomassrea,
                     rtres,
                     self.__extentbio,
                     self.__params['EPSG'])                                    
        evior.close()
        
        evi = rasterio.open(rtres)
        arrevi=evi.read(1)
        arrevi=arrevi.astype(np.float32)
        arrevi[arrevi[:,:]<=0]=np.nan
        
        # Open LAI (250 m) target resolution ----------------------------------
        
        rt=self.__workspace+'data'+delim+'downscaling'+delim+'lai_original.tif'
        laior = rasterio.open(rt)
        
        arrlaior=laior.read(1)
        arrlaior=arrlaior.astype(np.float32)
        arrlaior[arrlaior[:]<=0]=np.nan
        
        # Resample LAI to 1 km spatial resolution (biomass resolution) -------
        
        rtres=self.__workspace+'data'+delim+'downscaling'+delim+'lai_res.tif'
        resample_ras(laior,                  # origin evi (fine) - target biomass (coarse)
                      self.__biomassrea,
                      rtres,
                      self.__extentbio,
                      self.__params['EPSG'])                                   
        laior.close()
        
        lai = rasterio.open(rtres)
        arrlai=lai.read(1)
        arrlai=arrlai.astype(np.float32)
        arrlai[arrlai[:,:]<=0]=np.nan
        
        # ******************* RESAMPLIG TO COARSE RESOUTION *******************
        
        ndvichs=np.reshape(arrndvi,newshape=self.__params['rows']*self.__params['columns'],order='C')
        evichs=np.reshape(arrevi,newshape=self.__params['rows']*self.__params['columns'],order='C')
        laichs=np.reshape(arrlai,newshape=self.__params['rows']*self.__params['columns'],order='C')
        
        # Transform data according to the type of regression
        
        nvbles=3
        x=np.zeros((ndvichs.shape[0],nvbles),dtype=np.float32)
        y=np.zeros((ndvichs.shape[0]),dtype=np.float32)
        
        if(self.__params['down_id']==1):            
            print('Linear regression, no data transform')
            x[:,0]=ndvichs[:]
            x[:,1]=evichs[:]
            x[:,2]=laichs[:]
            
        elif(self.__params['down_id']==2):            
            print('Exponential regression, logarithmic transform of biomass')
            logbio=np.log10(self.__arrbiosim)
            x[:,0]=ndvichs[:]
            x[:,1]=evichs[:]
            x[:,2]=laichs[:]
            
        elif(self.__params['down_id']==3):            
            print('Power regression, logarithmic transform of biomass, ndvi, evi, and lai')
            logbio=np.log10(self.__arrbiosim)
            x[:,0]=np.log10(ndvichs)
            x[:,1]=np.log10(evichs)
            x[:,2]=np.log10(laichs)
            
        else:
            raise Exception('Invalid downscaling scheme')
            sys.exit()
        
        arrbiodown=np.zeros((self.__params['sgsim_reas'],nyndvior,nxndvior),dtype=np.float64)
        
        for rea in range(self.__params['sgsim_reas']):
            
            print('Downscaling realization '+str(rea+1)+' **********************************')
            
            if(self.__params['down_id']==1):                
                y[:]=np.reshape(self.__arrbiosim[rea,:,:],
                                newshape=self.__params['rows']*self.__params['columns'],
                                order='C')
                
            elif(self.__params['down_id']==2):                
                y[:]=np.reshape(logbio[rea,:,:],
                                newshape=self.__params['rows']*self.__params['columns'],
                                order='C')
            
            elif(self.__params['down_id']==3):                
                y[:]=np.reshape(logbio[rea,:,:],
                                newshape=self.__params['rows']*self.__params['columns'],
                                order='C')
                
            else:                
                raise Exception('Invalid downscaling scheme')
                sys.exit()
            
            # ******************* REGRESSION *********************************    
            
            model = LinearRegression()
            model.fit(x,y)
            r_sq = model.score(x,y)
            print(f"intercept: {model.intercept_}")
            print(f"slope: {model.coef_}")
            
            # Linear regression
            intercept=model.intercept_
            slopes=model.coef_
            
            if(self.__params['down_id']==2):
                biout=intercept+slopes[0]*arrndvi+slopes[1]*arrevi+slopes[2]*arrlai
                biop=np.power(10.0,biout)
            elif(self.__params['down_id']==3):
                biout=intercept+slopes[0]*np.log10(arrndvi)+slopes[1]*np.log10(arrevi)+slopes[2]*np.log10(arrlai)
                biop=np.power(10.0,biout)
            else:
                biout=intercept+slopes[0]*arrndvi+slopes[1]*arrevi+slopes[2]*arrlai
                biop=biout
            
            # ******************* BIAS ESTIMATION ****************************    
            
            resarr=self.__arrbiosim[rea,:,:]-biop[:,:]
            
            rtres=self.__workspace+'data'+delim+'downscaling'+delim+'errb_original.tif'
            meta=self.__metaSRC.copy()
            meta.update(count = 1)
            with rasterio.open(rtres, 'w', **meta) as dst:
                dst.write(resarr[:,:],1)
            dst.close()
            
            # Open residual and resample to original ndvi, evi, lai resolution
            
            residual=rasterio.open(rtres)
            resarr=residual.read(1)
            resarr=resarr.astype(np.float32)
                
            rtresidualres=self.__workspace+'data'+delim+'downscaling'+delim+'errb_res.tif'
            resample_ras(residual,ndvior,rtresidualres,extentndvior,'EPSG:3116')       #target/source    
            residual.close()
            
            residualres = rasterio.open(rtresidualres)
            resarr=residualres.read(1)
            resarr=resarr.astype(np.float32)
            residualres.close()
            
            # ****************REGRESSION + BIAS CORRECTION *******************
            
            # Apply regression models to the data at the fine resolution and add error
            
            if(self.__params['down_id']==2):
                bioast=intercept+slopes[0]*arrndvior+slopes[1]*arrevior+slopes[2]*arrlaior     #Exponential Regresssion
                bioast=np.power(10.0,bioast)+resarr
            elif(self.__params['down_id']==3):
                bioast=intercept+slopes[0]*np.log10(arrndvior)+slopes[1]*np.log10(arrevior)+slopes[2]*np.log10(arrlaior) #Power Regresssion
                bioast=np.power(10.0,bioast)+resarr
            else:
                bioast=intercept+slopes[0]*arrndvior+slopes[1]*arrevior+slopes[2]*arrlaior+resarr     #Linear Regresssion
            
            arrbiodown[rea,:,:]=bioast[:,:]
        
        ndvior.close()
        
        # Save Biomass
        
        rt=self.__workspace+'data'+delim+'fdas'+delim+'biosim_down.tif'
            
        # Write to TIFF
        kwargs = ndvior.meta
        kwargs.update(
            dtype=rasterio.float64,
            count=arrbiodown.shape[0],
            compress='lzw')    
        with rasterio.open(rt, 'w', **kwargs) as dst:
            dst.write(arrbiodown)
            
        del x,y,bioast,biop,biout,logbio,resarr,arrndvior,arrevior,arrlaior
        del ndvichs,evichs,laichs,arrndvi,arrevi,arrlai
        
        gc.collect()
        
        return arrbiodown
        
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def getBiomassSRC(self):
        """
        Returns the biomass rasterio object simulated in sgsim

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            biomassrea : rasterio multiband raster
                raster storing the biomass simulated via sgsim
        """
            
        #code begins here
        return self.__biomassrea
    
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def getBiomassArray(self):
        """
        Returns the biomass array without performing downscaling

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            arrbiosim : numpy array (float 32)
                array of simulated biomass realization
        """
            
        #code begins here
        return self.__arrbiosim
    
    #End of method ------------------------------------------------------------