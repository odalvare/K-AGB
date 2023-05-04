# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 08:48:23 2023

@author: oscar

This class reads the primary execution parameter for K-AGB.

"""

import sys
import os
import platform
import gc

#Detect operaring system and apply folder delimiter
ost=platform.system()
delim='/'
if(ost=='Windows'):
    delim='//'

import numpy as np
import rasterio

class bldKAGB:
    
    #Begining of constructor --------------------------------------------------
    def __init__(self,project):
        """
        Construction of an instance for reading K-AGB basic execution parameters

        Parameters:
        ----------
            workspace : string
                The folder location of the project
        
        Returns:
        -------
            null
            
        """
        #Constructor begins here
        self.__workspace=project
        self.__parameters={}
        self.__metaref={}
        
        self.readkagbdata()
        
        # Read secondary information and saves array and metadata
        
        rt=self.__workspace+delim+'data'+delim+'sgsim'+delim+'sbiomass.tif'
        biomass = rasterio.open(rt)
        self.__metaref=biomass.meta
        
        nxbio=biomass.width                     # Number of columns
        nybio=biomass.height                    # Number of rows
        rxbio,rybio=biomass.res
        r=rxbio
        
        biomass_bnds = biomass.bounds
        minxbio=biomass_bnds[0]
        minybio=biomass_bnds[1]
        
        kwargs=biomass.meta.copy()
        crs=kwargs['crs']
        code=''
        if crs.is_epsg_code:
            code = int(crs['init'].lstrip('epsg:'))
        self.__agbarray=biomass.read(1)
        self.__agbarray=self.__agbarray.astype(np.float32)
        self.__agbarray[self.__agbarray[:]<=0]=np.nan
        biomass.close()
        
        if(nxbio!=self.__parameters['columns']):
            raise Exception('Invalid number of columns in biomass data')
            sys.exit()
        if(nybio!=self.__parameters['rows']):
            raise Exception('Invalid number of rows in biomass data')
            sys.exit()
        if(minxbio-self.__parameters['minX']>1.0e-3):
            raise Exception('Invalid most western coordinate in biomass data')
            sys.exit()
        if(minybio-self.__parameters['minY']>1.0e-3):
            raise Exception('Invalid most southern coordinate in biomass data')
            sys.exit()
        if(r-self.__parameters['res']>1.0e-3):
            raise Exception('Invalid resolution in biomass data')
            sys.exit()
        epsg='EPSG:'+str(code)
        print(epsg)
        if(epsg!=self.__parameters['EPSG']):
            raise Exception('Invalid coordinate system in biomass data')
            sys.exit()
        
    #End of constructor -------------------------------------------------------
    
    #Begining of method -------------------------------------------------------    
    def readkagbdata(self,verbose=False): 
        """
        Method for reading the KAGB basic settings
        
        Parameters
        ----------        
            verbose: logical
                flag to detect error in reading input data
        
        Returns
        -------        
            parameters : dictionary
                required data for running NS EnKF.
                keys of parameters dictionary:
                    name: Project's name (str)
                    EPSG: Geographic reference of the project (str)
                    minX: Most Western coordinate in the working area (integer)
                    minY: Most Southern coordinate in the working area (integer)
                    res: Pixel size (float)
                    rows: Number of rows on which the domain is discretized (integer)
                    cols: Number of columns on which the domain is discretized (integer)
                    n_years: Number of years with land cover data (integer)
                    n_categories: Number of categories used in land cover (integer)
                    years: Years with land cover data (int array)
                    categories: Names of land cover categories (str list)
                    sgsim_reas: Number of realization to execute sgsim (int)
                    down_id: Type of multivariate bias corrected regression (int)
                        1: Linear regression
                        2: Exponential regression
                        3. Power regression
                    mc_reas: Number of realization to execute MC simulations (int)
                    n_vbles: Number of variables in multivariate bias corrected regression
                    variables: Name of the variables in multivariate bias corrected regression (str list)
                
        """
        
        #code begins here
        
        flnm=self.__workspace+'data'+delim+'basic'+delim+'config_kagb.dat'
        
        try:
            
            with open(flnm, 'rb') as in_file:
                # A dictionary to contain info about the grid                
                if verbose :
                    print(('\n    Reading file: "{0}"'.format(self.__flnm)))
                    
                print('*--- Preparing simulation, reading K-AGB information ---*')
                
                # Project name
                line=in_file.readline()
                line=in_file.readline()                
                line_split=line.split()
                self.__parameters['name']=line_split[0]
                print('Name='+str(self.__parameters['name'])) 
                
                # Geographic reference
                line=in_file.readline()
                line=in_file.readline()
                line=in_file.readline()
                line_split=line.split()
                c=str(line_split[0])
                self.__parameters['EPSG']=c[2:len(c)-1]
                print('EPSG='+str(self.__parameters['EPSG'])) 
                
                # Geographic reference
                line=in_file.readline()
                line=in_file.readline()
                line=in_file.readline()
                line_split=line.split()
                self.__parameters['minX']=float(line_split[1])
                print('minX='+str(self.__parameters['minX'])) 
                line=in_file.readline()
                line_split=line.split()
                self.__parameters['minY']=float(line_split[1])
                print('minY='+str(self.__parameters['minY'])) 
                line=in_file.readline()
                line_split=line.split()
                self.__parameters['res']=float(line_split[1])
                print('Res='+str(self.__parameters['res'])) 
                line=in_file.readline()
                line_split=line.split()
                self.__parameters['rows']=int(line_split[1])
                print('Rows='+str(self.__parameters['rows'])) 
                line=in_file.readline()
                line_split=line.split()
                self.__parameters['columns']=int(line_split[1])
                print('Columns='+str(self.__parameters['columns']))
                
                # Reading the simulation drivers
                line=in_file.readline()
                line=in_file.readline()
                line=in_file.readline()
                line_split=line.split()
                self.__parameters['n_years']=int(line_split[1])
                print('n_years='+str(self.__parameters['n_years']))
                line=in_file.readline()
                line_split=line.split()
                self.__parameters['n_categories']=int(line_split[1])
                print('n_categories='+str(self.__parameters['n_categories']))
                
                # Reading the years with landcover data
                line=in_file.readline()
                line=in_file.readline()
                years=np.zeros((self.__parameters['n_years']),dtype=np.int64)
                for i in range(self.__parameters['n_years']):
                    line=in_file.readline()
                    line_split=line.split()
                    years[i]=int(line_split[0])
                self.__parameters['years']=years
                
                # Reading the landcover categories
                line=in_file.readline()
                line=in_file.readline()
                categories=[]
                for i in range(self.__parameters['n_categories']):
                    line=in_file.readline()
                    line_split=line.split()
                    categories.append(line_split[0])
                self.__parameters['categories']=categories
                
                # Reading the SGSIM realizations
                line=in_file.readline()
                line=in_file.readline()
                line=in_file.readline()
                line_split=line.split()
                self.__parameters['sgsim_reas']=int(line_split[0])
                print('sgsim_reas='+str(self.__parameters['sgsim_reas']))
                
                # Reading the MonteCarlo realizations
                line=in_file.readline()
                line=in_file.readline()
                line=in_file.readline()
                line_split=line.split()
                self.__parameters['mc_reas']=int(line_split[0])
                print('mc_reas='+str(self.__parameters['mc_reas']))
                
                # Reading the AGB downscalint alternatives
                line=in_file.readline()
                line=in_file.readline()
                line=in_file.readline()
                line_split=line.split()
                self.__parameters['down_id']=int(line_split[0])
                print('down_id='+str(self.__parameters['down_id']))
                
                # Reading the number of variables for regression
                line=in_file.readline()
                line=in_file.readline()
                line=in_file.readline()
                line_split=line.split()
                self.__parameters['n_vbles']=int(line_split[0])
                print('n_vbles='+str(self.__parameters['n_vbles']))
                
                # Reading the downscaling variables
                line=in_file.readline()
                line=in_file.readline()
                vbles=[]
                for i in range(self.__parameters['n_vbles']):
                    line=in_file.readline()
                    line_split=line.split()
                    vbles.append(line_split[0])
                self.__parameters['variables']=vbles
                
        except IOError:
            print(('    Error reading file "{0}"'.format(flnm)))
            print('    Check if the file exists...')
            sys.exit(0)
            
        in_file.close()

        print('*---Succesfull reading of the execution data in file '+flnm+'---*')        
    
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def getSecAGB(self):
        '''
        
        This method return the above ground biomass secondary data and stores in RAM.
        Check compatibility with the project's metadata
        
        Parameters
        ----------        
            null        

        Returns
        -------
            agbarray: float array
                Above-ground biomass array

        '''
        
        return self.__agbarray
    
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def writeSgsimVariogram(self,variodict):
        """
        Returns the variogram model to execute sgsim

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            dfbiodata : Data Frame
                Biomass data from secondary information
        """
            
        #code begins here
        
        flnm=self.__workspace+'data'+delim+'sgsim'+delim+'vario.dat'
        params=variodict['params']
        
        try:
            
            with open(flnm, 'w') as f:
                
                line='#Variogram model for sgsim'  
                f.writelines(line)
                
                line='nug:		'+str(variodict['nug'])
                f.writelines('\n'+line)
                line='nst:		'+str(variodict['nst'])
                f.writelines('\n'+line)
                for m in range(variodict['nst']):
                    line='it'+str(m+1)+':		'+str(variodict['it'][m])  
                    f.writelines('\n'+line)
                    line='cc'+str(m+1)+':		'+str(params[4*m])
                    f.writelines('\n'+line)
                    line='azi'+str(m+1)+':		'+str(params[m*4+1])
                    f.writelines('\n'+line)
                    line='hmaj'+str(m+1)+':		'+str(params[m*4+2])
                    f.writelines('\n'+line)
                    line='hmin'+str(m+1)+':		'+str(params[m*4+3])
                    f.writelines('\n'+line)
                    
        except IOError:
            print(('    Error writing file "{0}"'.format(flnm)))
            print('    Check if the file exists...')
            sys.exit(0)
            
        f.close()
        
        del params
        gc.collect()
                
                
    #End of method ------------------------------------------------------------            
    
    #Begining of method -------------------------------------------------------
    def getParameters(self):
        """
        Returns the parameters dictionary

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            parameters : dictionary
                parameters of the kagb object
        """
            
        #code begins here
        return self.__parameters
    
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def getMetadateRef(self):
        """
        Returns the metadata dictionary

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            parameters : dictionary
                parameters of the biomass raster metadata
        """
            
        #code begins here
        return self.__metaref
    
    #End of method ------------------------------------------------------------
