# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 11:42:24 2023

@author: oscar

This class runs the sequential gaussian simulation algorithm implemented 
in geostatspy

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
import pandas as pd
import geostatspy.geostats as geostats
import geostatspy.GSLIB as GSLIB

class rdSGSIM:
    
    #Begining of constructor --------------------------------------------------
    def __init__(self,project,parameters,bioarray,verbose=False):
        """
        Construction of an instance for reading K-AGB basic execution parameters

        Parameters:
        ----------
            project : string
                The folder location of the project
            parameters : dictionary
                Basic parameters of the project
            bioarray: numpy float 32 array
                Biomass secondary information
        
        Returns:
        -------
            null
            
        """
        
        #Constructor begins here
        self.__workspace=project
        self.__sgsimparams={}
        self.__dfbiodata=0.0
        
        # Read the biomass data and transform it into a DataFrame
        
        rows=parameters['rows']
        columns=parameters['columns']
        r=parameters['res']
        minx=parameters['minX']
        miny=parameters['minY']
        
        self.__dfbiodata=pd.DataFrame(columns=['East','North','Value'])
        
        x = np.linspace(minx+r/2.0,minx+r*columns,num=columns)
        y = np.linspace(miny+r*rows-r/2.0,miny,num=rows)
        xx, yy = np.meshgrid(x, y)
        xar=np.reshape(xx,newshape=(rows*columns),order='C')
        yar=np.reshape(yy,newshape=(rows*columns),order='C')
        bar=np.reshape(bioarray,newshape=(rows*columns),order='C')
        mask=~np.isnan(bar)
        xar=xar[mask]
        yar=yar[mask]
        bar=bar[mask]
        
        self.__dfbiodata['East'] = xar.tolist()
        self.__dfbiodata['North'] = yar.tolist()
        self.__dfbiodata['Value'] = bar.tolist()
        
        self.__dfbiodata['Nbio'], tvBio, tnsBio = geostats.nscore(self.__dfbiodata,'Value')

        del x,y,xx,yy,xar,yar,bar,mask
        gc.collect()
        
        # Read sgsim parameters from file
        # SGSIM parameter are:
        '''
            df: Dataframe stroring conditioning dat
            xcol: Column in df storing x coordinates
            ycol: Column in df storing y coordinates
            vcol: Column in df storing conditioning values
            wcol: Column in df storing declustering weight (-1 means no weight)
            scol: Column in df storing second conditioning values (-1 means no secondary variable)
            tmin: Minimum allowable values for conditioning variable
            tmax: Maximum allowable values for conditioning variable
            itrans: if set to 0 then no transformation will be performed; 
                    the variable is assumed already standard normal (the simulation 
                    results will also be left unchanged). If itrans=1, transformations are performed.
            ismooth: if set to 0, then the data histogram, possibly with declustering 
                    weights is used for transformation, if set to 1, then the data are transformed 
                    according to the values in another file (perhaps from histogram smoothing)
            dftrans:
            tcol, twtcol: columns for the variable and the declustering weight 
            zmin, zmax: the minimum and maximum allowable data values. These are 
                        used in the back transformation procedure.
            ltail,ltpar: specify the back transformation implementation in the 
                         lower tail of the distribution: ltail=1 implements linear 
                         interpolation to the lower limit zmin, and ltail=2 implements 
                         power model interpolation, with w=ltpar, to the lower limit zmin.
                         The middle class interpolation is linear.
            utail, utpar: specify the back transformation implementation in the 
                          upper tail of the distribution: utail=1 implements linear 
                          interpolation to the upper limit zmax, utail=2 implements 
                          power model interpolation, with w=utpar, to the upper 
                          limit zmax, and utail=4 implements hyperbolic model 
                          extrapolation with w=utpar. The hyperbolic tail 
                          extrapolation is limited by zmax.
            nsim: the number of simulations to generate            
            nx, xmn, xsiz: definition of the grid system (x axis).
            ny, ymn, ysiz: definition of the grid system (y axis).
            seed: random number seed (a large odd integer)
            ndmin, ndmax: the minimum and maximum number of original data that 
                          should be used to simulate a grid node. If there are 
                          fewer than ndmin data points the node is not simulated.
            nodmax:
            mults: a multiple grid simulation will be performed if this is set 
                   to 1 (otherwise a standard spiral search for previously 
                         simulated nodes is considered)
            nmult: he number of multiple grid refinements to consider 
                    (used only if multgrid is set to 1).
            noct: the number of original data to use per octant. If this parameter 
                  is set less than or equal to 0, then it is not used; otherwise, 
                  it overrides the ndmax parameter and the data is partitioned 
                  into octants and the closest noct data in each octant is 
                  retained for the simulation of a grid node.
            radius: Major search radio
            radius1: Minor search radio
            sang1: the angle parameters that describe the orientation of the search ellipsoid
            mxctx=10,
            mxcty=10,
            ktype: the kriging type (0 = simple kriging, 1 = ordinary kriging, 
                                     2 = simple kriging with a locally varying mean, 
                                     3 = kriging with an external drift, or 
                                     4 = collocated cokriging with one secondary variable) 
                   used throughout the loop over all nodes. SK is required by theory; 
                   only in cases where the number of original data found in the 
                   neighborhood is large enough can OK be used without the 
                   risk of spreading data values beyond their range of influence
            colocorr: correlation coefficient to use for collocated cokriging (used only if ktype = 4)
            sec_map: Secondary variable map for colocated cokriging
            vario: Variogram model to perform the simulations
            
        '''
        
        #code begins here
        
        flnm=self.__workspace+'data'+delim+'sgsim'+delim+'config_sgsim.dat'
        
        try:
            
            with open(flnm, 'rb') as in_file:
                # A dictionary to contain info about the sgsim parameters                
                if verbose :
                    print(('\n    Reading file: "{0}"'.format(self.__flnm)))
                    
                print('*--- Preparing simulation, reading SGSIM parameters ---*')
                
                # Location of data in DataFrame
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['wcol']=int(line_split[1])
                print('wcol='+str(self.__sgsimparams['wcol'])) 
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['scol']=int(line_split[1])
                print('scol='+str(self.__sgsimparams['scol'])) 
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['tmin']=float(line_split[1])
                print('tmin='+str(self.__sgsimparams['tmin'])) 
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['tmax']=float(line_split[1])
                print('tmax='+str(self.__sgsimparams['tmax'])) 
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['itrans']=int(line_split[1])
                print('itrans='+str(self.__sgsimparams['itrans'])) 
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['ismooth']=int(line_split[1])
                print('ismooth='+str(self.__sgsimparams['ismooth']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['dftrans']=int(line_split[1])
                print('dftrans='+str(self.__sgsimparams['dftrans']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['tcol']=int(line_split[1])
                print('tcol='+str(self.__sgsimparams['tcol']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['twtcol']=int(line_split[1])
                print('twtcol='+str(self.__sgsimparams['twtcol']))
                
                # NS back transform parameters
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['zmin']=float(line_split[1])
                print('zmin='+str(self.__sgsimparams['zmin']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['zmax']=float(line_split[1])
                print('zmax='+str(self.__sgsimparams['zmax']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['ltail']=int(line_split[1])
                print('ltail='+str(self.__sgsimparams['ltail']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['ltpar']=float(line_split[1])
                print('ltpar='+str(self.__sgsimparams['ltpar']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['utail']=int(line_split[1])
                print('utail='+str(self.__sgsimparams['utail']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['utpar']=float(line_split[1])
                print('utpar='+str(self.__sgsimparams['utpar']))
                
                # Simulation parameters
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['nsim']=int(line_split[1])
                print('nsim='+str(self.__sgsimparams['nsim']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['seed']=int(line_split[1])
                print('seed='+str(self.__sgsimparams['seed']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['ndmin']=int(line_split[1])
                print('ndmin='+str(self.__sgsimparams['ndmin']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['ndmax']=int(line_split[1])
                print('ndmax='+str(self.__sgsimparams['ndmax']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['nodmax']=int(line_split[1])
                print('nodmax='+str(self.__sgsimparams['nodmax']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['mults']=int(line_split[1])
                print('mults='+str(self.__sgsimparams['mults']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['nmult']=int(line_split[1])
                print('nmult='+str(self.__sgsimparams['nmult']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['noct']=int(line_split[1])
                print('noct='+str(self.__sgsimparams['noct']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['radius']=float(line_split[1])
                print('radius='+str(self.__sgsimparams['radius']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['radius1']=float(line_split[1])
                print('radius1='+str(self.__sgsimparams['radius1']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['sang1']=int(line_split[1])
                print('sang1='+str(self.__sgsimparams['sang1']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['mxctx']=int(line_split[1])
                print('mxctx='+str(self.__sgsimparams['mxctx']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['mxcty']=int(line_split[1])
                print('mxcty='+str(self.__sgsimparams['mxcty']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['ktype']=int(line_split[1])
                print('ktype='+str(self.__sgsimparams['ktype']))
                
                line=in_file.readline()                
                line_split=line.split()
                self.__sgsimparams['colocorr']=float(line_split[1])
                print('colocorr='+str(self.__sgsimparams['colocorr']))
                
                self.__sgsimparams['sec_map']=0
                
        except IOError:
            print(('    Error reading file "{0}"'.format(flnm)))
            print('    Check if the file exists...')
            sys.exit(0)
        
        in_file.close()

        print('*---Succesfull reading of the sgsim parameters in file '+flnm+'---*')
        
    #End of constructor -------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def readSgsimVariogram(self):
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
        
        nug=0
        nst=0
        it=[]
        variop=[]        
        vario={}
        
        try:
            
            with open(flnm, 'rb') as in_file:                
                    
                print('*--- Preparing simulation, reading SGSIM variogram ---*')
                
                # Location of data in DataFrame
                
                line=in_file.readline()  
                
                line=in_file.readline()                
                line_split=line.split()
                nug=float(line_split[1])
                vario['nug']=nug
                print('nug='+str(vario['nug'])) 
                
                line=in_file.readline()                
                line_split=line.split()
                nst=int(line_split[1])
                vario['nst']=nst
                print('nst='+str(vario['nst']))
                
                for k in range(nst):
                    
                    line=in_file.readline()                
                    line_split=line.split()
                    it.append(int(line_split[1]))
                    print('it'+str(k+1)+'='+str(it[k]))
                    
                    line=in_file.readline()                
                    line_split=line.split()
                    cc=float(line_split[1])
                    variop.append(cc)
                    print('cc'+str(k+1)+'='+str(cc))
                    
                    line=in_file.readline()                
                    line_split=line.split()
                    azi=float(line_split[1])
                    variop.append(azi)
                    print('azi'+str(k+1)+'='+str(azi))
                    
                    line=in_file.readline()                
                    line_split=line.split()
                    hmaj=float(line_split[1])
                    variop.append(hmaj)
                    print('hmaj'+str(k+1)+'='+str(hmaj))
                    
                    line=in_file.readline()                
                    line_split=line.split()
                    hmin=float(line_split[1])
                    variop.append(hmin)
                    print('hmin'+str(k+1)+'='+str(hmin))
                    
                vario['it']=it
                vario['params']=variop
                
                
        except IOError:
            print(('    Error reading file "{0}"'.format(flnm)))
            print('    Check if the file exists...')
            sys.exit(0)
        
        in_file.close()

        print('*---Succesfull reading of the execution data in file '+flnm+'---*')
        
        return vario
    
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def getBiomassdf(self):
        """
        Returns the biomass dataframe

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            dfbiodata : Data Frame
                Biomass data from secondary information
        """
            
        #code begins here
        return self.__dfbiodata
    
    #End of method ------------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def getSgsimparams(self):
        """
        Returns the sgsim algoritm parameters dictionary

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            sgsimparams : Dictionary
                sgsim algoritm parameters dictionary, key are explained in constructor
        """
            
        #code begins here
        return self.__sgsimparams
    
    #End of method ------------------------------------------------------------
    
    