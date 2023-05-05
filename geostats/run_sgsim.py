# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 10:55:12 2023

@author: oscar
"""

import platform
import gc

import rasterio

#Detect operaring system and apply folder delimiter
ost=platform.system()
delim='/'
if(ost=='Windows'):
    delim='//'

import numpy as np
import geostatspy.geostats as geostats
import geostatspy.GSLIB as GSLIB

class runSGSIM:
    
    #Begining of constructor --------------------------------------------------
    def __init__(self,project,dfbio,bparams,bmeta,variodict,sgsimparams):
        """
        Construction of an instance for running sgsim via geostatspy

        Parameters:
        ----------
            project : string
                The folder location of the project
            parameters : dictionary
                Basic parameters of the project
            dfbio : pandas DataFrame
                Biomass secondary information stored in a DataFrame
            variodict : dictionary
                Dictionary with the variogram
            sgsimparams : dictionary
                Dictionary wit the paraneters of geostatspy sgsim routine
        
        Returns:
        -------
            null
            
        """
        
        #Constructor begins here
        self.__workspace=project
        self.__dfbio=dfbio
        self.__params=bparams
        self.__dfbio=dfbio
        self.__variodict=variodict
        self.__sgsimparams=sgsimparams
        self.__bmeta=bmeta
        
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
        
    #End of constructor -------------------------------------------------------
    
    #Begining of method -------------------------------------------------------
    def runSimulation(self):
        """
        Runs SGSIM using the current settings

        Parameters:
        ----------        
            null
        
        Returns:
        -------
            sims : numpy float array
                Simulated biomass array
                
        """
            
        #code begins here
        
        rows=self.__params['rows']
        columns=self.__params['columns']
        r=self.__params['res']
        minx=self.__params['minX']
        miny=self.__params['minY']
        
        it=self.__variodict['it']
        variop=self.__variodict['params']
        
        if(self.__variodict['nst']==1):
            print('Variogram with nested structure, 1 anisotropic models')
            vario = GSLIB.make_variogram(nug=self.__variodict['nug'],
                                         nst=self.__variodict['nst'],
                                         it1=it[0],cc1=variop[0],azi1=variop[1],hmaj1=variop[2],hmin1=variop[3])
            print(vario)
        elif(self.__variodict['nst']==2):
            print('Variogram with nested structure, 2 anisotropic models')
            vario = GSLIB.make_variogram(nug=self.__variodict['nug'],
                                         nst=self.__variodict['nst'],
                                         it1=it[0],cc1=variop[0],azi1=variop[1],hmaj1=variop[2],hmin1=variop[3],
                                         it2=it[1],cc2=variop[4],azi2=variop[5],hmaj2=variop[6],hmin2=variop[7])
            print(vario)
        
        sims=np.zeros((self.__sgsimparams['nsim'],rows,columns),dtype=np.float32)
        
        for i in range(self.__sgsimparams['nsim']):
            
            seed=np.random.randint(100000)
            
            print('  ')
            print('***************************************************************')
            print('***************************************************************')
            print('Simulation '+str(i+1))
            print('Seed='+str(seed))
            
            sim_sk = geostats.sgsim(self.__dfbio,
                                    'East',
                                    'North',
                                    'Value',
                                    wcol=self.__sgsimparams['wcol'],
                                    scol=self.__sgsimparams['scol'],
                                    tmin=self.__sgsimparams['tmin'],
                                    tmax=self.__sgsimparams['tmax'],
                                    itrans=self.__sgsimparams['itrans'],
                                    ismooth=self.__sgsimparams['ismooth'],
                                    dftrans=self.__sgsimparams['dftrans'],
                                    tcol=self.__sgsimparams['tcol'],
                                    twtcol=self.__sgsimparams['twtcol'],
                                    zmin=self.__sgsimparams['zmin'],
                                    zmax=self.__sgsimparams['zmax'],
                                    ltail=self.__sgsimparams['ltail'],
                                    ltpar=self.__sgsimparams['ltpar'],
                                    utail=self.__sgsimparams['utail'],
                                    utpar=self.__sgsimparams['utpar'],
                                    nsim=1,
                                    nx=columns,
                                    xmn=minx,
                                    xsiz=r,
                                    ny=rows,
                                    ymn=miny,
                                    ysiz=r,
                                    seed=seed,
                                    ndmin=self.__sgsimparams['ndmin'],
                                    ndmax=self.__sgsimparams['ndmax'],
                                    nodmax=self.__sgsimparams['nodmax'],
                                    mults=self.__sgsimparams['mults'],
                                    nmult=self.__sgsimparams['nmult'],
                                    noct=self.__sgsimparams['noct'],
                                    radius=self.__sgsimparams['radius'],
                                    radius1=self.__sgsimparams['radius1'],
                                    sang1=self.__sgsimparams['sang1'],
                                    mxctx=self.__sgsimparams['mxctx'],
                                    mxcty=self.__sgsimparams['mxcty'],
                                    ktype=self.__sgsimparams['ktype'],
                                    colocorr=self.__sgsimparams['colocorr'],
                                    sec_map=self.__sgsimparams['sec_map'],
                                    vario=vario)
        
            sims[i,:,:]=sim_sk[:,:]
            del sim_sk
            
            gc.collect()
        
        # Write sgsim simulations
        
        # Update meta to reflect the number of layers
        self.__bmeta.update(count = self.__sgsimparams['nsim'])
        
        rt=self.__workspace+'data'+delim+'downscaling'+delim+'sgsim_biosim.tif'
        
        # Read each layer and write it to stack
        with rasterio.open(rt, 'w', **self.__bmeta) as dst:
            for i in range(self.__sgsimparams['nsim']):
                dst.write_band(i+1, sims[i,:,:])
        dst.close()
        
        return sims        
    
    #End of method ------------------------------------------------------------