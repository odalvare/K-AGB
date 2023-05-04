# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 14:09:58 2023

Auxiliar functions in the downscaling process

@author: oscar
"""

import rasterio
from rasterio.enums import Resampling

import numpy as np

'''
###############################################################################
# FUNCTIONS IMPLEMENTING THE RESAMPLING PROCEDURE VIA RASTERIO
###############################################################################
'''

#Begining of function ---------------------------------------------------------
def resample_ras(r_src,r_trg,rtoutput,extent,crs):
    '''
     This functions  resamples the r_src map to another map having the 
     resolution of s_trg raster. Both input maps are represented using 
     np.float32 data type

    Parameters
    ----------
    r_src : rasterio raster (FLOAT 32)
        Source raster, it will be resampled to r_trg resolution
    r_trg : rasterio raster (FLOAT 32)
        Target raster. The resampling procedure uses its metadata for resampling
    rtoutput : str
        Destination file for writing resampled r_src
    extent : list
        Geographic coordinates for the resampled map
    crs : str
        Geographic system of the output resampled raster

    Returns
    -------
    None. The rasampled raster is saved in disk

    '''
    
    # The function starts here
    
    nxsrc=r_src.width                     # Number of columns
    nysrc=r_src.height                    # Number of rows
    rxsrc,ryndvi=r_src.res
    bndssrc=r_src.bounds
    
    nxtrg=r_trg.width                     # Number of columns
    nytrg=r_trg.height                    # Number of rows
    rxtrg,ryndvi=r_trg.res
    bndstrg=r_trg.bounds
    
    upscale_factorh=nytrg/nysrc    #target/source
    upscale_factorw=nxtrg/nxsrc
    r_srcres = r_src.read(
        out_shape=(
            r_src.count,
            int(r_src.height*upscale_factorh),
            int(r_src.width*upscale_factorw)
        ),
        resampling=Resampling.cubic_spline
        # resampling=Resampling.bilinear
    )    
    transform = r_src.transform * r_src.transform.scale(
        (r_src.width / r_srcres.shape[-1]),
        (r_src.height / r_srcres.shape[-2])
    )
    
    # Read each layer and write it to stack
    meta=r_trg.meta    
    meta.update(count = 1)
    with rasterio.open(rtoutput, 'w', **meta) as dst:
        dst.write(r_srcres[0,:,:],1)
    dst.close()
    
#End of function --------------------------------------------------------------
    
def resample_ras_int(r_src,r_trg,rtoutput,extent,crs):
    '''
     This functions  resamples the r_src map to another map having the 
     resolution of s_trg raster. Both input maps are represented using 
     np.int64 data type

    Parameters
    ----------
    r_src : rasterio raster (INT64)
        Source raster, it will be resampled to r_trg resolution
    r_trg : rasterio raster (INT64)
        Target raster. The resampling procedure uses its metadata for resampling
    rtoutput : str
        Destination file for writing resampled r_src
    extent : list
        Geographic coordinates for the resampled map
    crs : str
        Geographic system of the output resampled raster

    Returns
    -------
    None. The rasampled raster is saved in disk

    '''
    
    nxsrc=r_src.width                     # Number of columns
    nysrc=r_src.height                    # Number of rows
    rxsrc,ryndvi=r_src.res
    bndssrc=r_src.bounds
    
    nxtrg=r_trg.width                     # Number of columns
    nytrg=r_trg.height                    # Number of rows
    rxtrg,ryndvi=r_trg.res
    bndstrg=r_trg.bounds
    
    upscale_factorh=nytrg/nysrc           # target/source
    upscale_factorw=nxtrg/nxsrc
    r_srcres = r_src.read(
        out_shape=(
            r_src.count,
            int(r_src.height*upscale_factorh),
            int(r_src.width*upscale_factorw)
        ),
        resampling=Resampling.nearest
    )    
    transform = r_src.transform * r_src.transform.scale(
        (r_src.width / r_srcres.shape[-1]),
        (r_src.height / r_srcres.shape[-2])
    )
    
    # Read each layer and write it to stack
    meta=r_trg.meta    
    meta.update(count = 1)
    meta.update(dtype = 'int32')
    with rasterio.open(rtoutput, 'w', **meta) as dst:
        dst.write(r_srcres[0,:,:],1)
    dst.close()
    
#End of function --------------------------------------------------------------
