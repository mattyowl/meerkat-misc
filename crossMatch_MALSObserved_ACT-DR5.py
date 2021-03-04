"""

Cross match Neeraj's MALS .csv files with the ACT DR5 calalog

"""

import os, sys
import numpy as np
import astropy.table as atpy
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
import IPython

#------------------------------------------------------------------------------------------------------------
def crossMatch(refCatalog, matchCatalog, radiusArcmin = 2.5):
    """Cross matches matchCatalog onto refCatalog for objects found within some angular radius 
    (specified in arcmin).
    
    Args:
        refCatalog (:obj:`astropy.table.Table`): The reference catalog.
        matchCatalog (:obj:`astropy.table.Table`): The catalog to match onto the reference catalog.
        radiusArcmin (float, optional): Cross-match radius in arcmin.
    
    Returns:
        Cross-matched reference catalog, matchCatalog, and array of angular separation in degrees, for 
        objects in common within the matching radius. The cross matched columns are sorted such that rows in
        each correspond to the matched objects.
    
    """
    
    inTab=refCatalog
    outTab=matchCatalog
    RAKey1, decKey1=getTableRADecKeys(inTab)
    RAKey2, decKey2=getTableRADecKeys(outTab)
    cat1=SkyCoord(ra = inTab[RAKey1].data, dec = inTab[decKey1].data, unit = 'deg')
    xMatchRadiusDeg=radiusArcmin/60.
    cat2=SkyCoord(ra = outTab[RAKey2].data, dec = outTab[decKey2].data, unit = 'deg')
    xIndices, rDeg, sep3d = match_coordinates_sky(cat1, cat2, nthneighbor = 1)
    mask=np.less(rDeg.value, xMatchRadiusDeg)  
    matched_outTab=outTab[xIndices]
        
    inTab=inTab[mask]
    matched_outTab=matched_outTab[mask]
    rDeg=rDeg.value[mask]
    
    return inTab, matched_outTab, rDeg

#------------------------------------------------------------------------------------------------------------
def removeCrossMatched(refCatalog, matchCatalog, radiusArcmin = 2.5):
    """Cross matches matchCatalog onto refCatalog for objects found within some angular radius 
    (specified in arcmin), and returns refCatalog with the matching entries removed.
    
    Args:
        refCatalog (:obj:`astropy.table.Table`): The reference catalog.
        matchCatalog (:obj:`astropy.table.Table`): The catalog to match onto the reference catalog.
        radiusArcmin (float, optional): Cross-match radius in arcmin.
    
    Returns:
        Cross-matched reference catalog (:obj:`astropy.table.Table`) with matches to matchCatalog removed.
        
    """
        
    inTab=refCatalog
    outTab=matchCatalog
    RAKey1, decKey1=getTableRADecKeys(inTab)
    RAKey2, decKey2=getTableRADecKeys(outTab)
    cat1=SkyCoord(ra = inTab[RAKey1].data, dec = inTab[decKey1].data, unit = 'deg')
    xMatchRadiusDeg=radiusArcmin/60.
    cat2=SkyCoord(ra = outTab[RAKey2].data, dec = outTab[decKey2].data, unit = 'deg')
    xIndices, rDeg, sep3d = match_coordinates_sky(cat1, cat2, nthneighbor = 1)
    mask=np.greater(rDeg.value, xMatchRadiusDeg)  
    inTab=inTab[mask]
    
    return inTab
    
#------------------------------------------------------------------------------------------------------------
def getTableRADecKeys(tab):
    """Returns the column names in the table in which RA, dec coords are stored, after trying a couple of 
    variations.
    
    Args:
        tab (:obj:`astropy.table.Table`): The table to search.
        
    Returns:
        Name of RA column, name of dec. column
    
    """
    RAKeysToTry=['ra', 'RA', 'RADeg']
    decKeysToTry=['dec', 'DEC', 'decDeg', 'Dec']
    RAKey, decKey=None, None
    for key in RAKeysToTry:
        if key in tab.keys():
            RAKey=key
            break
    for key in decKeysToTry:
        if key in tab.keys():
            decKey=key
            break
    if RAKey is None or decKey is None:
        raise Exception("Couldn't identify RA, dec columns in the supplied table.")
    
    return RAKey, decKey

#------------------------------------------------------------------------------------------------------------
# Main

xMatchRadiusArcmin=30

# These will be grafted from MALS onto the ACT DR5 catalog
wantedCols=['source_id', 'ra', 'dec', 'flux', 'z_unique']
    
act=atpy.Table().read("DR5_cluster-catalog_v1.1.fits")

actXmals=[]
fileNames=["MALS-NVSS-observed-2021Feb26.csv", "MALS-SUMSS-observed-2021Feb26.csv"]
for f in fileNames:
    mals=atpy.Table().read(f)
    # Don't know why some MALS objects have 'null' where coordinates should be
    mals=mals[mals['ra'] != 'null']
    mals['ra']=np.array(mals['ra'], dtype = float)
    mals['dec']=np.array(mals['dec'], dtype = float)
    xMatch_act, xMatch_mals, rDeg=crossMatch(act, mals, radiusArcmin = xMatchRadiusArcmin)
    for key in wantedCols:
        if key not in xMatch_act:
            xMatch_act[key]=xMatch_mals[key]
    actXmals.append(xMatch_act)

actXmals=atpy.vstack(actXmals)
actXmals.write("ACTDR5-x-MALS-2021Feb26.fits", overwrite = True)

