#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This module generate a cloudfree Sentinel2 ImageCollection.
It is a adaptation of the javascript code:
http://bit.ly/Sen2CloudFree.

Functions
----------
  - exportCloudFreeSen2: Generate cloud free image for a specific
                         ee.ImageCollection.

How to cite
----------------
Schmitt, M., Hughes, L. H., Qiu, C., and Zhu, X. X.:
AGGREGATING CLOUD-FREE SENTINEL-2 IMAGES WITH GOOGLE
EARTH ENGINE, ISPRS Ann. Photogramm. Remote Sens.
Spatial Inf. Sci., IV-2/W7, 145–152,
https://doi.org/10.5194/isprs-annals-IV-2-W7-145-2019,
2019.
"""
import ee
import math
from preprocessing.s2.cloudscore import dilatedErossion, create_computeS2CloudScore
from preprocessing.s2.cloudstat import create_calcCloudStats
from preprocessing.s2.computeQualityScore import computeQualityScore_with_shadow, computeQualityScore_without_shadow

def exportCloudFreeSen2(ic, roi,
                        erosion=1.5,
                        dilation=3,
                        cloud_thresh=0.2,
                        shadow=True,
                        shadow_ndviThresh=-0.1,
                        shadow_irSumThresh=0.3,
                        mosaic_roicloudthresh=5
                        ):
    """ Generate a cloudfree Sentinel2 ImageCollection
    Args:
        ic (ee.ImageCollection): 
        roi (ee.Geometry.Polygon): Region Of Interest
        erosion (float): The erosion kernel radius.
        dilation (float): The dilation kernel radius.
        cloud_thresh: Ranges from 0-1. Parameter used in estimated the
                      cloud score, lower value will mask more pixels out. 
                      Generally 0.1-0.3 works well with 0.2 being used most
                      commonly.                              
        shadow_ndviThresh: The limit threshold to identify water bodies. By default -0.1
        shadow_irSumThresh: Sum of IR bands to include as shadows within TDOM and 
                            the shadow shift method (lower number masks out less)    
        mosaic_thresh:
    Returns:
        output (ee.Image): Original image with a band added.
        - shadowScore: A per-pixel score of shadowiness 
    """
    computeS2CloudScore = create_computeS2CloudScore(erosion=erosion, 
                                                     dilation=dilation,
                                                     cloud_thresh=cloud_thresh)
    calcCloudStats = create_calcCloudStats(roi = roi)    
    if shadow is True:                                                     
        projectShadows = create_projectShadows(erosion=erosion,
                                            dilation=dilation,
                                            shadow_ndviThresh=shadow_ndviThresh,
                                            shadow_irSumThresh=shadow_irSumThresh)
        ic = ee.ImageCollection("COPERNICUS/S2")\
            .filterBounds(roi)\
            .map(lambda x:x.clip(roi))\
            .map(computeS2CloudScore)\
            .map(calcCloudStats)\
            .map(projectShadows)\
            .map(computeQualityScore_with_shadow)\
            .sort('CLOUDY_PERCENTAGE_ROI')                
    else:
        ic = ee.ImageCollection("COPERNICUS/S2")\
            .filterBounds(roi)\
            .map(lambda x:x.clip(roi))\
            .map(computeS2CloudScore)\
            .map(calcCloudStats)\
            .map(computeQualityScore_without_shadow)\
            .sort('CLOUDY_PERCENTAGE_ROI')
    cloudFree = mergeCollection(
        ic = ic,
        mosaic_roicloudthresh = mosaic_roicloudthresh
        )
    return cloudFree

def mergeCollection(ic, mosaic_roicloudthresh):
    """Create best mosaic
    """    
    # Select the best images, which are below the cloud free threshold,
    # sort them in reverse order (worst on top) for mosaicing
    best = ic.filterMetadata('CLOUDY_PERCENTAGE_ROI', 'less_than', mosaic_roicloudthresh)\
             .sort('CLOUDY_PERCENTAGE_ROI',False)
    # Add the quality mosaic to fill in any missing areas of the ROI 
    # which aren't covered by good images             
    filtered = ic.qualityMosaic('cloudShadowScore')
    new_ic = ee.ImageCollection.fromImages([filtered, best.mosaic()])
    return ee.Image(new_ic.mosaic())
