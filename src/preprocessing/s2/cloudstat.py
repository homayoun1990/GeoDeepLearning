#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This module calculate cloud statistics for a specific 
Region Of Interest. It is a adaptation of the javascript
code: http://bit.ly/Sen2CloudFree.

Functions
----------
  - calcCloudStats: Calculates a mask for clouds in the image.

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
def create_calcCloudStats(roi):
  def calcCloudStats(img):
    ''' Calculates clouds statistic in an ee.Image.
    Args:  
      img (ee.Image): Image from image collection with a valid mask layer    
    Returns:  
      output (ee.Image): Original image with added stats.
        - CLOUDY_PERCENTAGE: The percentage of the image area affected by clouds
        - cloudScore: A per pixel score of cloudiness
    '''
    cloudAreaImg = img.select(['cloudMask']).multiply(ee.Image.pixelArea())   
    stats = cloudAreaImg.reduceRegion(
      reducer = ee.Reducer.sum(),
      geometry = roi,
      scale = 10,
      #maxPixels= 10**12,
      bestEffort = True)
    cloudPercentROI = ee.Number(stats.get('cloudMask'))\
                        .divide(roi.area())\
                        .multiply(100)
    img = img.set('CLOUDY_PERCENTAGE_ROI', cloudPercentROI)
    return img
  return calcCloudStats