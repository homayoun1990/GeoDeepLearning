#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Author: Gennadii Donchyts
License: Apache 2.0

This module generate a shadowScore band in Sentinel-2
images. It is a adaptation of the javascript code:
http://bit.ly/Sen2CloudFree.

Functions
----------
- projectShadows: Compute shadowScore in a Sentinel-2 ee.Image

How to cite
----------------
Schmitt, M., Hughes, L. H., Qiu, C., and Zhu, X. X.:
AGGREGATING CLOUD-FREE SENTINEL-2 IMAGES WITH GOOGLE
EARTH ENGINE, ISPRS Ann. Photogramm. Remote Sens.
Spatial Inf. Sci., IV-2/W7, 145–152,
https://doi.org/10.5194/isprs-annals-IV-2-W7-145-2019,
2019.
"""

import math
import ee
from preprocessing.s2.cloudscore import dilatedErossion

def create_projectShadows(erosion=1.5,
                          dilation=3,
                          shadow_ndviThresh=-0.1,
                          shadow_irSumThresh=0.3):
  ''' Compute shadowScore in a Sentinel-2 ee.Image
  Args:
    img (ee.Image): The Image to compute the shadowScore.
    shadow_ndviThresh: The limit threshold to identify water bodies. By default -0.1
    shadow_irSumThresh: Sum of IR bands to include as shadows within TDOM and 
                        the shadow shift method (lower number masks out less)
  
  Returns:
    output (ee.Image): Original image with a band added.
      - shadowScore: A per-pixel score of shadowiness 
  '''
  def projectShadows(img):
    cloudHeights = ee.List.sequence(200,10000,250)
    meanAzimuth = img.get('MEAN_SOLAR_AZIMUTH_ANGLE')
    meanZenith = img.get('MEAN_SOLAR_ZENITH_ANGLE')

    cloudMask = img.select(['cloudMask'])

    #Find dark pixels
    darkPixelsImg = img.select(['B8','B11','B12'])\
                      .divide(10000)\
                      .reduce(ee.Reducer.sum())

    ndvi = img.normalizedDifference(['B8','B4'])
    waterMask = ndvi.lt(shadow_ndviThresh)
    darkPixels = darkPixelsImg.lt(shadow_irSumThresh)

    # Get the mask of pixels which might be shadows excluding water
    darkPixelMask = darkPixels.And(waterMask.Not())
    darkPixelMask = darkPixelMask.And(cloudMask.Not())

    # Find where cloud shadows should be based on solar geometry
    # Convert to radians
    azR = ee.Number(meanAzimuth).add(180).multiply(math.pi).divide(180.0)
    zenR = ee.Number(meanZenith).multiply(math.pi).divide(180.0)

    # Find the shadows
    def fxshadows(cloudHeight):
      cloudHeight = ee.Number(cloudHeight)
      # Distance shadow is cast
      shadowCastedDistance = zenR.tan().multiply(cloudHeight)
      # X distance of shadow
      x = azR.sin().multiply(shadowCastedDistance).multiply(-1)
      # Y distance of shadow
      y = azR.cos().multiply(shadowCastedDistance).multiply(-1)
      return img.select(['cloudScore'])\
                .displace(ee.Image.constant(x).addBands(ee.Image.constant(y)))

    shadows = cloudHeights.map(fxshadows)  
    shadowMasks = ee.ImageCollection.fromImages(shadows)
    shadowMask = shadowMasks.mean()

    # Create shadow mask
    shadowMask = dilatedErossion(shadowMask.multiply(darkPixelMask))
    shadowScore = shadowMask.reduceNeighborhood(
      reducer = ee.Reducer.max(),
      kernel = ee.Kernel.square(1))
    img = img.addBands(shadowScore.rename(['shadowScore']))
    return img
  return projectShadows