#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This module generate a qualityScore band  (merging cloud and 
shadow score) in Sentinel-2 images. It is a adaptation of the
javascript code: http://bit.ly/Sen2CloudFree.

Functions
----------
  - computeQualityScore:Compute qualityScore in a Sentinel-2 ee.Image  

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

def computeQualityScore_with_shadow(img):
    '''Compute cloudShadowScore in a Sentinel-2 ee.Image
    Arg:
        img (ee.Image): The Image to compute the cloudShadowScore.
    '''
    score = img.select(['cloudScore']).max(img.select(['shadowScore']))
    score = score.reduceNeighborhood(
        reducer = ee.Reducer.mean(),
        kernel = ee.Kernel.square(5)
    )
    score = score.multiply(-1)
    return img.addBands(score.rename('cloudShadowScore'))

def computeQualityScore_without_shadow(img):
    '''Compute cloudShadowScore in a Sentinel-2 ee.Image
    Arg:
        img (ee.Image): The Image to compute the cloudShadowScore.
    '''
    score = img.select(['cloudScore'])
    score = score.multiply(-1)
    return img.addBands(score.rename('cloudShadowScore'))    