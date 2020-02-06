#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This module generate a cloudScore band in Sentinel-2
images. It is a adaptation of the javascript code:
http://bit.ly/Sen2CloudFree.

Functions
----------
  - rescale: scale images considering a range
  - dilatedErossion: applied morphological opening
                   and closing to images.
  - computeS2CloudScore: Compute cloudScore in a Sentinel-2 ee.Image

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

def rescale(img, exp, thresholds):
  ''' Rescale values in an ee.Image

  Scales the input so that the range of input values [low, high]
  becomes [0,1]. Values outside the range are NOT clamped. This
  algorithm always produces floating point pixels.

  Args:
    img (ee.Image): The Image to scale
    exp ('str'): The expression to apply
    thresholds (list): The value mapped to 0 and 1.
  '''
  return img.expression(exp, {'img': img})\
            .unitScale(thresholds[0], thresholds[1])
     


def dilatedErossion(score, erosion=1.5, dilation=3):
  ''' Applied morphological opening and closing to images.

  The opening operation is applied first in order to
  remove single pixels with a high cloud score. The
  closing operation is then used to fill any holes
  which occur in areas with a high cloud score and
  to ensure edge regions of clouds are correctly scored

  Args:
    score (ee.Image): The image to applied morphological operators.
    erosion (float): The erosion kernel radius.
    dilation (float): The dilation kernel radius.
  ''' 
  # Perform opening on the cloud scores
  #.reproject(crs = 'EPSG:4326',
  #                       crsTransform = None,
  #                       scale = 20)\
  return score.focal_min(radius = erosion,
                         kernelType = 'circle',
                         iterations = 3)\
              .focal_max(radius = dilation,
                         kernelType = 'circle',
                         iterations =3)


def create_computeS2CloudScore(erosion=1.5, dilation=3, cloud_thresh=0.2):
  '''Compute cloudScore in a Sentinel-2 ee.Image
    Args:
      img (ee.Image): The Image to compute the cloudScore.
      erosion (float): The erosion kernel radius.
      dilation (float): The dilation kernel radius.
      cloud_thresh: Ranges from 0-1. Lower value will mask more pixels out. 
                    Generally 0.1-0.3 works well with 0.2 being used most
                    commonly.
  '''
  def computeS2CloudScore(img):
    toa = img.select(['^(B|QA60).*']).divide(10000)

    # Compute several indicators of cloudyness and take the minimum of them.
    score = ee.Image(1)

    # Clouds are reasonably bright in the blue and cirrus bands.
    # (img - Rmin)/(Rmax - Rmin) --> ee.Image.unitscale
    score = score.min(rescale(toa, 'img.B2', [0.1, 0.5]))  
    score = score.min(rescale(toa, 'img.B1', [0.1, 0.3]))
    score = score.min(rescale(toa, 'img.B1 + img.B10', [0.15, 0.2]))
    
    # Clouds are reasonably bright in all visible bands.
    score = score.min(rescale(toa, 'img.B4 + img.B3 + img.B2',
                              [0.2,0.8]))

    # Clouds are moist
    ndmi = img.normalizedDifference(['B8', 'B11'])
    score = score.min(rescale(ndmi, 'img', [-0.1, 0.1]))
    
    # However, clouds are not snow.
    ndsi = img.normalizedDifference(['B3', 'B11'])
    score = score.min(rescale(ndsi, 'img', [-0.8, 0.6]))
    
    # Clip the lower end of the score
    score = score.max(ee.Image(0.001))
  
    # Remove small regions and clip the upper bound
    dilated = dilatedErossion(score=score, erosion=dilation,
                              dilation=dilation).min(ee.Image(1.0))

    # score = score.multiply(dilated)
    score = score.reduceNeighborhood(reducer=ee.Reducer.mean(),
                                     kernel=ee.Kernel.square(5))    
    img = img.addBands(score.rename('cloudScore'))
    cloud_mask = img.select('cloudScore').gt(cloud_thresh).rename('cloudMask')
    img = img.addBands(cloud_mask)
    return img
  return computeS2CloudScore