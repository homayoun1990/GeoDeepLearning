import ee
from display.mapdisplay import mapdisplay, embed_map
from preprocessing.s2.cloudfreeS2 import exportCloudFreeSen2
ee.Initialize()

# viz params
saga_palette = ["#000180", "#0075FD", "#6CFB93", "#F99D05", "#A70700"]

cloudscore_viz = {'palette': saga_palette, 'max':0.2, 'min':0}
cloud_shadow_score_viz = {'palette': saga_palette, 'max':0, 'min':-0.2}
s2viz = {'bands':['B4','B3','B2'],'max':4000,'min':0}

# 1. Define study area (Loreto Peru)
roi=ee.Geometry.Polygon([[-73.36532592773438,-3.866994949019476],
                         [-73.24859619140625,-3.866994949019476],
                         [-73.24859619140625,-3.751208225628627],
                         [-73.36532592773438,-3.751208225628627],
                         [-73.36532592773438,-3.866994949019476]])

# 2. Create Sentinel-2 Image Collection
ic = ee.ImageCollection("COPERNICUS/S2")\
        .filterBounds(roi)\
        .filterDate('2019-01-01','2019-12-31')\
        .filterMetadata('CLOUDY_PIXEL_PERCENTAGE','less_than',0.1)

# 3. Intermediate masks
center = roi.centroid().coordinates().getInfo()
mosaic_s2 = exportCloudFreeSen2(ic,roi,shadow=False)

folium_obj = mapdisplay(center = center,
                        dicc = {
                                'sentinel2':mosaic_s2.getMapId(s2viz),
                                'cloud_percentage':mosaic_s2.select('cloudScore').getMapId(cloudscore_viz),
                                #'shadow_percentage':mosaic_s2.select('shadowScore').getMapId(cloudscore_viz),
                                #'cloudshadow_percentage':mosaic_s2.select('cloudShadowScore').getMapId(cloud_shadow_score_viz)
                                },
                        zoom_start=12)
embed_map(folium_obj)
#mosaic_s2.get('CLOUDY_PERCENTAGE_ROI').getInfo()