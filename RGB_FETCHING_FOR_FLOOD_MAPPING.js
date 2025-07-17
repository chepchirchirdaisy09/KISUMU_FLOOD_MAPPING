var before_start = '2017-07-15';
var before_end = '2019-08-10';
var after_start = '2019-08-10';
var after_end = '2023-03-23';

var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
          .filter(ee.Filter.date(before_start,before_end))
          .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',4))
          .filter(ee.Filter.bounds(roi));
          

print(s2.size());

var mosaic = s2.mosaic();

var clip = mosaic.clip(roi);

Map.addLayer(clip,{min:0.0,max:3000.0,bands:['B4','B3','B2']},'CLIP');
Map.centerObject(roi,13);

Export.image.toDrive({
  image: clip.select('B4','B3','B2'), 
  description: 'FLOOD_MAPPING_EXPORT', 
  folder: 'DAISY_PRIME', 
  fileNamePrefix: 'RGB_BEFORE', 
  region: roi, 
  scale: 10, 
  crs: 'EPSG:4326', 
  maxPixels: 1e10, 
  });
