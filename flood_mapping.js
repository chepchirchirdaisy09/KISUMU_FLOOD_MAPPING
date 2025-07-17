// ======================== DRAWN GEOMETRY USED ==========================
Map.centerObject(geometry3, 10);
// Map.addLayer(geometry, {color: 'grey'}, 'Drawn AOI');

// ======================== DEFINE DATES ================================
var beforeStart = '2017-07-15';
var beforeEnd   = '2019-08-10';
var afterStart  = '2019-08-10';
var afterEnd    = '2023-03-23';

// ======================== LOAD S1 (VV) DATA ============================
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(geometry)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .select('VV');

// Check available VV images
var beforeVVCol = s1.filterDate(beforeStart, beforeEnd);
var afterVVCol  = s1.filterDate(afterStart, afterEnd);

print('Before Flood VV Dates:', beforeVVCol.aggregate_array('system:index'));
print('After Flood VV Dates:', afterVVCol.aggregate_array('system:index'));

var beforeVV = beforeVVCol.mosaic().clip(g2);
var afterVV  = afterVVCol.mosaic().clip(g2);

// Handle missing data
beforeVV = ee.Algorithms.If(beforeVV.bandNames().size().gt(0), beforeVV, ee.Image(0).clip(g2));
afterVV  = ee.Algorithms.If(afterVV.bandNames().size().gt(0), afterVV, ee.Image(0).clip(g2));
beforeVV = ee.Image(beforeVV);
afterVV  = ee.Image(afterVV);

// ======================== SPECKLE FILTER FUNCTIONS =====================
function toNatural(img) {
  return ee.Image(10.0).pow(img.divide(10.0));
}
function toDB(img) {
  return ee.Image(img).log10().multiply(10.0);
}
function RefinedLee(img) {
  var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
  var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);
  var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
  var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);
  var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], 
                                [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);
  var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);
  var sample_mean = mean3.neighborhoodToBands(sample_kernel); 
  var sample_var = variance3.neighborhoodToBands(sample_kernel);
  var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));
  var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);
  var dir_mean = mean3;
  var dir_var = variance3;
  var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));
  var b = varX.divide(dir_var);
  var result = dir_mean.add(b.multiply(img.subtract(dir_mean)));
  return(result.arrayFlatten([['sum']]));
}

// Apply filters
var beforeFiltered = toDB(RefinedLee(toNatural(beforeVV)));
var afterFiltered  = toDB(RefinedLee(toNatural(afterVV)));

// ======================== FLOOD DETECTION =============================
var difference = afterFiltered.divide(beforeFiltered);
var floodThreshold = 1.0;
var flooded = difference.gt(floodThreshold).rename('Flooded').selfMask();

Map.addLayer(beforeFiltered, {min: -25, max: 0}, 'Before VV Filtered', false);
Map.addLayer(afterFiltered, {min: -25, max: 0}, 'After VV Filtered', false);
Map.addLayer(flooded, {min: 0, max: 1, palette: ['red']}, 'Flooded Areas');

// ======================== PERMANENT WATER MASK ========================
var gsw = ee.Image("JRC/GSW1_2/GlobalSurfaceWater").select('seasonality');
var permanentWater = gsw.gte(5).clip(g2);
var permanentWaterMask = permanentWater.unmask(0).not();
flooded = flooded.updateMask(permanentWaterMask);

// ======================== SLOPE MASKING ===============================
var slope = ee.Algorithms.Terrain(ee.Image("WWF/HydroSHEDS/03VFDEM").clip(g2)).select('slope');
var slopeMask = slope.lt(5);
flooded = flooded.updateMask(slopeMask);

// ======================== REMOVE NOISE (ISOLATED PIXELS) ==============
var connections = flooded.connectedPixelCount(25);
flooded = flooded.updateMask(connections.gt(8));

// ======================== SENTINEL-2 RGB FOR CONTEXT ==================
var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(g2)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 5));

var beforeS2 = s2.filterDate(beforeStart, beforeEnd).mean().clip(geometry3);
var afterS2  = s2.filterDate(afterStart, afterEnd).mean().clip(geometry3);

var visParams = {
  min: 0.0,
  max: 3000.0,
  bands: ['B4', 'B3', 'B2']
};

Map.addLayer(beforeS2, visParams, 'Sentinel-2 Before Flood');
Map.addLayer(afterS2, visParams, 'Sentinel-2 After Flood');

// ======================== EXPORT IMAGES ===============================
// Export Sentinel-2 Before
Export.image.toDrive({
  image: beforeS2.select('B4', 'B3', 'B2'),
  description: 'Kisumu_Sentinel2_Before_RGB_prime',
  folder: 'THE_RGB',
  fileNamePrefix: 'Kisumu_S2_Before_RGB_prime',
  region: geometry3,
  scale: 10,
  crs:'EPSG:4326',
  maxPixels: 1e10
});

// Export Sentinel-2 After
Export.image.toDrive({
  image: afterS2.select('B4', 'B3', 'B2'),
  description: 'Kisumu_Sentinel2_After_RGB_prime',
  folder: 'THE_RGB',
  fileNamePrefix: 'Kisumu_S2_After_RGB_prime',
  region: geometry3,
  scale: 10,
  crs:'EPSG:4326',
  maxPixels: 1e10
});

// Export Detected Flood Map
Export.image.toDrive({
  image: flooded,
  description: 'Kisumu_Flood_Map_VV',
  folder: 'EarthEngine',
  fileNamePrefix: 'Kisumu_Flooded_Area_VV',
  region: geometry3,
  crs:'EPSG:4326',
  scale: 10,
  maxPixels: 1e13
});

// Example: Exporting a FeatureCollection to a shapefile in Google Drive
// var featureCollection = ee.FeatureCollection(geometry2); // Replace with your FeatureCollection
// Export.table.toDrive({
//   collection: featureCollection,
//   description: 'ROI_EXPORT',
//   folder:'EarthEngine',
//   fileNamePrefix:'FINAL_ROI',
//   fileFormat: 'SHP',
// });
