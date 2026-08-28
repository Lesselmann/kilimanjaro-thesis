// ============================================================
// COMBINED GEE Export -- ALL Raster Layers (Master Script, FINAL)
// Kilimanjaro Thesis
// ============================================================
// Copy into the Google Earth Engine Code Editor:
// https://code.earthengine.google.com/
//
// Produces SEVEN export tasks in total:
//   1. kilimanjaro_ndvi_ndsi_ephemera.tif  -- NDVI_mean, NDSI, NDVI_var (3 bands)
//   2. kilimanjaro_canopy_height.tif       -- ETH Global Canopy Height (1 band)
//   3. kilimanjaro_lai.tif                 -- LAI_mean, LAI_var (2 bands)
//   4. kilimanjaro_rgb.tif                 -- dry-season RGB (3 bands)
//   5. kilimanjaro_rgb_seasonal.tif        -- 4 seasons x RGB (12 bands)
//   6. kilimanjaro_dsm_cop30.tif           -- Copernicus DEM GLO-30 DSM (1 band)
//   7. kilimanjaro_worldcover.tif          -- ESA WorldCover 10m v200 (1 band)
//
// Uses the UPDATED, wider bounding box covering the full 373-cell
// grid (including Community-Ecosystem cells) plus a 15km buffer
// for the 10km viewshed search radius.
//
// NOTE: Each export creates a SEPARATE task in the "Tasks" tab --
// you need to click "Run" on each of the 7 tasks individually.
// GEE does not support batch-confirming multiple exports at once.
//
// NOTE ON FABDEM: FABDEM (bare-earth) is NOT available in GEE --
// continue downloading it separately via the `fabdem` Python
// package, using the same bounds as below:
//   bounds = (36.8665, -3.5621, 37.8755, -2.7317)
// After downloading, it must be manually reprojected to EPSG:32737
// (see reproject_dsm_fabdem.R) -- unlike all GEE exports below,
// which are already in EPSG:32737 directly from the export itself.
// ============================================================

// ---------------------------------------------------------------
// 0. SHARED CONFIGURATION
// ---------------------------------------------------------------
var aoi = ee.Geometry.Rectangle([36.8665, -3.5621, 37.8755, -2.7317]);

var startYear = 2018;
var endYear = 2024;
var cloudThreshold = 20;

var seasons = [
  {name: 'dry1', startMonth: 6, endMonth: 10},   // main dry season
  {name: 'dry2', startMonth: 1, endMonth: 2},    // secondary dry season
  {name: 'wet1', startMonth: 3, endMonth: 5},    // main wet season
  {name: 'wet2', startMonth: 11, endMonth: 12}   // secondary wet season
];

Map.centerObject(aoi, 9);

function maskS2clouds(image) {
  var scl = image.select('SCL');
  var cloudMask = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  return image.updateMask(cloudMask).divide(10000)
              .copyProperties(image, ['system:time_start']);
}

function maskLaiQuality(image) {
  var qc = image.select('FparLai_QC');
  var goodQuality = qc.bitwiseAnd(3).eq(0);
  return image.updateMask(goodQuality);
}

// ================================================================
// 1. SENTINEL-2 NDVI / NDSI / NDVI_var
// ================================================================
var dryS2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(startYear + '-01-01', endYear + '-12-31')
  .filter(ee.Filter.calendarRange(6, 10, 'month'))
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloudThreshold))
  .map(maskS2clouds);

print('1. Number of scenes in the dry-season composite:', dryS2.size());

var dryComposite = dryS2.median();
var ndvi_mean = dryComposite.normalizedDifference(['B8', 'B4']).rename('NDVI_mean');
var ndsi = dryComposite.normalizedDifference(['B3', 'B11']).rename('NDSI');

var seasonalNDVIList = seasons.map(function(season) {
  var seasonCollection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi)
    .filterDate(startYear + '-01-01', endYear + '-12-31')
    .filter(ee.Filter.calendarRange(season.startMonth, season.endMonth, 'month'))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloudThreshold))
    .map(maskS2clouds);
  return seasonCollection.median().normalizedDifference(['B8', 'B4']).rename(season.name);
});

var ndvi_var = ee.Image.cat(seasonalNDVIList).reduce(ee.Reducer.stdDev()).rename('NDVI_var');
var sentinelStack = ndvi_mean.addBands(ndsi).addBands(ndvi_var).toFloat();

Map.addLayer(ndvi_mean, {min: -0.2, max: 0.9, palette: ['red', 'yellow', 'green']}, '1. NDVI mean');

Export.image.toDrive({
  image: sentinelStack,
  description: 'kilimanjaro_ndvi_ndsi_ephemera',
  folder: 'GEE_exports',
  region: aoi,
  scale: 10,
  crs: 'EPSG:32737',
  maxPixels: 1e10
});

// ================================================================
// 2. ETH GLOBAL CANOPY HEIGHT 2020
// ================================================================
var height = ee.Image('users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1')
  .clip(aoi)
  .rename('canopy_height_m');

print('2. Canopy Height Info:', height);
Map.addLayer(height, {min: 0, max: 40, palette: ['white', 'yellow', 'green', 'darkgreen']}, '2. Canopy Height');

Export.image.toDrive({
  image: height,
  description: 'kilimanjaro_canopy_height',
  folder: 'GEE_exports',
  region: aoi,
  scale: 10,
  crs: 'EPSG:32737',
  maxPixels: 1e10
});

// ================================================================
// 3. MODIS LAI_mean + LAI_var
// ================================================================
var laiCollection = ee.ImageCollection('MODIS/061/MCD15A3H')
  .filterBounds(aoi)
  .filterDate(startYear + '-01-01', endYear + '-12-31')
  .filter(ee.Filter.calendarRange(6, 10, 'month'))
  .map(maskLaiQuality);

print('3. Number of LAI scenes (dry season):', laiCollection.size());

var laiComposite = laiCollection.select('Lai').median().multiply(0.1).rename('LAI_mean');

var seasonalLaiList = seasons.map(function(season) {
  var seasonCollection = ee.ImageCollection('MODIS/061/MCD15A3H')
    .filterBounds(aoi)
    .filterDate(startYear + '-01-01', endYear + '-12-31')
    .filter(ee.Filter.calendarRange(season.startMonth, season.endMonth, 'month'))
    .map(maskLaiQuality);
  return seasonCollection.select('Lai').median().multiply(0.1).rename(season.name);
});

var laiVar = ee.Image.cat(seasonalLaiList).reduce(ee.Reducer.stdDev()).rename('LAI_var');
var laiStack = laiComposite.addBands(laiVar).toFloat();

Map.addLayer(laiComposite, {min: 0, max: 7, palette: ['white', 'yellow', 'green', 'darkgreen']}, '3. LAI mean');
Map.addLayer(laiVar, {min: 0, max: 2, palette: ['white', 'orange', 'red']}, '3. LAI variability');

Export.image.toDrive({
  image: laiStack,
  description: 'kilimanjaro_lai',
  folder: 'GEE_exports',
  region: aoi,
  scale: 500,
  crs: 'EPSG:32737',
  maxPixels: 1e10
});

// ================================================================
// 4. SENTINEL-2 RGB -- dry-season composite
// ================================================================
var rgb = dryComposite.select(['B4', 'B3', 'B2']).rename(['Red', 'Green', 'Blue']);
Map.addLayer(dryComposite, {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3}, '4. RGB True Color (dry season)');

Export.image.toDrive({
  image: rgb,
  description: 'kilimanjaro_rgb',
  folder: 'GEE_exports',
  region: aoi,
  scale: 10,
  crs: 'EPSG:32737',
  maxPixels: 1e10
});

// ================================================================
// 5. SENTINEL-2 RGB -- 4 seasonal composites (for ColorDiversity_var)
// ================================================================
var seasonalRgbBandsList = [];
seasons.forEach(function(season) {
  var seasonCollection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi)
    .filterDate(startYear + '-01-01', endYear + '-12-31')
    .filter(ee.Filter.calendarRange(season.startMonth, season.endMonth, 'month'))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloudThreshold))
    .map(maskS2clouds);

  var seasonComposite = seasonCollection.median();
  var seasonRgb = seasonComposite.select(['B4', 'B3', 'B2'])
    .rename([season.name + '_Red', season.name + '_Green', season.name + '_Blue']);
  seasonalRgbBandsList.push(seasonRgb);
});

var seasonalRgbStack = ee.Image.cat(seasonalRgbBandsList).toFloat();
print('5. Seasonal RGB stack bands:', seasonalRgbStack.bandNames());

Export.image.toDrive({
  image: seasonalRgbStack,
  description: 'kilimanjaro_rgb_seasonal',
  folder: 'GEE_exports',
  region: aoi,
  scale: 10,
  crs: 'EPSG:32737',
  maxPixels: 1e10
});

// ================================================================
// 6. COPERNICUS DEM GLO-30 (DSM)
// ================================================================
var dsm = ee.ImageCollection('COPERNICUS/DEM/GLO30')
  .filterBounds(aoi)
  .select('DEM')
  .mosaic()
  .clip(aoi)
  .rename('dsm_elevation_m');

print('6. DSM Info:', dsm);
Map.addLayer(dsm, {min: 700, max: 5900, palette: ['blue', 'green', 'yellow', 'brown', 'white']}, '6. DSM (Copernicus GLO-30)');

Export.image.toDrive({
  image: dsm,
  description: 'kilimanjaro_dsm_cop30',
  folder: 'GEE_exports',
  region: aoi,
  scale: 30,
  crs: 'EPSG:32737',
  maxPixels: 1e10
});

// ================================================================
// 7. ESA WORLDCOVER 10m v200
// ================================================================
var worldcover = ee.ImageCollection('ESA/WorldCover/v200')
  .first()
  .select('Map')
  .clip(aoi);

print('7. WorldCover Info:', worldcover);
Map.addLayer(worldcover, {}, '7. ESA WorldCover 10m v200');

Export.image.toDrive({
  image: worldcover,
  description: 'kilimanjaro_worldcover',
  folder: 'GEE_exports',
  region: aoi,
  scale: 10,
  crs: 'EPSG:32737',
  maxPixels: 1e10
});

// ---------------------------------------------------------------
// NOTES -- AFTER RUNNING THIS SCRIPT
// ---------------------------------------------------------------
// 1. Go to the "Tasks" tab (top right)
// 2. You will see 7 pending tasks:
//    - kilimanjaro_ndvi_ndsi_ephemera
//    - kilimanjaro_canopy_height
//    - kilimanjaro_lai
//    - kilimanjaro_rgb
//    - kilimanjaro_rgb_seasonal
//    - kilimanjaro_dsm_cop30
//    - kilimanjaro_worldcover
// 3. Click "Run" on EACH of the 7 tasks individually
// 4. All seven will appear in your Google Drive under "GEE_exports/"
//    once complete -- download each to your local machine
// 5. Still download FABDEM separately (see note at the top) --
//    it is not available as a GEE asset, and needs manual
//    reprojection to EPSG:32737 afterwards (unlike all 7 exports
//    above, which are already correctly projected)
// 6. Place all files in their respective local folders (DEM,
//    Indices, Satellite, Landcover) and re-run the R indicator
//    scripts (01-18) as needed
