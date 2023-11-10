var geom = ee.Geometry.Polygon([[
            [103,24],
            [116,24],
            [116,32],
            [103,32],
            [103,24],]])
Map.addLayer(geom,{},'study area')
//------------------------------------------------------------data prepare---------------------------------------------------------------------------------------------------
var Climate_collection = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE').filter(ee.Filter.date('2001-01-01', '2004-12-31')).select(['pr','tmmn','tmmx','srad'])
var MOD09GA_collection_fir2y = ee.ImageCollection('MODIS/061/MOD09GA').filter(ee.Filter.date('2001-01-01', '2002-08-08'))
var MOD09GA_collection_lateryears = ee.ImageCollection('MODIS/061/MOD09GA').filter(ee.Filter.date('2002-08-08', '2004-12-31'))
var MYD09GA_collection = ee.ImageCollection('MODIS/061/MYD09GA').filter(ee.Filter.date('2002-08-08', '2004-12-31'))
var Landcover_collection = ee.ImageCollection('MODIS/061/MCD12Q1').select('LC_Type1')

function reprojectandclip(image){
  return image.reproject('EPSG:4326',null,500).clip(geom)
}

function solar_radiation_qc(image){
  var srad = image.select('srad').multiply(0.1).multiply(2.592)
  var tmmn = image.select('tmmn').multiply(0.1)
  var tmmx = image.select('tmmx').multiply(0.1)
  var pr = image.select('pr')
  var tmean = tmmn.add(tmmx).divide(2).rename('tmean')
  var datestr = image.get('system:time_start')
  var year = ee.Date(datestr).get('year')
  var month = ee.Date(datestr).get('month')
  var I_month = tmean.divide(5).pow(1.514).rename('I_month')
  
  return image.addBands(srad,null,true).addBands(tmmn,null,true).addBands(tmmx,null,true).addBands(pr,null,true).addBands(tmean).addBands(I_month).set('year',year).set('month',month).set('yearmonth',ee.Date(datestr).format('YYYYMM'))
}

var Climate_collection = Climate_collection.map(solar_radiation_qc).map(reprojectandclip)
var MOD09GA_collection_fir2y = MOD09GA_collection_fir2y.map(reprojectandclip)
var MOD09GA_collection_lateryears = MOD09GA_collection_lateryears.map(reprojectandclip)
var MYD09GA_collection = MYD09GA_collection.map(reprojectandclip)
var Landcover_collection = Landcover_collection.map(reprojectandclip)
//--------------------------------------------------------------------data prepare end-------------------------------------------------------------------------------------
//--------------------------------------------------------------------function start---------------------------------------------------------------------------------------
//---------------------------------------------------------------------bit extract-----------------------------------------------------------------------------------------
var bitwiseExtract = function(input, fromBit, toBit) {
  var maskSize = ee.Number(1).add(toBit).subtract(fromBit)
  var mask = ee.Number(1).leftShift(maskSize).subtract(1)
  return input.rightShift(fromBit).bitwiseAnd(mask)
}
//---------------------------------------------------------------------bit extract end-------------------------------------------------------------------------------------
//---------------------------------------------------------------------quality control-------------------------------------------------------------------------------------
function MOD09GA_QC_FC(image){
    var qc_500_MOD = image.select('state_1km_MOD')
    var cloud_state_MOD = bitwiseExtract(qc_500_MOD,0,1).eq(0)
    var cloud_shadow_MOD = bitwiseExtract(qc_500_MOD,2,2).eq(0)
    var cirrus_MOD = bitwiseExtract(qc_500_MOD,8,9).lte(1)
    var cloud_algorithm_MOD = bitwiseExtract(qc_500_MOD,10,10).eq(0)
    var mask_MOD = cloud_state_MOD.and(cloud_shadow_MOD).and(cirrus_MOD).and(cloud_algorithm_MOD)
    var MODBands = image.select('sur_refl_MOD_b0.').updateMask(mask_MOD).multiply(0.0001)
    
    var qc_500_MYD = image.select('state_1km_MYD')
    var cloud_state_MYD = bitwiseExtract(qc_500_MYD,0,1).eq(0)
    var cloud_shadow_MYD = bitwiseExtract(qc_500_MYD,2,2).eq(0)
    var cirrus_MYD = bitwiseExtract(qc_500_MYD,8,9).lte(1)
    var cloud_algorithm_MYD = bitwiseExtract(qc_500_MYD,10,10).eq(0)
    var mask_MYD = cloud_state_MYD.and(cloud_shadow_MYD).and(cirrus_MYD).and(cloud_algorithm_MYD)
    var MYDBands = image.select('sur_refl_MYD_b0.').updateMask(mask_MYD).multiply(0.0001)
    
    return image.addBands(MODBands,null,true).addBands(MYDBands,null,true)
}

function MOD09GA_QC(image){
  var qc_500 = image.select('state_1km')
  var cloud_state = bitwiseExtract(qc_500,0,1).eq(0)
  var cloud_shadow = bitwiseExtract(qc_500,2,2).eq(0)
  var cirrus = bitwiseExtract(qc_500,8,9).lte(1)
  var cloud_algorithm = bitwiseExtract(qc_500,10,10).eq(0)
  var mask = cloud_state.and(cloud_shadow).and(cirrus).and(cloud_algorithm)
  var Bands = image.select('sur_refl_b0.').updateMask(mask).multiply(0.0001)
  
  return image.addBands(Bands,null,true)
}

//---------------------------------------------------------------------quality control end---------------------------------------------------------------------------------
//---------------------------------------------------------------------inner join and combined imagecollection-------------------------------------------------------------
var innerJoin = ee.Join.inner()
var filterTimeEq = ee.Filter.equals({
    'leftField':'system:time_start',
    'rightField':'system:time_start'
})

function combineimage(feature){
  return ee.Image.cat(feature.get('primary'), feature.get('secondary'))
}
//---------------------------------------------------------------------inner join and combined imagecollection end --------------------------------------------------------
//-------------------------------------------------------------------------calculate terra and aqua average value----------------------------------------------------------
var Cal_terra_aqua_average_value = function(terra,aqua){
  var data = terra.unmask(0).add(aqua.unmask(0))
  var mask = terra.mask().add(aqua.mask())
  var mask_for_multiply = mask.where(mask.eq(0),10).pow(-1)
  var mask_for_update = mask.where(mask.gt(0),1)
  var average_value = data.multiply(mask_for_multiply)
  var average_value = average_value.updateMask(mask_for_update)
  
  return average_value
}
//--------------------------------------------------------------------- calculate terra and aqua average value end----------------------------------------------------------
//----------------------------------------------------------function end-----------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------MOD09 process-----------------------------------------------------------------------------------
var MOD09GA = MOD09GA_collection_lateryears.select(['sur_refl_b01','sur_refl_b02','sur_refl_b03','sur_refl_b04','sur_refl_b05','sur_refl_b06','sur_refl_b07','state_1km'],['sur_refl_MOD_b01','sur_refl_MOD_b02','sur_refl_MOD_b03','sur_refl_MOD_b04','sur_refl_MOD_b05','sur_refl_MOD_b06','sur_refl_MOD_b07','state_1km_MOD'])
var MYD09GA = MYD09GA_collection.select(['sur_refl_b01','sur_refl_b02','sur_refl_b03','sur_refl_b04','sur_refl_b05','sur_refl_b06','sur_refl_b07','state_1km'],['sur_refl_MYD_b01','sur_refl_MYD_b02','sur_refl_MYD_b03','sur_refl_MYD_b04','sur_refl_MYD_b05','sur_refl_MYD_b06','sur_refl_MYD_b07','state_1km_MYD'])

var innerJoinedMOD09 = innerJoin.apply(MOD09GA, MYD09GA, filterTimeEq)
var MODIS_09 = ee.ImageCollection(innerJoinedMOD09.map(function(feature){
  return ee.Image.cat(feature.get('primary'), feature.get('secondary'))
}))

var MODIS09GA_qc = MODIS_09.map(MOD09GA_QC_FC).select(['sur_refl_MOD_b01','sur_refl_MOD_b02','sur_refl_MOD_b03','sur_refl_MOD_b04','sur_refl_MOD_b05','sur_refl_MOD_b06','sur_refl_MOD_b07','state_1km_MOD','sur_refl_MYD_b01','sur_refl_MYD_b02','sur_refl_MYD_b03','sur_refl_MYD_b04','sur_refl_MYD_b05','sur_refl_MYD_b06','sur_refl_MYD_b07'])
var MOD09GA_collection_fir2y_qc = MOD09GA_collection_fir2y.map(MOD09GA_QC).select(['sur_refl_b01','sur_refl_b02','sur_refl_b03','sur_refl_b04','sur_refl_b05','sur_refl_b06','sur_refl_b07'],['MODIS_Band1','MODIS_Band2','MODIS_Band3','MODIS_Band4','MODIS_Band5','MODIS_Band6','MODIS_Band7'])

function sameday_average(image){
  var MOD_Band1 = image.select('sur_refl_MOD_b01')
  var MYD_Band1 = image.select('sur_refl_MYD_b01')
  var MOD_Band2 = image.select('sur_refl_MOD_b02')
  var MYD_Band2 = image.select('sur_refl_MYD_b02')
  var MOD_Band3 = image.select('sur_refl_MOD_b03')
  var MYD_Band3 = image.select('sur_refl_MYD_b03')
  var MOD_Band4 = image.select('sur_refl_MOD_b04')
  var MYD_Band4 = image.select('sur_refl_MYD_b04')
  var MOD_Band5 = image.select('sur_refl_MOD_b05')
  var MYD_Band5 = image.select('sur_refl_MYD_b05')
  var MOD_Band6 = image.select('sur_refl_MOD_b06')
  var MYD_Band6 = image.select('sur_refl_MYD_b06')
  var MOD_Band7 = image.select('sur_refl_MOD_b07')
  var MYD_Band7 = image.select('sur_refl_MYD_b07')
  
  var averaged_Band1 = Cal_terra_aqua_average_value(MOD_Band1,MYD_Band1)
  var averaged_Band2 = Cal_terra_aqua_average_value(MOD_Band2,MYD_Band2)
  var averaged_Band3 = Cal_terra_aqua_average_value(MOD_Band3,MYD_Band3)
  var averaged_Band4 = Cal_terra_aqua_average_value(MOD_Band4,MYD_Band4)
  var averaged_Band5 = Cal_terra_aqua_average_value(MOD_Band5,MYD_Band5)
  var averaged_Band6 = Cal_terra_aqua_average_value(MOD_Band6,MYD_Band6)
  var averaged_Band7 = Cal_terra_aqua_average_value(MOD_Band7,MYD_Band7)
  
  return image.addBands(averaged_Band1.rename('MODIS_Band1')).addBands(averaged_Band2.rename('MODIS_Band2')).
            addBands(averaged_Band3.rename('MODIS_Band3')).addBands(averaged_Band4.rename('MODIS_Band4')).
            addBands(averaged_Band5.rename('MODIS_Band5')).addBands(averaged_Band6.rename('MODIS_Band6')).
            addBands(averaged_Band7.rename('MODIS_Band7')).toFloat()
}

var MODIS09GA_average = MODIS09GA_qc.map(sameday_average).select(['MODIS_Band1','MODIS_Band2','MODIS_Band3','MODIS_Band4','MODIS_Band5','MODIS_Band6','MODIS_Band7'])

function avoid_data_type_error(image){
  return image.toFloat()
}

var MOD09GA_collection_fir2y_qc = MOD09GA_collection_fir2y_qc.map(avoid_data_type_error)
var MODIS09GA_Daily_averaged = MOD09GA_collection_fir2y_qc.merge(MODIS09GA_average)
//print(MODIS09GA_Daily_averaged.limit(50))
//--------------------------------------------------------------------filling gap year-----------------------------------------------------------------------------------
var landcover2022 = ee.Image('MODIS/061/MCD12Q1/2021_01_01').select('LC_Type1').
                      set('system:time_start',1640995200000).set('system:time_end',1672531200000)//.set()
var landcover_final = Landcover_collection.merge(landcover2022)
//--------------------------------------------------generate LUE NDVImin NDVImax SRmin SRmax--------------------------------------------------------------------------
function generate_LUE(image){
  var image = image.where(image.eq(3),0.485).where(image.eq(1),0.389).
                where(image.eq(4),0.692).where(image.eq(2),0.985).
                where(image.eq(5),0.768).where(image.eq(6),0.429).
                where(image.eq(7),0.429).where(image.eq(8),0.542).
                where(image.eq(9),0.542).where(image.eq(10),0.542).
                where(image.eq(11),0.542).where(image.eq(12),0.542).
                where(image.eq(13),0.542).where(image.eq(14),0.542).
                where(image.eq(15),0.542).where(image.eq(16),0.542).where(image.eq(17),0.542)
  return image.toFloat()
}

function set_NDVImin(image){
  var image = image.where(image.gt(0),0.023)
  return image.toFloat()
}

function set_NDVImax(image){
  var image = image.where(image.eq(3),0.738).where(image.eq(1),0.647).
                where(image.eq(2),0.676).where(image.eq(4),0.747).
                where(image.eq(6),0.636).where(image.eq(7),0.636).
                where(image.eq(8),0.636).where(image.eq(9),0.636).
                where(image.gt(9),0.634).where(image.eq(5),0.702)
  return image.toFloat()
}

function set_SRmin(image){
  var image = image.where(image.gt(0),1.05)
  return image.toFloat()
}

function set_SRmax(image){
  var image = image.where(image.eq(3),6.63).where(image.eq(1),4.67).
                where(image.eq(2),5.17).where(image.eq(4),6.91).
                where(image.eq(6),4.49).where(image.eq(7),4.49).
                where(image.eq(8),4.49).where(image.eq(9),4.49).
                where(image.gt(9),4.44).where(image.eq(5),5.85)
  return image.toFloat()
}

var LUE = landcover_final.map(generate_LUE).select(['LC_Type1'],['LUE'])
var NDVImin = landcover_final.map(set_NDVImin).select(['LC_Type1'],['NDVI_min'])
var NDVImax = landcover_final.map(set_NDVImax).select(['LC_Type1'],['NDVI_max'])
var SRmin = landcover_final.map(set_SRmin).select(['LC_Type1'],['SR_min'])
var SRmax = landcover_final.map(set_SRmax).select(['LC_Type1'],['SR_max'])
//-------------------------------------------------------------------------VI------------------------------------------------------------------------------------------------------
function Cal_index(image){
  var red = image.select('MODIS_Band1')
  var Nir = image.select('MODIS_Band2')
  var Blue = image.select('MODIS_Band3')
  var Green = image.select('MODIS_Band4')
  var Swir = image.select('MODIS_Band6')
  
  var NDVI = Nir.subtract(red).divide(Nir.add(red))
  var NDVI = NDVI.where(NDVI.gt(1),1).where(NDVI.lt(-1),-1)
  var NDWI = Green.subtract(Nir).divide(Green.add(Nir))
  var NDWI = NDWI.where(NDWI.gt(1),1).where(NDWI.lt(-1),-1)
  var NDBI = Swir.subtract(Nir).divide(Swir.add(Nir))
  var NDBI = NDVI.where(NDBI.gt(1),1).where(NDBI.lt(-1),-1)
  var NIRV = NDVI.multiply(Nir)
  var SR = ee.Image(1).add(NDVI).divide(ee.Image(1).subtract(NDVI))
  
  return image.addBands(NDVI.rename('NDVI')).addBands(NDWI.rename('NDWI')).addBands(NDBI.rename('NDBI')).addBands(NIRV.rename('NIRV')).addBands(SR.rename('SR'))
}
//Map.addLayer(MODIS09GA_Daily_averaged.select('MODIS_Band2').first(),{},'test')
var MODIS09_indices = MODIS09GA_Daily_averaged.map(Cal_index).select(['NDVI','NDWI','NDBI','NIRV','SR'])

//var Indices_LST = Indices_LST.filter(ee.Filter.date('2001-01-01', '2005-12-31'))
//-------------------------------------------------------------------combine NDVImin max SRmin max---------------------------------------------------------------------------------
function setyear(image){
  var year_str = image.get('system:time_start')
  var year = ee.Date(year_str).get('year')
  return image.set({'year':year})
}
var LUE = LUE.map(setyear)
var NDVImin = NDVImin.map(setyear)
var NDVImax = NDVImax.map(setyear)
var SRmin = SRmin.map(setyear)
var SRmax = SRmax.map(setyear)
var MODIS09_indices = MODIS09_indices.map(setyear)

var filterTimeEqminmax = ee.Filter.equals({
    'leftField':'year',
    'rightField':'year'
})
var innerJoinedIL_LUE = innerJoin.apply(MODIS09_indices, LUE, filterTimeEqminmax)
var Indices_LUE = ee.ImageCollection(innerJoinedIL_LUE.map(function(feature){
  return ee.Image.cat(feature.get('primary'), feature.get('secondary'))
}))
var innerJoinedILL_Nm = innerJoin.apply(Indices_LUE, NDVImin, filterTimeEqminmax)
var Indices_LLM = ee.ImageCollection(innerJoinedILL_Nm.map(function(feature){
  return ee.Image.cat(feature.get('primary'), feature.get('secondary'))
}))
var innerJoinedILL_NN = innerJoin.apply(Indices_LLM, NDVImax, filterTimeEqminmax)
var Indices_LLM = ee.ImageCollection(innerJoinedILL_NN.map(function(feature){
  return ee.Image.cat(feature.get('primary'), feature.get('secondary'))
}))
var innerJoinedILL_NNS = innerJoin.apply(Indices_LLM, SRmin, filterTimeEqminmax)
var Indices_LLM = ee.ImageCollection(innerJoinedILL_NNS.map(function(feature){
  return ee.Image.cat(feature.get('primary'), feature.get('secondary'))
}))
var innerJoinedILL_NNSR = innerJoin.apply(Indices_LLM, SRmax, filterTimeEqminmax)
var Indices_LLM = ee.ImageCollection(innerJoinedILL_NNSR.map(function(feature){
  return ee.Image.cat(feature.get('primary'), feature.get('secondary'))
}))
//print(innerJoinedILL_NNSR)
//print(Indices_LLM)
//print('NDVI_LST',NDVI_LST)
//----------------------------------------------------------------Monthly max composition------------------------------------------------------------------------------------------
var start = ee.Date('2001-01-01')
//print(start.millis())
var end = ee.Date('2005-01-31')

var n_months = end.difference(start, 'month').subtract(1)
var months = ee.List.sequence(0, n_months).map(function (n) {return start.advance(n, 'month')})
//------------------------------------------------------mONTHLY MEAN ndvilst-------------------------------------------
print(Climate_collection)
var T = Climate_collection.filter(ee.Filter.eq('month', 7)).select(['tmean'],['T'])
function getTopt(image){
  var tmean = image.select('T')
  
  var Topt_mean = ee.Number(tmean.reduceRegion({
            reducer: ee.Reducer.mean(),
            geometry :geom,
            scale: 500,
            maxPixels: 1e9
            }).values().get(0));
  
  var Topt = ee.Image(Topt_mean).rename('Topt')

  return image.addBands(Topt).clip(geom)
}
var T = T.map(getTopt).select('Topt')


var I = Climate_collection.select('I_month')
var start = ee.Date('2001-01-01')
var end = ee.Date('2005-01-31')
var n_years = end.difference(start, 'year').subtract(1)
var years = ee.List.sequence(0, n_years).map(function (n){return start.advance(n, 'year')})
function makeYearlyComposite_I(date){
  var date = ee.Date(date)
  return I.filterDate(date, date.advance(1, 'year')).sum().set('system:index',date.format("YYYY")).set('year',date.get('year')).set('system:time_start',date.millis()).set('system:time_end',date.advance(1, 'year').millis())
}
var I = ee.ImageCollection.fromImages(years.map(makeYearlyComposite_I)).select(['I_month'],['I'])


function makeMonthlyComposite_MOD09(date){
  var date = ee.Date(date)
  return Indices_LLM.filterDate(date, date.advance(1, 'month')).max().set("system:index", date.format("YYYYMM")).set("year", date.format("YYYY")).set("system:time_start", date.millis()).set("system:time_end", date.advance(1, 'month').millis())
}
var Indices_LLM = ee.ImageCollection.fromImages(months.map(makeMonthlyComposite_MOD09))
var MODIS_SR = Climate_collection.combine(Indices_LLM)
var innerJoinedTopt = innerJoin.apply(MODIS_SR, T, filterTimeEqminmax)
var MODIS_SR = ee.ImageCollection(innerJoinedTopt.map(function (feature){return ee.Image.cat(feature.get('primary'),feature.get('secondary'))}))

var innerJoinedClim = innerJoin.apply(MODIS_SR, I, filterTimeEqminmax)
var MODIS_SR = ee.ImageCollection(innerJoinedClim.map(function (feature){return ee.Image.cat(feature.get('primary'),feature.get('secondary'))}))

function Calstress(image){
  var I = image.select('I')
  var tmean = image.select('tmean')
  var pre = image.select('pr')
  var Topt = image.select('Topt')
  
  var a = ee.Image().expression(
    '(0.6751 * I ** 3 - 77.1 * I ** 2 + 17920 * I + 482390) / 1000000',{I:I}
    ).rename('a')
  var Ep0 = ee.Image().expression(
    '16 * (10 * d / b) ** c',{d:tmean,b:I,c:a}
    ).rename('Ep0')
  var Rn = ee.Image().expression(
    '((a * b) ** 0.5) * (0.369 + 0.598 * (a / b) ** 0.5)',{a:Ep0,b:pre}
    )
  var EET = ee.Image().expression(
    'a * b * (a ** 2 + b ** 2 + a + b) / ((a + b) * (a ** 2 + b ** 2))',{a:pre,b:Rn}
    )
  var PET = ee.Image().expression(
    '(a + b) / 2',{a:EET,b:Ep0}
    )
  var Wstress = ee.Image().expression(
    '0.5 + 0.5 * (a / b)',{a:EET,b:PET}
    ).rename('Wstress')
  var Tstress1 = ee.Image().expression(
    '0.8 + 0.02 * a - 0.0005 * a ** 2',{a:Topt}
    )
  var Tstress1 = Tstress1.where(tmean.lt(-10),0).rename('Tstress1')
  var Tsigma2 = ee.Image().expression(
    '1.184 / (1 + exp(0.2 * (a - 10 - b))) * 1 / (1 + exp(0.3 * (b - a - 10)))',{a:Topt,b:tmean}
    )
  var Tsigma3 = ee.Image(1.184).divide(ee.Image(1).add(ee.Image(-2).exp())).multiply(ee.Image(1).divide(ee.Image(1).add(ee.Image(-3).exp())))
  
  var diff = tmean.subtract(Topt)
  var diff1 = diff.where(diff.lt(-13),0).where(diff.gte(-13).and(diff.lte(10)),1).where(diff.gt(10),0)
  var diff2 = diff.where(diff.gte(-13).and(diff.lte(10)),1).where(diff.gt(10),0).where(diff.lt(-13),0)
  var Tstress2 = Tsigma2.multiply(diff1).add(Tsigma3.multiply(diff2)).rename('Tstress2')
  var date = ee.Date(image.get('system:time_start'))

  return image.addBands(Wstress).addBands(Tstress1).addBands(Tstress2).clip(geom)//.addBands(Ep0).clip(geom)
}

var MODIS_SR = MODIS_SR.map(Calstress)

function CalFPAR(image){
  var NDVI = image.select('NDVI')
  var SR = image.select('SR')
  var NDVImin = image.select('NDVI_min')
  var NDVImax = image.select('NDVI_max')
  var SRmin = image.select('SR_min')
  var SRmax = image.select('SR_max')
  
  var FPAR1 = NDVI.subtract(NDVImin).divide(NDVImax.subtract(NDVImin)).multiply(ee.Image(0.949)).add(0.001)
  var FPAR2 = SR.subtract(SRmin).divide(SRmax.subtract(SRmin)).multiply(ee.Image(0.949)).add(0.001)
  var FPAR = FPAR1.multiply(0.5).add(FPAR2.multiply(0.5))
  var FPAR = FPAR.where(FPAR.gt(0.95),0.95).where(FPAR.lt(0.05),0.05).rename('FPAR')
  
  return image.addBands(FPAR).toFloat()
}
var MODIS_SR = MODIS_SR.map(CalFPAR)

function CalAPAR(image){
  var FPAR = image.select('FPAR')
  var Sol = image.select('srad')
  
  
  var APAR = FPAR.multiply(Sol).multiply(ee.Image(0.5)).rename('APAR')
  var LUE = image.select('LUE')
  var Wstress = image.select('Wstress')
  var Tstress1 = image.select('Tstress1')
  var Tstress2 = image.select('Tstress2')
  var yupson = Tstress1.multiply(Tstress2).multiply(Wstress).multiply(LUE).rename('yupson')
  var NPP = APAR.multiply(yupson).rename('NPP')
  
  return image.addBands(NPP).addBands(yupson).addBands(APAR).toFloat()
}
var MODIS_SR = MODIS_SR.map(CalAPAR)
print(MODIS_SR)
var gppVis = {
  min: 0.0,
  max: 2000.0,
  palette: ['bbe029', '0a9501', '074b03'],
};
var nppVis = {
  min: 0.0,
  max: 10.0,
  palette: ['bbe029', '0a9501', '074b03'],
};
//Map.addLayer(MODIS_SR.select('NPP').filter(ee.Filter.date('2001-01-01', '2002-01-01')).sum(),gppVis,'NPP')
// Map.addLayer(MODIS_SR.select('yupson').filter(ee.Filter.date('2001-01-01', '2002-01-01')).sum(),nppVis,'yupson')
// Map.addLayer(MODIS_SR.select('NPP').filter(ee.Filter.date('2001-01-01', '2002-01-01')).sum(),gppVis,'NPP')
var NPP_month = MODIS_SR.select('NPP')
print(NPP_month)

function makeYearlyComposite_NPP(date){
  var date = ee.Date(date)
  return NPP_month.filterDate(date, date.advance(1, 'year')).sum().set('system:index',date.format("YYYY")).set('year',date.get('year')).set('system:time_start',date.millis()).set('system:time_end',date.advance(1, 'year').millis())
}
var NPP_Year = ee.ImageCollection.fromImages(years.map(makeYearlyComposite_NPP))
Map.addLayer(NPP_Year,gppVis,'NPP')