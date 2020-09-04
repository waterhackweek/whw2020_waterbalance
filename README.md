# Water Hack Week 2020: Irrigation Water Balance

## Collaborators 
- Project Lead: Ayman Alafifi
- Data Science Lead: Scott Black
- Science Team: Hemendra Kumar, Takhellambam Bijoychandra Singh,Parisa Shahbaz, Brooke Mason, Santosh Subhash Palmate, Deen Dayal, Takhellambam Bijoychandra Singh, Caitlin Eger, Bart Nijssen 

## Background Statement
- Precipitation, and in particular its effective portion, provides part of the water crops need to satisfy their transpiration requirements. The soil, acting as a buffer, stores part of the precipitation water and returns it to the crops in times of deficit. 
- In humid climates, this mechanism is sufficient to ensure satisfactory growth in rainfed agriculture. In arid climates or during extended dry seasons, irrigation is necessary to compensate for the evapotranspiration (crop transpiration and soil evaporation) deficit due to insufficient or erratic precipitation.
- Irrigation consumptive water use is defined as the volume of water needed to compensate for the deficit between potential evapotranspiration on the one side and effective precipitation over the crop growing period and change in soil moisture content on the other side. It varies considerably with climatic conditions, seasons, crops and soil types.

## Project Goal
Map irrigation consumptive water use in the Lower Yakima Basin, Washington using the crop water demand equation:

ICU = ETc - P – DS

where:
ICU = irrigation consumptive water use needed to satisfy crop water demand (mm)
ETc = potential crop evapotranspiration (mm)
P = effective precipitation (mm)
DS = change in soil moisture (mm)

## Data Sources
- Potential crop evapotranspiration data:
- Preciptation data: PM_3IMERGDF: GPM IMERG Final Precipitation L3 1 day 0.1 degree x 0.1 degree V06
https://disc.gsfc.nasa.gov/datasets/GPM_3IMERGDF_06/summary?keywords=imerg
- Crop data: USDA National Agricultural Statistics Service CropScape - Cropland Data Layer
https://nassgeodata.gmu.edu/CropScape/
- Crop coefficients: Determining Crop Water Use - Crop Coefficients
https://farmwest.com/node/932
- Soil moisture data: NLDAS_MOS0125_M: NLDAS Mosaic Land Surface Model L4 Monthly 0.125 x 0.125 degree V002 https://disc.gsfc.nasa.gov/datasets/NLDAS_MOS0125_M_002/summary?keywords=NLDAS_MOS0125_M_002

## Python Packages
- folium
- numpy
- pandas
- geopandas
- pyproj
- xarray
- rioxarray
- regionmask
- shapely 
- pynhd
- pygeoogc
- pygeoutils
- hydrodata 
- datetime
- warnings
- netCDF4 as nc
- pathlib
- contextily
- dask
- pygeoogc
-requests

## Additional Readings
- AQUASTAT - FAO's Global Information System on Water and Agriculture: http://www.fao.org/aquastat/en/data-analysis/irrig-water-use/irrig-water-requirement
- Chapter 2 - FAO Penman-Monteith equation: http://www.fao.org/3/X0490E/x0490e06.htm
