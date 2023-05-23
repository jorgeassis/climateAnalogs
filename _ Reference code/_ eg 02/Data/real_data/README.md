## Files

* HOT_bottle_data_mod.txt - Data from Hawaii Ocean Time-series Station Aloha, missing or bad data coded as -9, headers were cleaned post download. Downloaded from https://hahana.soest.hawaii.edu/hot/hot-dogs/bextraction.html
	* columns (units in parentheses):  
	"botid" - unique sample identified  
	"date(mmddyy)" - data of collection  
	"time(hhmmss)" - time of day sample collected  
	"press(dbar)" - pressure  
	"theta(ITS-90)" - potential temperature  
	"sigma(kg/m3)" - potential density  
	"temp(ITS-90)" - CTD measured temperature  
	"csal(PSS-78)" - CTD salinity  
	"coxy(umol/kg)" - CTD oxygen  
	"bsal(PSS-78)" - bottle salinity  
	"boxy(umol/kg)" - bottle oxygen  
	"dic(umol/kg)" - dissolved inorganic carbon  
	"ph" - pH  
	"alk(ueq/kg)" - alkalinity  
	"phos(umol/k" - phosphate  

* BBH.csv - Sea surface temperature data from Booth Bay Harbour, Maine, USA. Downloaded from https://www.maine.gov/dmr/science-research/weather-tides/bbhenv.html
	* columns
	"COLLECTION_DATE" - date of collection  
	"Sea Surface Temp Max C" - maximum measured sea surface temperature in degrees Celsius  
	"Sea Surface Temp Min C" - minimum measured sea surface temperature in degrees Celsius  
	"Sea Surface Temp Ave C" - mean measured sea surface temperature in degrees Celsius  

* CML.csv - Sea surface temperature and pH data from the Coastal Marine Laboratory, New Hampshire, USA. Downloaded from http://www.neracoos.org/erddap/tabledap/UNH_CML.html
	* columns 
	"time_UTC" - collection time  
	"depth_m" - sampling depth  
	"latitude_dN" - latitude in decimal degrees North  
	"longitude_dE" - longitude in decimal degrees East  
	"temperature_C" - water temperature in degrees Celsius  
	"pH" - pH 

* NH70W_43N_20##.csv - Sea surface temperature and pH data collected from NOAA mooring NH_70W_43N. Data collected for each year between 2011 and 2016 was downloaded as a separate file. Data downloaded from https://www.nodc.noaa.gov/archive/arc0062/0115402/8.8/data/0-data/ 
	* column (units in parentheses)
	"Mooring Name" - unique ID of mooring where sample was taken
	"Latitude" 
	"Longitude"
	"Date" - sampling date  
	"Time" - sampling time   
	"xCO2  SW (wet) (umol/mol)" - mole fraction of CO2 in air in equilibrium with the sea water at sea surface temperature and measured humidity  
	"CO2 SW QF" - quality flag for xCO2 sw (wet)  
	"H2O SW (mmol/mol)" - mole fraction of H2O in air from equilibrator
	"xCO2  Air (wet) (umol/mol)  
	"CO2 Air QF" - quality flag for xCO2 Air (wet)  
	"H2O Air (mmol/mol)" - mole fraction of H2O in air from airblock, 4 feet above the sea surface  
	"Licor Atm Pressure (hPa)" - atmospheric pressure at the airblock, 4 feet above the seasurface  
	"Licor Temp (C)" - temperature of the Infrared Licor 820  
	"MAPCO2 %O2" - percent oxygen of the surface seawater divided by the percent oxygen of the atmosphere at 4 feet above the sea surface
	"SST (C)" - sea surface temperature  
	"Salinity (PSU)" - sea surface salinity
	"xCO2  SW (dry) (umol/mol)" - mole fraction of CO2 in air in equilibrium with the seawaterat sea surface temperature (dry air)  
	"xCO2  Air (dry) (umol/mol)" - mole fraction of CO2 in air at the airblock, 4 feet above the sea surface (dry air)  
	"fCO2  SW (sat) uatm" - fugacity of CO2 in air in equilibrium with the seawater at seasurface temperature (100% humidity)  
	"fCO2  Air (sat) uatm" - fugacity of CO2 in air at the airblock, 4 feet above the seasurface (100% humidity)  
	"dfCO2" - difference of the fugacity of the CO2 in seawater and the fugacity ofthe CO2 in air (fCO2 SW - fCO2 Air)  
	"pCO2 SW (sat) uatm" - partial Pressure of CO2 in air in equilibrium with the seawater at sea surface temperature (100% humidity)  
	"pCO2 Air (sat) uatm" - partial Pressure of CO2 in air at the airblock, 4 feet above thesea surface (100% humidity)  
	"dpCO2" - difference of the partial pressure of CO2 in seawater and air (pCO2SW - pCO2 Air)  
	"pH_SW" - sea water pH
	"pH_QF" - quality flag of sea water pH
