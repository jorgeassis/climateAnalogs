The `model` folder contains the model data from Li-Qing Jiang that was used to do the analysis. 
Each data point is the decadal mean for an ocean gridpoint at that month.
	* `Temp_Arag_1800_2000.txt` Model data between 1800 and 2000
	* `Temp_Arag_2070_2100_RCP85.txt` Model data between 2070 and 2100 for RCP 8.5
	* `Temp_Arag_2070_2100_RCP45.txt` Model data between 2070 and 2100 for RCP 4.5
	* `ICV_Temp_Arag_pH_1960_2020.txt` ICV of model data used for each gridpoint
	
These four files have the same headers:
* `No` Station number 
* `Lon`Longitude
* `Lat` Latitude
* `Year` Year
* `Month` Month
* `SST` Sea surface temperature (C)
* `Arag` Aragonite saturation state
* `pH` pH

The `real_data` folder contains the empirical data used for comparison to model data.
See README within that folder.