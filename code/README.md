# The Code

All the tasks were completely done in Julia. The code folder is divided into 3 subfolders, each one corresponding to a different task. In this document there is a simple description of the content of each folder and the instructions to run the code. The first time the code is run on a machine, it will take a while to download the necessary packages. The enviroment and all the packages needed to run the code will be automatically downloaded and installed.

## Task 15
This folder contains the script to execute task 15 ( `task15.jl` ) and a "Plots" folder for the output plots. To run the code these are the instructions:
1. Make sure you are either in the root directory of the repository or in the code folder: so `PoCS_Project`.
2. From the bash terminal, run the following command: `julia code/task_15/task15.jl`. This will execute the script.

## Task 34
This folder contains the script to execute task 34 ( `task34.jl` ) and a "Plots" folder for the output plots. To run the code these are the instructions:
1. Make sure you are either in the root directory of the repository or in the code folder: so `PoCS_Project`.
2. From the bash terminal, run the following command: `julia code/task_34/task34.jl`. This will execute the script.

## Task 44
This folder contains the script to execute task 44 ( `task44_1.jl` and `task44_2.jl` ), the "Plots" folder for the output plots, the "Output" folder for the output data (Edges lists and Node list) and the "Data" folder where the Data needed for the task is stored. As per instructions on the course moodle, the data is not included in the repository. To run the code additional files are needed, all the necessary information is contained in the final report.
The file needed are:
- `gadm1_nuts3_counties_levels.csv`: Provided by the course link.
- `gadm1_nuts3_counties-gadm1_nuts3_counties - FB Social Connectedness Index - October 2021.tsv`: Provided by the course link.
- `NUTS_RG_60M_2016_4326.geojson`: Obtained on the [Eurostat website](https://ec.europa.eu/eurostat/web/gisco/geodata/statistical-units/territorial-units-statistics). In the options, select the NUTS 2016 year, the GeoJSON file forma, the Polygon geometry type, the 60M resolution and the EPSG 4326 coordinate reference system.
- The "GADM_DATA" folder contains several files needed for the GADM data format analysis. The files can be obtained by clicking this link [GADM website](https://biogeo.ucdavis.edu/data/gadm2.8/gadm28.shp.zip). The link will download a zip file that will need decompression and renaming to "GADM_DATA".
- `graph_data.csv`: This file is provided with the repo and contains the pre-computed analysis of the graphs.


Given that all data files are present in the "Data" folder, to run the code these are the instructions:
1. Make sure you are either in the root directory of the repository or in the code folder: so `PoCS_Project`.
2. In the bash terminal, run the following command: `julia code/task_44/task44_1.jl`. This will execute the first part of the script where the raw data is processed and the output files are generated.
3. In the bash terminal, run the following command: `julia code/task_44/task44_2.jl`. This will execute the second part of the script where the analysis is performed and the output plots are generated.
