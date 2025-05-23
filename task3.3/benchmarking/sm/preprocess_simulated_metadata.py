#%%
import pathlib as pl
import pandas as pd
import pathlib as pl
import warnings
import numpy as np
import xarray as xr
import glob 

#simulation_directory = pl.Path("/home/l.avila/DATA/Stars4Water/data/simulations")
simulation_directory = pl.Path("/p/data1/slts/avila2/Stars4Water/WP3/benchmarking/data/simulations")

output_directory = pl.Path("../saves/simulations/")

simulation_patterns = [dir.stem for dir in simulation_directory.iterdir() if dir.is_dir()] # GOEframe has its own benchmark and parflowclm_hres has too few simulation dates
simulation_fields = [pattern.split("_") for pattern in simulation_patterns]
simulation_patterns = pd.DataFrame(data = simulation_fields, index = simulation_patterns)
simulation_patterns.columns = ["model", "meteo", "region", "resolution"]

meta=pd.read_csv("../observations/meta_region.csv")

regions=meta['region'].unique()

for region in regions:
    
    print(region)

    meta = pd.read_parquet("../saves/observations/{}/meta.parquet".format(region))
    meta=meta.iloc[:,1:]

    patterns=simulation_patterns.index 

    for pattern in patterns:
        
        if pattern=='wflowsbm_era5_europe_30sec':
            continue

        pattern_directory = pl.Path("{}/{}".format(simulation_directory, pattern))

        moisture_files = np.array([file for file in pattern_directory.iterdir() if file.is_file() and file.stem.split("_")[2] == "sm"])
        moisture_files = np.sort(moisture_files)

        moisture_file = moisture_files[0]

        pattern_meta = meta.copy()

        dataset=xr.open_dataset(moisture_files[0]).sm.isel(time=0).drop_vars(["time"])
        resolution = dataset["lon"][1] - dataset["lon"][0]
        tol=0.4
        
        for index,row in meta.iterrows():   
            
            lon=row['lon']
            lat=row['lat']         
            
            if 'wflowsbm' in pattern:
                
                try:
                    nearest = dataset.where((dataset.lon<=lon+tol) &
                                        (dataset.lon>=lon-tol) &
                                        (dataset.lat<=lat+tol) &
                                        (dataset.lat>=lat-tol),drop=True)


                    distance = np.sqrt((nearest.lon - lon)**2 + (nearest.lat - lat)**2)
                    valid_distances = xr.where(~np.isnan(nearest), distance, np.inf)

                    nearest_index = np.unravel_index(np.nanargmin(valid_distances.values), nearest.shape)

                    nearest_value = nearest.values[nearest_index]
                    
                    if np.isnan(nearest_value)==False:

                        nearest_coords = {dim: nearest[dim][i].item() for dim, i in zip(nearest.dims, nearest_index)}
                        
                        pattern_meta.at[index, "simulated_lat"] = nearest_coords['lat']
                        pattern_meta.at[index, "simulated_lon"] = nearest_coords['lon']   
                except:
                    pass     
            else: 
                 try:
                
                    moisture_point = dataset.sel(lon = row["lon"],
                                                            lat = row["lat"],
                                                            method = "nearest",tolerance=resolution / 2 + 1e-10)   

                    pattern_meta.at[index, "simulated_lat"] = moisture_point.coords["lat"]
                    pattern_meta.at[index, "simulated_lon"] = moisture_point.coords["lon"]
    
                 except:
                    pass

        if 'simulated_lon' in pattern_meta:
            
            pattern_meta = pattern_meta.dropna(subset=["simulated_lat", "simulated_lon"])

            meta_out = pl.Path("../saves/simulations/{}/{}/meta.csv".format(region,pattern))
            meta_out.parent.mkdir(parents=True, exist_ok=True)
            pattern_meta.to_csv(meta_out)
            
            print('\t',pattern)
        
            
      
    


