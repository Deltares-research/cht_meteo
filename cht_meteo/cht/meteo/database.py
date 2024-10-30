import os
# import yaml
import toml

from .dataset import MeteoDataset
from .coamps_tc_forecast_s3 import MeteoDatasetCOAMPSTCForecastS3
from .gfs_forecast_0p25_v3 import MeteoDatasetGFSForecast0p25
from .gfs_anl_0p50_05 import MeteoDatasetGFSAnalysis0p50

class MeteoDatabase():
    def __init__(self):
        self.path = None
        self.dataset = {}
        # self.source = {}
        # self.set_sources()

    def print_datasets(self):
        # Print a list of all datasets
        for dataset_name, dataset in self.dataset.items():
            print(dataset_name + " - source : " + dataset.source_name)

    def list_sources(self):
        # Returns a list the available sources
        return ["gfs_forecast", "coamps_tc_forecast"]    

    def add_dataset(self, dataset_name, source_name, **kwargs):

        # Add a dataset to the database
        dataset_path = os.path.join(self.path, dataset_name)

        if source_name is not None:
            if source_name == "coamps_tc_forecast":
                md =  MeteoDatasetCOAMPSTCForecastS3(name=dataset_name, path=dataset_path, **kwargs)
            elif source_name == "gfs_forecast_0p25":
                md = MeteoDatasetGFSForecast0p25(name=dataset_name, path=dataset_path, **kwargs)
            elif source_name == "gfs_analysis_0p50":
                md = MeteoDatasetGFSAnalysis0p50(name=dataset_name, path=dataset_path, **kwargs)
            else:
                md = MeteoDataset(name=dataset_name)
                print(f"Error while reading meteo database : source {source_name} not recognized")

        else:
            # Use generic meteo dataset (this does not have download functionality)
            md = MeteoDataset(name=dataset_name)

        # Add to database        
        self.dataset[dataset_name] = md

        return md

    def read_datasets(self, filename=None):
        """Read the datasets from a toml file"""

        if filename is None:
            filename = os.path.join(self.path, "meteo_database.toml")
        else:
            self.path = os.path.dirname(filename)

        # Check if the file exists
        if not os.path.exists(filename):
            print(f"Error while reading meteo database : file {filename} does not exist")
            return    

        # Read the toml file
        with open(filename) as f:
            contents = toml.load(f)

        dataset_list = contents["meteo_dataset"]
        # Loop through datasets and add them to the database
        for meteo_dataset in dataset_list:
            if "x_range" in meteo_dataset:
                lon_range = meteo_dataset["x_range"]
            else:
                lon_range = None
            if "y_range" in meteo_dataset:
                lat_range = meteo_dataset["y_range"]
            else:
                lat_range = None
            if "tau" in meteo_dataset:
                tau = meteo_dataset["tau"]
            else:
                tau = 0

            self.add_dataset(meteo_dataset["name"],
                            source_name=meteo_dataset["source"],
                            lon_range=lon_range,
                            lat_range=lat_range,
                            tau=tau)
