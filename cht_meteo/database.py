import os

# import yaml
import toml

from .coamps_tc_forecast_s3 import MeteoDatasetCOAMPSTCForecastS3
from .dataset import MeteoDataset
from .gfs_anl_0p50 import MeteoDatasetGFSAnalysis0p50
from .gfs_forecast_0p25 import MeteoDatasetGFSForecast0p25
from .gfs_forecast_0p25_ncar_archive import MeteoDatasetGFSForecast0p25NCARArchive
from .ecmwf_forecast_0p25 import MeteoDatasetECMWFForecast0p25
from .matroos_forecast import MeteoDatasetMatroos


class MeteoDatabase:
    def __init__(self, path=None):
        # Initialize the database
        self.path = path
        self.dataset = {}
        # self.source = {}
        # self.set_sources()

    def print_datasets(self):
        # Print a list of all datasets
        for dataset_name, dataset in self.dataset.items():
            print(dataset_name + " - source : " + dataset.source_name)

    def list_sources(self):
        # Returns a list the available sources
        return ["gfs_forecast_0p25", "gfs_analysis_0p50", "coamps_tc_forecast", "ecmwf_forecast_0p25", "matroos", "custom"]

    def list_dataset_names(self):
        # Returns a list the available sources
        lst = []
        for dataset_name, dataset in self.dataset.items():
            lst.append(dataset_name)
        return lst

    def add_dataset(self, dataset_name, source_name, **kwargs):
        # Add a dataset to the database
        dataset_path = os.path.join(self.path, dataset_name)

        if source_name is not None:
            if source_name == "coamps_tc_forecast":
                md = MeteoDatasetCOAMPSTCForecastS3(
                    name=dataset_name, path=dataset_path, **kwargs
                )
            elif source_name == "gfs_forecast_0p25":
                md = MeteoDatasetGFSForecast0p25(
                    name=dataset_name, path=dataset_path, **kwargs
                )
            elif source_name == "gfs_forecast_0p25_ncar_archive":
                md = MeteoDatasetGFSForecast0p25NCARArchive(
                    name=dataset_name, path=dataset_path, **kwargs
                )
            elif source_name == "gfs_analysis_0p50":
                md = MeteoDatasetGFSAnalysis0p50(
                    name=dataset_name, path=dataset_path, **kwargs
                )
            elif source_name == "ecmwf_forecast_0p25":
                md = MeteoDatasetECMWFForecast0p25(
                    name=dataset_name, path=dataset_path, **kwargs
                )
            elif source_name == "matroos":
                md = MeteoDatasetMatroos(name=dataset_name, path=dataset_path, **kwargs)
            elif source_name == "custom":
                # Import class from specific python file location
                if "filepath" not in kwargs:
                    print("Error while reading meteo database : for source 'custom' the filepath to the python file must be provided")
                    return
                filepath = kwargs.pop("filepath")
                if not os.path.exists(filepath):
                    print(f"Error while reading meteo database : file {filepath} does not exist")
                    return
                import importlib.util
                spec = importlib.util.spec_from_file_location("custom_module", filepath)
                custom_module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(custom_module)
                md = custom_module.MeteoDatasetCustom(name=dataset_name, path=dataset_path, **kwargs)
            else:
                md = MeteoDataset(name=dataset_name, path=dataset_path, **kwargs)
                print(
                    f"Error while reading meteo database : source {source_name} not recognized"
                )

                # #TODO do we want to add download functionality in a flexible way to the generic MeteoDatasets? E.g. using a file with download function only.
                # if "download_forecast_cycle" in kwargs:
                #     md.download_forecast_cycle = download

        else:
            # Use generic meteo dataset (this does not have download functionality)
            md = MeteoDataset(name=dataset_name, path=dataset_path, **kwargs)

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
            print(
                f"Error while reading meteo database : file {filename} does not exist"
            )
            return

        # Read the toml file
        with open(filename) as f:
            contents = toml.load(f)

        dataset_list = contents["meteo_dataset"]
        # Loop through datasets and add them to the database
        for meteo_dataset in dataset_list:
            # pop most common variables from dict, rest is parsed through kwargs
            dataset_name = meteo_dataset.pop("name")
            source_name = meteo_dataset.pop("source")
            lon_range = meteo_dataset.pop("x_range", None)
            lat_range = meteo_dataset.pop("y_range", None)
            tau = meteo_dataset.pop("tau", 0)

            self.add_dataset(
                dataset_name=dataset_name,
                source_name=source_name,
                lon_range=lon_range,
                lat_range=lat_range,
                tau=tau,
                **meteo_dataset,
            )
