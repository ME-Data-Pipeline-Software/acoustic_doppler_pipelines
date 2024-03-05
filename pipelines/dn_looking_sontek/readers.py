from typing import Dict, Union
from pydantic import BaseModel, Extra
import numpy as np
import pandas as pd
import xarray as xr
from tsdat import DataReader


class ReadSontek(DataReader):
    """---------------------------------------------------------------------------------
    Custom DataReader that can be used to read data from a specific format.

    Built-in implementations of data readers can be found in the
    [tsdat.io.readers](https://tsdat.readthedocs.io/en/latest/autoapi/tsdat/io/readers)
    module.

    ---------------------------------------------------------------------------------"""

    class Parameters(BaseModel, extra=Extra.forbid):
        """If your CustomDataReader should take any additional arguments from the
        retriever configuration file, then those should be specified here.

        e.g.,:
        custom_parameter: float = 5.0

        """

    parameters: Parameters = Parameters()
    """Extra parameters that can be set via the retrieval configuration file. If you opt
    to not use any configuration parameters then please remove the code above."""

    def read(self, input_key: str) -> Union[xr.Dataset, Dict[str, xr.Dataset]]:
        vel_file = input_key
        summary_file = input_key[:-4] + ".sum"

        df = pd.read_csv(
            vel_file,
            header=0,
            index_col="Date/Time",
            sep=",",
            infer_datetime_format=True,
        )
        dataset = df.to_xarray()

        dfsum = pd.read_csv(
            summary_file,
            header=0,
            index_col="Date/Time",
            sep=",",
            infer_datetime_format=True,
        )
        dataset = xr.merge((dataset, dfsum.to_xarray()))

        # Compress row of variables in input into variables dimensioned by time and height
        raw_categories = [
            # "Sample #",
            # "Date/Time",
            # "Frequency (MHz)",
            # "Profile Type",
            # "Depth (m)",
            # "Cell Size (m)",
            # "Cell Start (m)",
            "Location (m)",
            "Ve (m/s)",
            "Vn (m/s)",
            "Vu (m/s)",
            "Vd (m/s)",
            "Spd (m/s)",
            "Dir (deg)",
        ]
        output_var_names = [
            # "number",
            # "time",
            # "frequency",
            # "profile_type",
            # "depth",
            # "cell_size",
            # "blank_dist",
            "bin_depth",
            "vel_east",
            "vel_north",
            "vel_up",
            "vel_err",
            "U_mag",
            "U_dir",
        ]

        n_cells = len([v for v in dataset.data_vars if "Location" in v])
        for category, output_name in zip(raw_categories, output_var_names):
            var_names = [f"Cell{i} {category}" for i in range(1, n_cells + 1)]
            var_data = [dataset[name].data for name in var_names]
            dataset = dataset.drop_vars(var_names)
            var_data = np.array(var_data)  # type: ignore
            dataset[output_name] = xr.DataArray(var_data, dims=["range", "time"])

        return dataset
