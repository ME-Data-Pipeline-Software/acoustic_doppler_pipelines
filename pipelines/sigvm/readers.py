from typing import Dict, Union
from pydantic import BaseModel, Extra
import xarray as xr
import os
from io import BytesIO

from tsdat import DataReader
import dolfyn
from dolfyn.adp import api


class DownFacingADCPReader(DataReader):
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

        depth_offset: float = 0.5
        salinity: float = 35
        magnetic_declination: float = 0
        correlation_threshold: float = 30

    parameters: Parameters = Parameters()
    """Extra parameters that can be set via the retrieval configuration file. If you opt
    to not use any configuration parameters then please remove the code above."""

    def read(self, input_key: str) -> Union[xr.Dataset, Dict[str, xr.Dataset]]:
        """-------------------------------------------------------------------
        SigVM datafiles are a zip folder containing two files, a .anpp file
        and a .ad2cp file. This reader skips the first .anpp file and reads
        the raw data from the .ad2cp file.

        Args:
            filename (str): The path to the ADCP file to read in.

        Returns:
            xr.Dataset: An xr.Dataset object
        -------------------------------------------------------------------"""

        # dolfyn requires a file, not a bytesIO object
        buffer: BytesIO = input_key  # type: ignore
        with open("data.ad2cp", "wb") as f:
            f.write(buffer.getvalue())

        ds = dolfyn.read("data.ad2cp")

        # Remove extra files
        os.remove("data.ad2cp")
        os.remove("data.ad2cp.index")

        api.clean.set_range_offset(ds, self.parameters.depth_offset)
        api.clean.find_surface_from_P(ds, salinity=self.parameters.salinity)

        # Locate surface using pressure data and remove data above it
        ds["depth"] = ds.h_deploy + ds["alt_dist"]
        ds = api.clean.nan_beyond_surface(ds)

        # Rotate to Earth coordinates
        dolfyn.set_declination(ds, self.parameters.magnetic_declination)
        dolfyn.rotate2(ds, "earth")

        ds = api.clean.correlation_filter(
            ds, thresh=self.parameters.correlation_threshold
        )

        # Dropping the detailed configuration stats because netcdf can't save it
        for key in list(ds.attrs.keys()):
            if "config" in key:
                ds.attrs.pop(key)
        ds.attrs.pop("rotate_vars")

        return ds
