from typing import Dict, Union
from pydantic import BaseModel, Extra
import xarray as xr
import mhkit.dolfyn as dolfyn
from mhkit.dolfyn.adp import api
from tsdat import DataReader


class DnFacingADCPReader(DataReader):
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

        ds = dolfyn.read(input_key)

        # Set depth below water surface

        api.clean.set_range_offset(ds, self.parameters.depth_offset)
        api.clean.find_surface_from_P(ds, salinity=self.parameters.salinity)
        ds["depth"] = ds.h_deploy + ds["dist_bt"].mean("beam")

        # Rotate to Earth coordinates
        dolfyn.set_declination(ds, self.parameters.magnetic_declination)
        dolfyn.rotate2(ds, "earth")

        return ds
