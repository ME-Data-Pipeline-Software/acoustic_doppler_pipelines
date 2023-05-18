from typing import Dict, Union
from pydantic import BaseModel, Extra
import numpy as np
import xarray as xr
import pandas as pd
import os
from io import BytesIO
import pynmea2

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

        ds = dolfyn.io.nortek2.read_signature("data.ad2cp", rebuild_index=True)

        # Remove extra files
        os.remove("data.ad2cp")
        os.remove("data.ad2cp.index")

        # Set depth below water surface
        api.clean.set_range_offset(ds, self.parameters.depth_offset)
        api.clean.find_surface_from_P(ds, salinity=self.parameters.salinity)

        # Rotate to Earth coordinates
        dolfyn.set_declination(ds, self.parameters.magnetic_declination)
        dolfyn.rotate2(ds, "earth")

        return ds


class NMEAReader(DataReader):
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

        buffer: BytesIO = input_key  # type: ignore
        with open("data.nmea", "wb") as f:
            f.write(buffer.getvalue())

        df = pd.read_csv("data.nmea", header=None, sep="\\")
        d = {}
        for index, row in df.iterrows():
            try:
                time = float(row[1][2:16])
                d[time] = pynmea2.parse(row[2])
            except:
                pass

        new_dict: dict = {
            "gga_time": [],
            "gga_lat": [],
            "gga_lon": [],
            "gga_gps_qual": [],
            "gga_num_sats": [],
            "gga_alt": [],
            "gga_hdop": [],
            "rmc_time": [],
            "vtg_time": [],
            "vtg_course_deg": [],
            "vtg_spd_over_grnd_kt": [],
            "hdt_time": [],
            "hdt_heading": [],
        }

        for t in d:
            try:  # skip strange sentences
                d[t].sentence_type
            except:
                continue

            if "GGA" in d[t].sentence_type:
                new_dict["gga_time"].append(float(t))
                new_dict["gga_lat"].append(d[t].latitude)
                new_dict["gga_lon"].append(d[t].longitude)
                new_dict["gga_gps_qual"].append(int(d[t].gps_qual))
                new_dict["gga_num_sats"].append(int(d[t].num_sats))
                new_dict["gga_alt"].append(float(d[t].altitude))
                new_dict["gga_hdop"].append(float(d[t].horizontal_dil))

            elif "RMC" in d[t].sentence_type:
                new_dict["rmc_time"].append(float(t))

            elif "VTG" in d[t].sentence_type:
                new_dict["vtg_time"].append(float(t))
                new_dict["vtg_course_deg"].append(float(d[t].true_track))
                new_dict["vtg_spd_over_grnd_kt"].append(float(d[t].spd_over_grnd_kts))

            elif "HDT" in d[t].sentence_type:
                new_dict["hdt_time"].append(float(t))
                new_dict["hdt_heading"].append(float(d[t].heading))

        os.remove("data.nmea")

        for ky in new_dict:
            tg = ky[:3]
            if "time" in ky:
                new_dict[ky] = dolfyn.time.epoch2dt64(np.array(new_dict[ky]))
            new_dict[ky] = {"dims": (tg + "_time"), "data": new_dict[ky]}

        ds = xr.Dataset.from_dict(new_dict)
        ds = ds.interp(
            gga_time=ds.rmc_time,
            vtg_time=ds.rmc_time,
            hdt_time=ds.rmc_time,
            method="nearest",
        ).rename({"rmc_time": "time"})

        return ds
