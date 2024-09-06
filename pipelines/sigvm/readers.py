from typing import Dict, Union
from pydantic import BaseModel, Extra
import numpy as np
import xarray as xr
import pandas as pd
import os
import warnings
from io import BytesIO
import pynmea2
import mhkit.dolfyn as dolfyn
from mhkit.dolfyn.adp import api
from tsdat import DataReader
from tsdat.qc.checkers import CheckMonotonic


class SigVMReader(DataReader):
    """---------------------------------------------------------------------------------
    Custom DataReader that can be used to read data from a specific format.

    Built-in implementations of data readers can be found in the
    [tsdat.io.readers](https://tsdat.readthedocs.io/en/latest/autoapi/tsdat/io/readers)
    module.

    ---------------------------------------------------------------------------------"""

    class Parameters(BaseModel, extra=Extra.forbid):
        depth_offset: float = 0.5
        magnetic_declination: float = 0

    parameters: Parameters = Parameters()

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

        # Trim certain SigVM files with extra header lines - this will likely change
        idx = buffer.getvalue().find(b"GETCLOCKSTR") - 11
        if idx:
            buffer_data = buffer.getvalue()[idx:]  # Should get firmware version
        else:
            buffer_data = buffer.getvalue()

        # Save temp file to read with dolfyn
        with open("data.ad2cp", "wb") as f:
            f.write(buffer_data)

        ds = dolfyn.io.nortek2.read_signature("data.ad2cp", rebuild_index=True)

        # Remove temp files
        os.remove("data.ad2cp")
        os.remove("data.ad2cp.index")

        # Set range given ADCP depth
        api.clean.set_range_offset(ds, self.parameters.depth_offset)

        # Rotate to Earth coordinates
        dolfyn.set_declination(ds, self.parameters.magnetic_declination)
        dolfyn.rotate2(ds, "earth")

        # Kill repeated timestamps
        time_stamps = [c for c in ds.coords if "time" in c]
        checker = CheckMonotonic()
        bad = checker.run(ds, "time")
        idx = np.argwhere(~bad).squeeze()
        for t in time_stamps:
            ds = ds.isel({t: idx})

        # Add time_gps since tsdat can't add extra coordinates
        ds = ds.assign_coords({"time_gps": ds["time"]})

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

        def get(item, key):
            out = getattr(item, key)
            if out is None:
                out = np.nan
            return out

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
                new_dict["gga_lat"].append(get(d[t], "latitude"))
                new_dict["gga_lon"].append(get(d[t], "longitude"))
                new_dict["gga_gps_qual"].append(int(get(d[t], "gps_qual")))
                new_dict["gga_num_sats"].append(int(get(d[t], "num_sats")))
                new_dict["gga_alt"].append(float(get(d[t], "altitude")))
                new_dict["gga_hdop"].append(float(get(d[t], "horizontal_dil")))

            elif "RMC" in d[t].sentence_type:
                new_dict["rmc_time"].append(float(t))

            elif "VTG" in d[t].sentence_type:
                new_dict["vtg_time"].append(float(t))
                new_dict["vtg_course_deg"].append(float(get(d[t], "true_track")))
                new_dict["vtg_spd_over_grnd_kt"].append(
                    float(get(d[t], "spd_over_grnd_kts"))
                )

            elif "HDT" in d[t].sentence_type:
                new_dict["hdt_time"].append(float(t))
                new_dict["hdt_heading"].append(float(get(d[t], "heading")))

        os.remove("data.nmea")

        for ky in new_dict:
            tg = ky[:3]
            if "time" in ky:
                new_dict[ky] = dolfyn.time.epoch2dt64(np.array(new_dict[ky]))
            new_dict[ky] = {"dims": (tg + "_time"), "data": new_dict[ky]}

        ds = xr.Dataset.from_dict(new_dict)
        if any(ds["rmc_time"]):
            ds = ds.interp(
                gga_time=ds["rmc_time"],
                vtg_time=ds["rmc_time"],
                hdt_time=ds["rmc_time"],
                method="nearest",
            ).rename({"rmc_time": "time_gps"})
        else:
            try:  # I hate pandas
                warnings.warn("RMC timestamp not detected. Using GGA timestamp.")
                ds = ds.interp(
                    vtg_time=ds["gga_time"],
                    hdt_time=ds["gga_time"],
                    method="nearest",
                ).rename({"gga_time": "time_gps"})
            except:
                try:
                    warnings.warn("GGA timestamp failed. Using VTG timestamp.")
                    ds = ds.interp(
                        gga_time=ds["vtg_time"],
                        hdt_time=ds["vtg_time"],
                        method="nearest",
                    ).rename({"vtg_time": "time_gps"})
                except:
                    warnings.warn("VTG timestamp failed. Using HDT timestamp.")
                    ds = ds.interp(
                        gga_time=ds["hdt_time"],
                        vtg_time=ds["hdt_time"],
                        method="nearest",
                    ).rename({"hdt_time": "time_gps"})

        return ds
