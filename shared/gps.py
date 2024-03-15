import numpy as np
import xarray as xr
from scipy.signal import convolve2d
from mhkit.dolfyn.time import dt642epoch
import shared.gis as gis


def center_diff(x, axis=0):
    """Centered differences, along axis.
    except at edges where first-differences are used."""
    if axis != 0:
        x = np.moveaxis(x, axis, 0)
    out = np.zeros_like(x)
    out[1:-1] = (x[2:] - x[:-2]) / 2
    out[0] = x[1] - x[0]
    out[-1] = x[-1] - x[-2]
    if axis != 0:
        out = np.moveaxis(out, 0, axis)
    return out


def calc_gps_vel(lonlat, dt_sec):
    prj0 = gis.simple_proj(lonlat[:, 0])
    xy = prj0(lonlat)
    vel = center_diff(xy, axis=1)
    vel /= dt_sec
    # add 0's for vertical velocity
    vel = np.append(vel, np.zeros((1, vel.shape[1])), axis=0)
    return vel


def process_gps_data(ds):
    lonlat = xr.DataArray(
        np.stack((ds.longitude_gps, ds.latitude_gps)),
        dims=["gps", "time_gps"],
    )
    if lonlat[0, 0] == 0:
        # firstgd = np.nonzero()
        lonlat[:, [0]] = lonlat[:, [1]] - np.diff(lonlat[:, 1:3], axis=1)

    # Difference between GPS sampling frequency and velocity sampling frequency
    fs = getattr(
        ds,
        "fs",
        1 / np.round(np.mean(np.diff(dt642epoch(ds["time"]))), 4),
    )
    dt_diff = np.round(np.mean(np.diff(dt642epoch(ds["time_gps"]))), 4) / (1 / fs)

    # Calculate GPS velocity
    vel_gps = xr.DataArray(calc_gps_vel(lonlat, dt_diff), dims=["earth", "time_gps"])
    window = np.hanning(7)[None, :]
    vel_gps.values = convolve2d(vel_gps.values, window, "same")

    return vel_gps.fillna(0).astype("float32")
