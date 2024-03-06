import numpy as np
import pandas as pd
import xarray as xr
from scipy.signal import convolve2d, medfilt
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mhkit import dolfyn
from mhkit.dolfyn.adp import api
from mhkit.dolfyn.time import dt642epoch
from tsdat import IngestPipeline, FileSystem

import pipelines.sigvm.gis as gis
from shared.writers import MatlabWriter


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


class SigVM(IngestPipeline):
    """---------------------------------------------------------------------------------
    This is an example ingestion pipeline meant to demonstrate how one might set up a
    pipeline using this template repository.

    ---------------------------------------------------------------------------------"""

    def hook_customize_dataset(self, dataset: xr.Dataset) -> xr.Dataset:
        # (Optional) Use this hook to modify the dataset before qc is applied

        # Set depth to be the distance to the seafloor, and remove data beyond it
        dist = None
        # Try using the altimeter data
        if getattr(dataset, "use_alt_depth", 0):
            dist = dataset["le_dist_alt"]
        # Or the bottom track data
        elif getattr(dataset, "use_bt_depth", 0):
            dist_bt = dataset["dist_bt"].mean("beam")
            dist = dist_bt.interp(
                {"time_bt": dataset["time"]}, kwargs=dict(fill_value="extrapolate")
            )
        # And set depth (otherwise don't)
        if dist is not None:
            dataset["depth"] = dataset.h_deploy + dist
            api.clean.nan_beyond_surface(dataset, inplace=True)

        # Calculate GPS velocity
        if "speed_over_grnd_gps" in dataset:
            vel_E = dataset["speed_over_grnd_gps"] * np.sin(
                np.deg2rad(dataset["dir_over_grnd_gps"])
            )
            vel_N = dataset["speed_over_grnd_gps"] * np.cos(
                np.deg2rad(dataset["dir_over_grnd_gps"])
            )
            dataset["vel_gps"].loc[{"earth": "E"}] = vel_E
            dataset["vel_gps"].loc[{"earth": "N"}] = vel_N
            dataset["vel_gps"].loc[{"earth": "U"}] = vel_N * 0
        elif "latitude_gps" in dataset:
            dataset["vel_gps"] = process_gps_data(dataset)

        # Correct velocity with bottom track
        if getattr(dataset, "vel_bt_correction", 0):
            # Subtract bottom track from water velocity
            vel_corrected = dataset["vel"] - dataset["vel_bt"].interp(
                {"time_bt": dataset["time"]}, kwargs=dict(fill_value="extrapolate")
            )
            dataset["vel"].values = vel_corrected.values

        # Correct velocity with GPS
        elif getattr(dataset, "vel_gps_correction", 0):
            # Must ADD GPS velocity to remove from water velocity
            tmp = np.zeros((4, dataset["time"].size), dtype=np.float32)
            tmp[:3] = dataset["vel_gps"].interp(time_gps=dataset["time"]).values
            dataset["vel"].values += tmp[:, None]

        # Speed and Direction
        dataset["U_mag"].values = dataset.velds.U_mag
        dataset["U_dir"].values = dataset.velds.U_dir

        return dataset

    def hook_finalize_dataset(self, dataset: xr.Dataset) -> xr.Dataset:
        # (Optional) Use this hook to modify the dataset after qc is applied
        # but before it gets saved to the storage area

        # Save mat file
        storage = FileSystem()
        storage.handler.writer = MatlabWriter()

        storage.save_data(dataset)

        return dataset

    def hook_plot_dataset(self, ds: xr.Dataset):
        def add_colorbar(ax, plot, label):
            cb = plt.colorbar(plot, ax=ax, pad=0.01)
            cb.ax.set_ylabel(label, fontsize=12)
            cb.outline.set_linewidth(1)
            cb.ax.tick_params(size=0)
            cb.ax.minorticks_off()
            return cb

        date = pd.to_datetime(ds["time"].values)

        y_max = int(ds["depth"].mean() * 1.5)

        with plt.style.context("shared/styling.mplstyle"):
            # Current
            fig, ax = plt.subplots(
                nrows=2, ncols=1, figsize=(14, 8), constrained_layout=True
            )

            velE = ax[0].pcolormesh(
                date,
                -ds["range"],
                ds["vel"][0],
                cmap="coolwarm",
                shading="nearest",
            )
            ax[0].plot(date, -ds["depth"])
            ax[0].set_xlabel("Time (UTC)")
            ax[0].set_ylabel(r"Range [m]")
            ax[0].set_ylim([-y_max, 0])
            add_colorbar(ax[0], velE, r"Velocity East [m/s]")
            velE.set_clim(-3, 3)

            velN = ax[1].pcolormesh(
                date,
                -ds["range"],
                ds["vel"][1],
                cmap="coolwarm",
                shading="nearest",
            )
            ax[1].plot(date, -ds["depth"])
            ax[1].set_xlabel("Time (UTC)")
            ax[1].set_ylabel(r"Range [m]")
            ax[1].set_ylim([-y_max, 0])
            add_colorbar(ax[1], velN, r"Velocity North [m/s]")
            velN.set_clim(-3, 3)

            plot_file = self.get_ancillary_filepath(title="current")
            fig.savefig(plot_file)
            plt.close(fig)

            # Amplitude
            fig, ax = plt.subplots(
                nrows=ds.n_beams, ncols=1, figsize=(14, 8), constrained_layout=True
            )

            for beam in range(ds.n_beams):
                amp = ax[beam].pcolormesh(
                    date, -ds["range"], ds["amp"][beam], shading="nearest"
                )
                ax[beam].set_title("Beam " + str(beam + 1))
                ax[beam].set_xlabel("Time (UTC)")
                ax[beam].set_ylabel(r"Range [m]")
                ax[beam].set_ylim([-y_max, 0])
                add_colorbar(ax[beam], amp, "Amplitude [dB]")

            plot_file = self.get_ancillary_filepath(title="amplitude")
            fig.savefig(plot_file)
            plt.close(fig)

            # Correlation
            fig, ax = plt.subplots(
                nrows=ds.n_beams, ncols=1, figsize=(14, 8), constrained_layout=True
            )

            for beam in range(ds.n_beams):
                corr = ax[beam].pcolormesh(
                    date,
                    -ds["range"],
                    ds["corr"][beam],
                    cmap="copper",
                    shading="nearest",
                )
                ax[beam].set_title("Beam " + str(beam + 1))
                ax[beam].set_xlabel("Time (UTC)")
                ax[beam].set_ylabel(r"Range [m]")
                ax[beam].set_ylim([-y_max, 0])
                add_colorbar(ax[beam], corr, "Correlation [%]")

            plot_file = self.get_ancillary_filepath(title="correlation")
            fig.savefig(plot_file)
            plt.close(fig)

            # Lat/lon
            fig, ax = plt.subplots()

            h = ax.scatter(
                ds["longitude_gps"],
                ds["latitude_gps"],
                c=ds["U_mag"].mean("range").interp(time=ds["time_gps"]).values,
                cmap="Blues",
                s=100,
            )
            fig.colorbar(h, ax=ax, label="Current Speed [m/s]")
            ax.quiver(
                ds["longitude_gps"],
                ds["latitude_gps"],
                ds["vel"][0].mean("range").interp(time=ds["time_gps"]).values,
                ds["vel"][1].mean("range").interp(time=ds["time_gps"]).values,
            )

            ax.set_title("")  # Remove bogus title created by xarray
            ax.set_ylabel("Latitude [deg N]")
            ax.set_xlabel("Longitude [deg E]")
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.4f"))
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.4f"))

            plot_file = self.get_ancillary_filepath(title="location")
            fig.savefig(plot_file)
            plt.close(fig)
