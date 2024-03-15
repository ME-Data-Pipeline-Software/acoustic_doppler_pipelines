import warnings
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mhkit import dolfyn
from mhkit.dolfyn.adp import api
from tsdat import IngestPipeline, FileSystem

from shared.gps import process_gps_data
from shared.writers import MatlabWriter


class DnLookingADCP(IngestPipeline):
    """---------------------------------------------------------------------------------
    This is an example ingestion pipeline meant to demonstrate how one might set up a
    pipeline using this template repository.

    ---------------------------------------------------------------------------------"""

    def hook_customize_dataset(self, dataset: xr.Dataset) -> xr.Dataset:
        # (Optional) Use this hook to modify the dataset before qc is applied

        # Set depth to be the distance to the seafloor, and remove data beyond it
        api.clean.nan_beyond_surface(dataset, inplace=True)

        # Realign ADCP based on GPS heading
        if "heading_gps" in dataset and not all(
            dataset["heading_gps"] == dataset["heading_gps"]._FillValue
        ):
            # Remove magnetic declination and set heading based on GPS
            dataset.attrs["rotate_vars"] = ["vel", "vel_bt"]
            dolfyn.set_declination(dataset, -dataset.declination)
            dolfyn.rotate2(dataset, "inst")

            # Manually realign beam 3 (y-axis in Nortek coordinate system) to GPS heading
            if "rdi" in dataset.inst_make.lower():
                warnings.warn(
                    "Assumed TRDI ADCP Y-axis rotated -45 degrees (to port). Switch sign in "
                    "dn_looking_adcp/pipeline.py to (+) if rotated to starboard."
                )
            else:
                warnings.warn(
                    "Assumed Nortek ADCP X-axis rotated -45 degrees (to starboard). Switch sign in "
                    "dn_looking_adcp/pipeline.py to (+) if rotated to port."
                )
            dataset.attrs["heading_misalign_deg"] = -45
            dataset["heading"] = ((dataset["heading_gps"] - 45) % 360).interp(
                time_gps=dataset["time"]
            )

            dataset = dataset.drop_vars(["orientmat"])
            dataset["orientmat"] = dolfyn.rotate.vector._euler2orient(
                dataset["time"], dataset["heading"], dataset["pitch"], dataset["roll"]
            )
            # Rotate to earth now in true coordinates
            dolfyn.rotate2(dataset, "earth")

        # Correct velocity with bottom track
        if getattr(dataset, "vel_bt_correction", 0):
            # Subtract bottom track from water velocity
            if "time_bt" in dataset.coords:
                vel_corrected = dataset["vel"] - dataset["time_bt"].interp(
                    {"time_bt": dataset["time"]}, kwargs=dict(fill_value="extrapolate")
                )
            else:
                vel_corrected = dataset["vel"] - dataset["vel_bt"]

            dataset["vel"].values = vel_corrected.values

        # Correct velocity with GPS
        elif getattr(dataset, "vel_gps_correction", 0):
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

        y_max = ds["depth"].max() * 1.1

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
                add_colorbar(ax[beam], amp, "Amplitude")

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
                add_colorbar(ax[beam], corr, "Correlation")

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
