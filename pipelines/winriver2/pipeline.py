import warnings
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mhkit import dolfyn
from mhkit.dolfyn.adp import api
from tsdat import IngestPipeline

from shared.gps import process_gps_data


class WinRiver2ADCP(IngestPipeline):
    """---------------------------------------------------------------------------------
    This pipeline runs motion correction of ADCP measurements collected from a vessel
    using the Teledyne WinRiver2 software.
    ---------------------------------------------------------------------------------"""

    def hook_customize_dataset(self, dataset: xr.Dataset) -> xr.Dataset:
        # (Optional) Use this hook to modify the dataset before qc is applied

        # Don't bother if less than a minute of data recorded
        if dataset["time"].size <= (60 * dataset.fs):
            raise AssertionError(
                "File is less than a minute long. Not running through pipeline."
            )

        # Remove GPS coords if a file doesn't exist
        if "time" in dataset["time_gps"].dims:
            dataset["time_gps"] = (
                dataset["time_gps"].swap_dims(time="time_gps").drop_vars("time")
            )

        ## Set depth to be the distance to the seafloor, and remove data beyond it
        dist = None
        # Try using the altimeter data
        if getattr(dataset, "use_rdi_depth"):
            dist = dataset["depth"]
        # Or the bottom track data
        elif getattr(dataset, "use_bt_depth"):
            dist_bt = (
                dataset["dist_bt"].min("beam") * 0.95
            )  # need to clip a little further
            dist = dist_bt
        # And set depth (otherwise don't)
        if dist is not None:
            dist = dist.where(dist > 0, np.nan)  # Remove less than 0
            dataset["depth"].values = dataset.attrs["range_offset"] + dist
            api.clean.remove_surface_interference(dataset, inplace=True)

        # Realign ADCP based on GPS heading
        if (
            "heading_gps" in dataset
            and not all(dataset["heading_gps"] == dataset["heading_gps"]._FillValue)
            and len(dataset["heading_gps"]) >= (0.5 * len(dataset["time"]))
            and (not dataset.coord_sys == "ship")
        ):
            # Remove magnetic declination and set heading based on GPS
            dataset.attrs["rotate_vars"] = ["vel", "vel_bt"]
            dolfyn.set_declination(dataset, -dataset.declination)
            dolfyn.rotate2(dataset, "inst")

            # Manually realign beam 3 (y-axis in Nortek coordinate system) to GPS heading
            offset = dataset.attrs["heading_misalign_deg"]
            warnings.warn(
                f"Using heading misalignment of {offset}. Update dn_looking_adcp/pipeline.py as necessary."
                "FYI: Heading misaligment should be 45 degrees if TRDI ADCP Y-axis (beam 3 notch) is"
                "rotated 45 degrees to port. It should be -45 if rotated 45 degrees to starboard."
            )
            # Change sign or value if necessary
            dataset["heading"] = (
                (dataset["heading_gps"] + dataset.attrs["heading_misalign_deg"]) % 360
            ).interp(time_gps=dataset["time"])

            # Recreate orientation matrix
            dataset = dataset.drop_vars(["orientmat"])
            dataset["orientmat"] = dolfyn.rotate.vector._euler2orient(
                dataset["time"], dataset["heading"], dataset["pitch"], dataset["roll"]
            )
            # Rotate to earth now in true coordinates
            dolfyn.rotate2(dataset, "earth")

        elif dataset.coord_sys == "ship":
            # If WinRiver saved in "ship" coordinates, all of the above is already applied
            dolfyn.rotate2(dataset, "earth")

        ## Motion Correction
        # Correct velocity with bottom track
        if getattr(dataset, "vel_bt_correction") and not getattr(
            dataset, "vel_gps_correction"
        ):
            # Subtract bottom track from water velocity
            vel_corrected = dataset["vel"] - dataset["vel_bt"]
            dataset["vel"].values = vel_corrected.values

        # Correct velocity with GPS or with both GPS and BT
        elif getattr(dataset, "vel_gps_correction"):
            # If GPS file is missing (time_gps is set as time)
            if (dataset["time"] == dataset["time_gps"]).all():
                if getattr(dataset, "vel_bt_correction"):
                    warnings.warn(
                        "No GPS recorded. Velocity correction completed with BT alone."
                    )
                    vel_gps = dataset["vel"][:, 0, :] * np.nan
                    dataset.attrs["vel_gps_correction"] = 0
                else:
                    raise Exception("No GPS data available.")
            # Otherwise calculate GPS velocity
            else:
                # Try to get VTG reading first
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
                # If not try to use the GGA reading
                else:
                    dataset["vel_gps"] = process_gps_data(dataset)

                # Must ADD GPS velocity to remove from water velocity
                vel_gps = np.zeros((4, dataset["time"].size), dtype=np.float32)
                vel_gps[:3] = dataset["vel_gps"].interp(time_gps=dataset["time"]).values

            # If using both GPS and BT, fill missing BT with GPS timestamps
            if getattr(dataset, "vel_bt_correction") and dist is not None:
                vel_bt = dataset["vel_bt"]
                # Double negative (negate vel_gps to put in same coordinate system as vel_bt,
                # negate again to add to velocity to stay consistent with this elif block)
                vel_cor = -vel_bt.where(
                    (
                        (abs(dataset["depth"]) < 9999)
                        & (abs(vel_bt) <= vel_bt.valid_max)
                    ),
                    -vel_gps,
                )
            elif dist is None:
                raise Exception(
                    "Cannot use 'vel_bt_correction' if 'depth' is unknown. Please set "
                    "'use_alt_depth' or 'use_bt_depth' to 1."
                )
            else:
                vel_cor = vel_gps[:, None]

            dataset["vel"].values = (dataset["vel"] + vel_cor).values

        # Speed and Direction
        dataset["U_mag"] = dataset.velds.U_mag
        dataset["U_dir"] = dataset.velds.U_dir

        return dataset

    def hook_finalize_dataset(self, dataset: xr.Dataset) -> xr.Dataset:
        # (Optional) Use this hook to modify the dataset after qc is applied
        # but before it gets saved to the storage area

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
            velE.set_clim(ds["vel"].attrs["valid_min"], ds["vel"].attrs["valid_max"])

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
            velN.set_clim(ds["vel"].attrs["valid_min"], ds["vel"].attrs["valid_max"])

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
            h.set_clim(0, ds["vel"].attrs["valid_max"])
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
