import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from tsdat import IngestPipeline, FileSystem

from shared.writers import MatlabWriter


class DnFacingSontek(IngestPipeline):
    """---------------------------------------------------------------------------------
    This is an example ingestion pipeline meant to demonstrate how one might set up a
    pipeline using this template repository.

    ---------------------------------------------------------------------------------"""

    def hook_customize_dataset(self, dataset: xr.Dataset) -> xr.Dataset:
        # (Optional) Use this hook to modify the dataset before qc is applied

        ## Remap variable cell sizes onto new range grid
        max_bd = np.round(dataset["blank_dist"].max(), 3)  # 0.35
        max_depth = np.round(dataset["depth"].max(), 3)  # 5.85
        new_range = xr.DataArray(
            np.arange(max_bd + 0.5, max_depth + 0.5, 0.5),
            dims=["range"],
        )

        vel = np.zeros((4, new_range.size, dataset.time.size), dtype="float32")
        # U_mag and U_dir?
        vel_vars = [v for v in dataset.data_vars if "vel" in v]  # type: ignore
        for i, var in enumerate(vel_vars):
            for k in range(dataset.time.size):
                vel[i, :, k] = np.interp(
                    new_range,
                    dataset["bin_depth"][:, k].values,
                    dataset[var][:, k].values,
                )

        U_mag = np.zeros((new_range.size, dataset.time.size), dtype="float32")
        U_dir = np.zeros((new_range.size, dataset.time.size), dtype="float32")
        for k in range(dataset.time.size):
            U_mag[:, k] = np.interp(
                new_range,
                dataset["bin_depth"][:, k].values,
                dataset["U_mag"][:, k].values,
            )
            U_dir[:, k] = np.interp(
                new_range,
                dataset["bin_depth"][:, k].values,
                dataset["U_dir"][:, k].values,
            )

        dataset = dataset.drop_vars(
            vel_vars
            + ["bin_depth", "cell_size", "blank_dist", "range", "U_mag", "U_dir"]
        )

        dataset["vel"] = xr.DataArray(
            vel.astype(np.float32),
            coords={
                "dir": ["E", "N", "U", "Err"],
                "range": new_range,
                "time": dataset.time,
            },
            dims=["dir", "range", "time"],
        )

        dataset["U_mag"] = xr.DataArray(
            U_mag.astype(np.float32),
            coords={
                "range": new_range,
                "time": dataset.time,
            },
            dims=["range", "time"],
        )

        dataset["U_dir"] = xr.DataArray(
            U_dir.astype(np.float32),
            coords={
                "range": new_range,
                "time": dataset.time,
            },
            dims=["range", "time"],
        )

        return dataset

    def hook_finalize_dataset(self, dataset: xr.Dataset) -> xr.Dataset:
        # (Optional) Use this hook to modify the dataset after qc is applied
        # but before it gets saved to the storage area

        # Save mat file
        storage = FileSystem()
        storage.handler.writer = MatlabWriter()

        for var in ["vel", "U_mag", "U_dir"]:
            dataset[var] = dataset[var].astype("float32")

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

        datastream: str = self.dataset_config.attrs.datastream
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

            # Lat/lon
            fig, ax = plt.subplots()

            h = ax.scatter(
                ds["longitude"],
                ds["latitude"],
                c=ds["U_mag"].mean("range").interp(time=ds["time"]).values,
                cmap="Blues",
                s=100,
            )
            fig.colorbar(h, ax=ax, label="Current Speed [m/s]")
            ax.quiver(
                ds["longitude"],
                ds["latitude"],
                ds["vel"][0].mean("range").interp(time=ds["time"]).values,
                ds["vel"][1].mean("range").interp(time=ds["time"]).values,
            )

            ax.set_title("")  # Remove bogus title created by xarray
            ax.set_ylabel("Latitude [deg N]")
            ax.set_xlabel("Longitude [deg E]")
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.4f"))
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.4f"))

            plot_file = self.get_ancillary_filepath(title="location")
            fig.savefig(plot_file)
            plt.close(fig)
