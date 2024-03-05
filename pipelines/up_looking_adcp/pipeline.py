import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mhkit.dolfyn.adp import api
from tsdat import IngestPipeline, FileSystem

from shared.writers import MatlabWriter


class UpLookingADCP(IngestPipeline):
    """---------------------------------------------------------------------------------
    This is an example ingestion pipeline meant to demonstrate how one might set up a
    pipeline using this template repository.

    ---------------------------------------------------------------------------------"""

    def hook_customize_dataset(self, dataset: xr.Dataset) -> xr.Dataset:
        # (Optional) Use this hook to modify the dataset before qc is applied

        # Locate surface using pressure data and remove data beyond it
        api.clean.nan_beyond_surface(dataset, inplace=True)

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

        datastream: str = self.dataset_config.attrs.datastream
        date = pd.to_datetime(ds["time"].values)

        plt.style.use("default")  # clear any styles that were set before
        plt.style.use("shared/styling.mplstyle")

        y_max = int(ds["depth"].max() * 1.1)

        with plt.style.context("shared/styling.mplstyle"):
            # Current
            fig, ax = plt.subplots(
                nrows=2, ncols=1, figsize=(14, 8), constrained_layout=True
            )

            velE = ax[0].pcolormesh(
                date,
                ds["range"],
                ds["vel"][0],
                cmap="coolwarm",
                shading="nearest",
            )
            ax[0].plot(date, ds["depth"])
            ax[0].set(xlabel="Time (UTC)", ylabel=r"Range [m]", ylim=(0, y_max))
            add_colorbar(ax[0], velE, r"Velocity East [m/s]")
            # velE.set_clim(-3, 3)

            velN = ax[1].pcolormesh(
                date,
                ds["range"],
                ds["vel"][1],
                cmap="coolwarm",
                shading="nearest",
            )
            ax[1].plot(date, ds["depth"])
            ax[1].set(xlabel="Time (UTC)", ylabel=r"Range [m]", ylim=(0, y_max))
            add_colorbar(ax[1], velN, r"Velocity North [m/s]")
            # velN.set_clim(-3, 3)

            plot_file = self.get_ancillary_filepath(title="current")
            fig.savefig(plot_file)
            plt.close(fig)

            # Amplitude
            fig, ax = plt.subplots(
                nrows=ds.n_beams, ncols=1, figsize=(14, 8), constrained_layout=True
            )

            for beam in range(ds.n_beams):
                amp = ax[beam].pcolormesh(
                    date, ds["range"], ds["amp"][beam], shading="nearest"
                )
                ax[beam].set_title("Beam " + str(beam + 1))
                ax[beam].set(xlabel="Time (UTC)", ylabel=r"Range [m]", ylim=(0, y_max))
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
                    ds["range"],
                    ds["corr"][beam],
                    cmap="copper",
                    shading="nearest",
                )
                ax[beam].set_title("Beam " + str(beam + 1))
                ax[beam].set(xlabel="Time (UTC)", ylabel=r"Range [m]", ylim=(0, y_max))
                add_colorbar(ax[beam], corr, "Correlation [%]")

            plot_file = self.get_ancillary_filepath(title="correlation")
            fig.savefig(plot_file)
            plt.close(fig)
