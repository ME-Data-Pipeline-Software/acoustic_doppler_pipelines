import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from tsdat import TransformationPipeline, FileSystem

from shared.writers import CSVWriter


class VapSigVM(TransformationPipeline):
    """---------------------------------------------------------------------------------
    This is an example pipeline meant to demonstrate how one might set up a
    pipeline using this template repository.

    ---------------------------------------------------------------------------------"""

    def hook_customize_dataset(self, dataset: xr.Dataset) -> xr.Dataset:
        # (Optional) Use this hook to modify the dataset before qc is applied
        return dataset

    def hook_finalize_dataset(self, dataset: xr.Dataset) -> xr.Dataset:
        # (Optional) Use this hook to modify the dataset after qc is applied
        # but before it gets saved to the storage area
        
        # Create additional storage handler to save CSV file(s)
        # All 2D variables should have been dropped
        storage = FileSystem()
        storage.parameters.storage_root = self.storage.parameters.storage_root
        storage.parameters.data_storage_path = self.storage.parameters.data_storage_path
        storage.parameters.data_filename_template = (
            self.storage.parameters.data_filename_template
        )
        storage.handler.writer = CSVWriter()

        # 1D vars
        csv_vars = [
            "depth",
            "heading",
            "pitch",
            "roll",
            "pressure",
            "temperature",
        ]
        ds_csv = dataset[csv_vars]
        
        # Convert -9999 fillvalue to nan for interp and mean
        lat = dataset["latitude_gps"].where(dataset["latitude_gps"]!=dataset["latitude_gps"]._FillValue)
        lon = dataset["longitude_gps"].where(dataset["longitude_gps"]!=dataset["longitude_gps"]._FillValue)
        ds_csv["latitude"] = lat.interp(time_gps=dataset["time"])
        ds_csv["longitude"] = lon.interp(time_gps=dataset["time"])
        ds_csv = ds_csv.drop("time_gps")

        # Add 2D vars
        U_mag = dataset["U_mag"].where(dataset["U_mag"]!=dataset["U_mag"]._FillValue)
        U_dir = dataset["U_dir"].where(dataset["U_dir"]!=dataset["U_dir"]._FillValue)
        ds_csv["U_mag"] = U_mag.mean("range")
        ds_csv["U_dir"] = U_dir.mean("range")

        storage.save_data(ds_csv)

        return dataset

    def hook_plot_dataset(self, ds: xr.Dataset):
        # (Optional, recommended) Create plots.
        def add_colorbar(ax, plot, label):
            cb = plt.colorbar(plot, ax=ax, pad=0.01)
            cb.ax.set_ylabel(label, fontsize=12)
            cb.outline.set_linewidth(1)
            cb.ax.tick_params(size=0)
            cb.ax.minorticks_off()
            return cb

        date = pd.to_datetime(ds["time"].values)

        y_max = 30

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
                nrows=4, ncols=1, figsize=(14, 8), constrained_layout=True
            )

            for beam in range(4):
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
                nrows=4, ncols=1, figsize=(14, 8), constrained_layout=True
            )

            for beam in range(4):
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
                s=200,
            )
            fig.colorbar(h, ax=ax, label="Current Speed [m/s]")
            ax.quiver(
                ds["longitude_gps"][0::10],
                ds["latitude_gps"][0::10],
                ds["vel"][0].mean("range").interp(time=ds["time_gps"]).values[0::10],
                ds["vel"][1].mean("range").interp(time=ds["time_gps"]).values[0::10],
            )

            ax.set_title("")  # Remove bogus title created by xarray
            ax.set_ylabel("Latitude [deg N]")
            ax.set_xlabel("Longitude [deg E]")
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.4f"))
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.4f"))

            plot_file = self.get_ancillary_filepath(title="location")
            fig.savefig(plot_file)
            plt.close(fig)
