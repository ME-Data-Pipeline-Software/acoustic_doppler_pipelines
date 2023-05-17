import xarray as xr
import numpy as np
from numpy.typing import NDArray
import scipy.stats as st
from pydantic import BaseModel, Extra

from tsdat import QualityChecker, QualityHandler
from dolfyn.adv.clean import GN2002, clean_fill


class CheckCorrelation(QualityChecker):
    """----------------------------------------------------------------------------
    Filters out velocity data where correlation is below a
    threshold in the beam correlation data.

    Parameters
    ----------
    ds : xarray.Dataset
      The adcp dataset to clean.
    corr_threshold : numeric
      The maximum value of correlation to screen, in counts or %

    ----------------------------------------------------------------------------"""

    class Parameters(BaseModel, extra=Extra.forbid):
        correlation_threshold: int = 30

    parameters: Parameters = Parameters()

    def run(self, dataset: xr.Dataset, variable_name: str) -> NDArray[np.bool8]:
        mask = dataset["corr"].values < self.parameters.correlation_threshold
        return mask


class CheckOutliers(QualityChecker):
    """---------------------------------------------------------------------------------
    Checks data for elements greater than `n_std` standard deviations away from the mean

    Built-in implementations of quality checkers can be found in the
    [tsdat.qc.checkers](https://tsdat.readthedocs.io/en/latest/autoapi/tsdat/qc/checkers)
    module.

    ---------------------------------------------------------------------------------"""

    class Parameters(BaseModel, extra=Extra.forbid):
        n_std: int = 3

    parameters: Parameters = Parameters()

    def run(self, dataset: xr.Dataset, variable_name: str) -> NDArray[np.bool8]:
        n_std = self.parameters.n_std

        std_dev = dataset[variable_name].std(dim="time", ddof=1)
        mean = dataset[variable_name].std(dim="time")
        mask = dataset[variable_name] > mean + std_dev * n_std

        return mask.data


class CheckGoringNikora2002(QualityChecker):
    """----------------------------------------------------------------------------
    The Goring & Nikora 2002 'despiking' method, with Wahl2003 correction.
    Returns a logical vector that is true where spikes are identified.

    Args:
        variable_name (str): array (1D or 3D) to clean.
        n_points (int) : The number of points over which to perform the method.

    Returns:
        mask [np.ndarray]: Logical vector with spikes labeled as 'True'

    ----------------------------------------------------------------------------"""

    class Parameters(BaseModel, extra=Extra.forbid):
        n_points: int = 5000

    parameters: Parameters = Parameters()

    def run(self, dataset: xr.Dataset, variable_name: str) -> NDArray[np.bool8]:
        return GN2002(dataset[variable_name], npt=self.parameters.n_points)


class CubicSplineInterp(QualityHandler):
    """----------------------------------------------------------------------------
    Interpolate over mask values in timeseries data using the specified method

    Parameters
    ----------
    variable_name : xarray.DataArray
        The dataArray to clean.
    mask : bool
        Logical tensor of elements to "nan" out and replace
    npt : int
        The number of points on either side of the bad values that
    interpolation occurs over
    method : string
        Interpolation scheme to use (linear, cubic, pchip, etc)
    max_gap : int
        Max number of consective nan's to interpolate across, must be <= npt/2

    Returns
    -------
    da : xarray.DataArray
        The dataArray with nan's filled in
    ----------------------------------------------------------------------------"""

    class Parameters(BaseModel, extra=Extra.forbid):
        n_points: int = 12
        method: str = "cubic"

    parameters: Parameters = Parameters()

    def run(
        self, dataset: xr.Dataset, variable_name: str, failures: NDArray[np.bool8]
    ) -> xr.Dataset:
        if failures.any():
            dataset[variable_name] = clean_fill(
                dataset[variable_name],
                mask=failures,
                npt=self.parameters.n_points,
                method=self.parameters.method,
            )
        return dataset
