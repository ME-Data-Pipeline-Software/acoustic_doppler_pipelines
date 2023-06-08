import xarray as xr
from pathlib import Path
from pydantic import BaseModel, Extra
from typing import Any, Optional

from tsdat import FileWriter
from dolfyn import save_mat


class MatlabWriter(FileWriter):
    """---------------------------------------------------------------------------------
    Writes the dataset to a parquet file.

    Converts a `xr.Dataset` object to a pandas `DataFrame` and saves the result to a
    parquet file using `pd.DataFrame.to_parquet()`. Properties under the
    `to_parquet_kwargs` parameter are passed to `pd.DataFrame.to_parquet()` as keyword
    arguments.

    ---------------------------------------------------------------------------------"""

    class Parameters(BaseModel, extra=Extra.forbid):
        dateNum: Optional[bool] = True

    parameters: Parameters = Parameters()
    file_extension: str = ".mat"

    def write(
        self, dataset: xr.Dataset, filepath: Optional[Path] = None, **kwargs: Any
    ) -> None:
        save_mat(dataset, str(filepath), self.parameters.dateNum)
