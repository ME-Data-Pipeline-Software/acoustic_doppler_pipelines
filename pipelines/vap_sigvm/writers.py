import xarray as xr
from pathlib import Path
from pydantic import BaseModel, Extra
from typing import Any, Optional, Dict
from tsdat import FileWriter


class CSVWriter(FileWriter):
    """---------------------------------------------------------------------------------
    Writes the dataset to a csv file.

    Converts a `xr.Dataset` object to a pandas `DataFrame` and saves the result to a
    parquet file using `pd.DataFrame.to_csv()`. Properties under the
    `to_csv_kwargs` parameter are passed to `pd.DataFrame.to_csv()` as keyword
    arguments.

    ---------------------------------------------------------------------------------"""

    class Parameters(BaseModel, extra=Extra.forbid):
        to_csv_kwargs: Dict[str, Any] = {}

    parameters: Parameters = Parameters()
    file_extension: str = "csv"

    def write(
        self, dataset: xr.Dataset, filepath: Optional[Path] = None, **kwargs: Any
    ) -> None:

        df = dataset.to_dataframe()
        df.to_csv(filepath, **self.parameters.to_csv_kwargs)  # type: ignore
