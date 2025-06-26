import xarray as xr
from pathlib import Path
from tsdat import PipelineConfig, assert_close


def test_winriver2_pipeline():
    config_path = Path("pipelines/winriver2/config/pipeline.yaml")
    config = PipelineConfig.from_yaml(config_path)
    pipeline = config.instantiate_pipeline()

    test_file = "pipelines/winriver2/test/data/input/winriver02_transect.PD0"
    expected_file = (
        "pipelines/winriver2/test/data/expected/pnnl.trdi1200.b1.20100923.210930.nc"
    )

    dataset = pipeline.run([test_file])
    expected: xr.Dataset = xr.open_dataset(expected_file)  # type: ignore
    assert_close(dataset, expected, check_attrs=False)
