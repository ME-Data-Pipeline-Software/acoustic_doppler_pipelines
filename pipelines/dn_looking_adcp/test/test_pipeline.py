import xarray as xr
from pathlib import Path
from tsdat import PipelineConfig, assert_close


def test_dn_looking_ADCP_pipeline():
    config_path = Path("pipelines/dn_looking_adcp/config/pipeline.yaml")
    config = PipelineConfig.from_yaml(config_path)
    pipeline = config.instantiate_pipeline()

    test_file = "pipelines/dn_looking_adcp/test/data/input/winriver02_transect.PD0"
    expected_file = "pipelines/dn_looking_adcp/test/data/expected/pnnl.trdi1200.b1.20100923.210930.nc"

    dataset = pipeline.run([test_file])
    expected: xr.Dataset = xr.open_dataset(expected_file)  # type: ignore
    assert_close(dataset, expected, check_attrs=False, atol=1e-5)
