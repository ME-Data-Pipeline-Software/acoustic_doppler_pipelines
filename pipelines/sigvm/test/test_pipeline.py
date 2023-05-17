import xarray as xr
from pathlib import Path
from tsdat import PipelineConfig, assert_close


def test_sigvm_pipeline():
    config_path = Path("pipelines/sigvm/config/pipeline.yaml")
    config = PipelineConfig.from_yaml(config_path)
    pipeline = config.instantiate_pipeline()

    test_file = "pipelines/sigvm/test/data/input/102582_20220518T212842UTC.SigVM"
    expected_file = "pipelines/sigvm/test/data/expected/asv.sigvm.b1.20220518.212846.nc"

    dataset = pipeline.run([test_file])
    expected: xr.Dataset = xr.open_dataset(expected_file)  # type: ignore
    assert_close(dataset, expected, check_attrs=False)
