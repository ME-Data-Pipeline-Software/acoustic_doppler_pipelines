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


def test_sigvm_pipeline_with_nmea_bt():
    config_path = Path("pipelines/sigvm/config/pipeline.yaml")
    config = PipelineConfig.from_yaml(config_path)
    pipeline = config.instantiate_pipeline()

    test_file = "pipelines/sigvm/test/data/input/102582_20230517T210850UTC.SigVM"
    expected_file = (
        "pipelines/sigvm/test/data/expected/asv.sigvm.b1.20230517.210854_bt.nc"
    )

    dataset = pipeline.run([test_file])
    expected: xr.Dataset = xr.open_dataset(expected_file)  # type: ignore
    assert_close(dataset, expected, check_attrs=False)


def test_sigvm_pipeline_with_nmea_gps():
    config_path = Path("pipelines/sigvm/config/pipeline.yaml")
    config = PipelineConfig.from_yaml(config_path)
    # Change for test
    config.dataset.attrs.vel_bt_correction = 0
    config.dataset.attrs.vel_gps_correction = 1

    pipeline = config.instantiate_pipeline()

    test_file = "pipelines/sigvm/test/data/input/102582_20230517T210850UTC.SigVM"
    expected_file = (
        "pipelines/sigvm/test/data/expected/asv.sigvm.b1.20230517.210854_gps.nc"
    )

    dataset = pipeline.run([test_file])
    expected: xr.Dataset = xr.open_dataset(expected_file)  # type: ignore
    assert_close(dataset, expected, check_attrs=False)
