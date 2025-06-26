import xarray as xr
from pathlib import Path
from tsdat import PipelineConfig, assert_close


def test_vmdas_pipeline():
    config_path = Path("pipelines/vmdas/config/pipeline.yaml")
    config = PipelineConfig.from_yaml(config_path)
    pipeline = config.instantiate_pipeline()

    test_file = "pipelines/vmdas/test/data/input/clllm_20250530T065119_009_000000.ENS"
    expected_file = (
        "pipelines/vmdas/test/data/expected/clllm.trdi300.a1.20250530.135604.nc"
    )

    dataset = pipeline.run([test_file])
    expected: xr.Dataset = xr.open_dataset(expected_file)  # type: ignore
    assert_close(dataset, expected, check_attrs=False)
