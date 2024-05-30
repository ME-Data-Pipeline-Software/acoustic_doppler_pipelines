# Down-looking ADCP Ingestion Pipeline

This pipeline reads in binary files output from a vessel-mounted ADCP, including those read by
WinRiver2 or VMDAS. It is currently set up to read in water-track data from a 4 beam TRDI ADCP. You 
may need to add or remove some variables if they aren't saved, as detailed below.

This README file contains instructions for running and testing this pipeline. Datafiles are saved under
"./storage/root/data" in netCDF4 and MATLAB file formats. Velocity, amplitude and correlation plots are 
saved in the corresponding "./storage/root/ancillary" folder.

## Data Collection
The ADCP should be mounted below water at a -45 degree angle. For a TRDI instrument, this means the 
Y-axis (beam 3 notch) is pointed to port when the ADCP is upside down. For a Nortek instrument, this
means the X-axis is pointed to port when the ADCP is upside down. A bin size of 0.5 m and a blank distance of 0.5 m are good general parameters. Enable bottom track if the water depth is within the ADCP's range, but turn it off if not. Sample as fast as possible, typically 2 Hz.

The GPS should be mounted directly above the ADCP and set it to output NMEA sentences at 2 Hz if 
bottom track is running and 10 Hz if not. Turn on the GGA (position) and VTG (speed) sentences, and HDT
(heading) if possible. If HDT is turned on, make sure to align the GPS heading longitudinally with the 
vessel bow and stern.

## Prerequisites

* Ensure that your development environment has been set up according to
[the instructions](../../README.md#development-environment-setup).

> **Windows Users** - Make sure to run your `conda` commands from an Anaconda prompt OR from a WSL shell with miniconda
> installed. If using WSL, see [this tutorial on WSL](https://tsdat.readthedocs.io/en/latest/tutorials/wsl.html) for
> how to set up a WSL environment and attach VS Code to it.

* Make sure to activate the adcp-pipelines anaconda environment before running any 
commands:  `conda activate adcp-pipelines`


## Editing pipeline data fields
This pipeline is set up to handle data created by an ADCP and GPS recording real-time data to a laptop.

1. There are a number of parameters listed in `retriever.yaml` in lines 15-18 that should be updated.

    a. "depth_offset" is the distance below the waterline that the ADCP transducers sit

    b. "magnetic_declination" is the current magnetic declination. You can look this up online.

2. In `dataset.yaml`, update the information under the `attrs` block per the data collection specifics.

    a. "vel_xx_correction" set to "1" to turn on either bottom track (bt) or GPS motion correction.

3. There is one parameter listed in `shared/quality.yaml` that can be updated:

    a. "correlation_threshold", on line 71, is for a QC test that removes velocity data below a certain % 
    acoustic signal correlation.

4. You may need to update the "trigger" listed in `config/pipeline.yaml` for your ADCP's binary file extension. It will need to be different from other pipelines listed in this repository.


## Running your pipeline
This section shows you how to run the ingest pipeline created by the template.  Note that `{ingest-name}` refers
to the pipeline name you typed into the template prompt, and `{location}` refers to the location you typed into
the template prompt.

1. Make sure to be in the `acoustic_Doppler_pipelines` folder

```bash
cd $PATH_TO_PIPELINE/acoustic_doppler_pipelines
conda create env
conda activate adcp-pipelines
```

2. Run the runner.py with your test data input file as shown below:

```bash
python runner.py ingest $PATH_TO_DATA/{filename}.{ext}
```

3. Finally, data and plots are saved in `acoustic_Doppler_pipelines/storage/root/`
