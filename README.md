**We rewrote the analysis pipeline, which can now be found as the proper [DeePhys package](https://github.com/hornauerp/DeePhys)**
**This repository is no longer maintained and we encourage you to use the [DeePhys repo](https://github.com/hornauerp/DeePhys) instead**




## Welcome to *DeePhys*
The package for **Deep electrophysiological phenotype characterization**:

<img src="https://github.com/hornauerp/EphysDopa/blob/2e43777e3fd1fe0e7467c4b3bf0aa25afb88b602/Figures/EphysDopaSchematic_v2202222.png" alt="Analysis schematic" style="width:500px;"/>

Created with [BioRender](BioRender.com)

## Overview
*DeePhys* was created to facilitate the analysis of extracellular recordings of neuronal cultures using high-density microelectrode arrays (HD-MEAs). MADEB allows users to easily:
- Extract electrophysiological features from spikesorted HD-MEA recordings
- Visualize differential developmental trajectories 
- Apply machine learning algorithms to classify different conditions
- Obtain biomarkers predictive of the respective condition
- Evaluate the effect of treatments

## Requirements
Currently *DeePhys* is only available on MATLAB, so a recent MATLAB installation (>2019b) is required. We plan on expanding *DeePhys* to Python in the near future.

If you want to use the [Notebooks](/Notebooks), you need to install [Jupyter lab or Jupyter notebook](https://jupyter.org/install) (`pip install jupyterlab` or `pip install notebook`) and the [MATLAB kernel](https://pypi.org/project/matlab-kernel/) (`pip install matlab_kernel`) + the [MATLAB API for Python](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html). 

However, all analysis scripts are also available as MATLAB live scripts, which do not require any additional software. 


## Installation
The package is ready-to-use right after cloning. 

## Usage
Code requires spikesorted data in the [phy format](https://github.com/cortex-lab/phy). For help with spikesorting check out the [Spikeinterface package](https://spikeinterface.readthedocs.io/en/latest/). 

Different parts of the analysis are subdevided into different analysis scripts:
- **Feature extraction** ([Jupyter notebook](/Notebooks/feature_extraction.ipynb), [MATLAB live script](/Notebooks/feature_extraction_single_recording.mlx))
- [Data exploration](/Notebooks/data_exploration.ipynb)
- [Classification analysis](/Notebooks/data_exploration.ipynb)
- [Treatment evaluation](/Notebooks/treatment_evaluation.ipynb)


## Citation
This package was published in "Electrophysiological classification of iPSC-derived dopaminergic neurons harbouring the SNCA-A53T mutation" and additionally contains code to replicate the figures used in the publication.

## Disclaimer
This package uses the `readNPY` function provided by the [npy-matlab package](https://github.com/kwikteam/npy-matlab), the `CCG` function provided by the [FMAToolbox](https://github.com/michael-zugaro/FMAToolbox), and the [`othercolor` function](https://ch.mathworks.com/matlabcentral/fileexchange/30564-othercolor).
