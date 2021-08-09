# ERGO: Efficient Recurrent Graph Optimized Emitter Density Estimation in Single Molecule Localization Microscopy
[![Build Status](https://travis-ci.com/bencardoen/ERGO.jl.svg?token=Vgy3r9mpYNd3nULfTYZT&branch=main)](https://travis-ci.com/bencardoen/ERGO.jl)

This repository holds the source accompanying our [IEEE Transactions in Medical Imaging 2020 paper](http://www.cs.sfu.ca/~hamarneh/ecopy/tmi2020.pdf).

### Repository organization

[ERGO](https://github.com/bencardoen/ERGO) has 2 stages:
  - Localization : [ERGO.jl](https://github.com/bencardoen/ERGO.py)
  - Counting/Density estimation : [ERGO.py](https://github.com/bencardoen/ERGO.py)

This repository holds the Julia code for the emitter localization stage (GO).

All files in this repository are licensed under **Affero GPL v 3**, copyright 2018-2021 Ben Cardoen.
The software was developed in a multidisciplinary collaboration between the labs of Prof. Ghassan Hamarneh, Prof. Ivan Robert Nabi, and Prof. Keng C. Chou.

### Repository organization
This project (ERGO) has 2 stages:
- Localization : ER**GO**.jl (this repository)
- Counting/Density estimation : **ER**GO.py https://github.com/bencardoen/ERGO.py

All files, with exception of the files in ./data, are licensed under [**Affero GPL v 3**](https://www.gnu.org/licenses/agpl-3.0.txt), copyright 2018-2021 [Ben Cardoen](https://orcid.org/0000-0001-6871-1165).

The software was developed in a multidisciplinary collaboration between the labs of [Prof. Ghassan Hamarneh](https://www.medicalimageanalysis.com/ghassans-bio), [Prof. Ivan Robert Nabi](https://www.bme.ubc.ca/person/ivan-nabi/), and [Prof. Keng C. Chou](https://www.chem.ubc.ca/keng-chou). This project could not have been realized without my other co-authors: Hanene Ben Yedder, and Anmol Sharma.

#### Visualization of algorithm

![Visualization of the ROI localization algorithm](figure3.png)


The manuscript explains the algorithm in more detail.

![ERGO.jl applied to a single frame.](figure4.png)
(Above) ERGO.jl applied to a single frame from a sequence. Left is the raw frame, and right is the processed frame. White bounding boxes indicate the ROIs ERGO centers around suspected emissions. Green are groundtruth locations of emissions (can be very faint).  Purple lines indicate the intensity weighting. Yellow are local maxima. Red is a 'noise' or 'negative' ROI, selected to contain with high confidence no emissions.

'Negative' ROI are useful for:

- Avoid selecting too much of the frame as signal
- Producing training data for the ER (deep learning emission density estimation) part.

Overlapping ROIs are intentional, to ensure each ROI is centered on a presumed emission. Postprocessing (density estimation, localization) can later process or reduce these ROIs.

Note how ERGO.jl only misses 1 emission, that is in intensity indistinguishable from background (top left corner).

#### Potential use cases:
- denoising : filter frame where density = 0
- enhance localization : knowing only k emitters are present, the localization problem becomes easier
- (real time) acquisition optimization : e.g. optimize 1 emitter per ROI for optimal localization & ideal signal separation

##### Performance
On a not so fast laptop with SSD, it's not unusual to get 250FPS in processing speed. Note that the code can be parallelized, and optimized in multiple places, drastically increasing FPS, but FPS of exceeding real time (60) is sufficient for our applications.


### Layout
- src
  - module source code
- notebooks
  - example usage (need jupyter)
To run the notebook examples, you need [Jupyter](https://jupyter.org/) and a [Julia kernel](https://github.com/JuliaLang/IJulia.jl).

### Data
You can download the data used at : http://bigwww.epfl.ch/smlm/challenge2016/datasets/MT0.N1.HD/Data/data.html

For your convenience, the data is provided (I make no claim or license to this data, rights belong to original owners) in folder ./data.

### Installation
You need [Julia 1.6.2](https://julialang.org/downloads/).
This project was originally developed on v1.4.x, but I recommend v1.6.x for speed and features.
Note that if you decide to use older version, please run the tests before you do anything else.

The below commands show the fastest way to get started.
```bash
mkdir newdir
cd newdir
julia> ]
> activate .
add https://github.com/bencardoen/ERGO.jl
instantiate
test ERGO
```
The ']' launches package manager, here we install in a new environment to ensure dependencies don't disturb your global Julia installation.
If the tests do not succeed, please make an issue with reproducible error.

The below screen recording shows how these instructions should look in execution

![Installation process](install.gif)

### Running
See notebooks/example.ipynb for a workflow used in the paper

Note that the notebook has extra dependencies (see notebook) not needed to run ERGO, but needed to display statistics, images, etc. Instructions to install those are at the top of the notebook.

Note: we run automated tests with travis on Windows and Linux, but Linux is the platform used to develop ERGO.jl.

If you use, or find this work useful, please cite the below paper:

```bibtex
@article{cardoen2019ergo,
  title={Ergo: efficient recurrent graph optimized emitter density estimation in single molecule localization microscopy},
  author={Cardoen, Ben and Yedder, Hanene Ben and Sharma, Anmol and Chou, Keng C and Nabi, Ivan Robert and Hamarneh, Ghassan},
  journal={IEEE transactions on medical imaging},
  volume={39},
  number={6},
  pages={1942--1956},
  year={2019},
  publisher={IEEE}
}

```


Note : reproduced figures are drafts of pre-submission versions, not published originals, and therefore under our copyright.
