# scratchIBM

Supplementary data, and code, for 

["Browning AP, Jin W, Plank MJ, Simpson MJ (2019) Identifying density-dependent interactions in collective cell behaviour" _J R Soc Interface_ __17__:20200143](https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0143#d3e4763).

[Preprint available on bioRxiv](https://www.biorxiv.org/content/10.1101/811257v1).


## Experimental data
Raw data in the `/Data` folder corresponds to experiment indices in the main document as follows:

| Index |  Folder  |
|:-----:|:--------:|
|   1   |  8000_E2 |
|   2   |  8000_G2 |
|   3   | 10000_H1 |
|   4   | 10000_H2 |
|   5   |  8000_F2 |
|   6   | 10000_A2 |
|   7   | 12000_H3 |
|   8   | 12000_F3 |
|   9   | 12000_D2 |

Experimental data is recorded 12 h after seeding. Hence, raw data indicated by filename ending corresponds to the following times:
 
 | Filename end | Time |
 |:------------:|:----:|
 |    "_12h"    |  0 h |
 |    "_30h"    | 18 h |
 |    "_48h"    | 36 h |


For more details, see `ExperimentalData.xlsx`


## Code

### Requirements
* MATLAB R2018a or later
* MEX and C compiler
* (optional) Parallel computing toolbox
      If you do not have this toolbox, you must replace `parfor` with `for` in the following files:
	1. `ABC_SimulateModel.m`
	2. `ABC_SMC.m`
  
  
### Installation

Before use, you must compile the MEX code contained in `IBM.c` using the 
`Install.m` script or typing:

    mex IBM.c  -R2018a
    
into the command window.


### Examples

#### Simulating IBM
See `Example_IBM.m` for an example IBM realisation.

#### ABC Rejection
* Parallel computing toolbox highly recommend
* This function can take a large amount of time to run for `Nsamples > 100`
  * We have set `Nsamples = 100` and `alpha = 0.1` for demonstrative purposes

Run `Example_ABCRejection.m` code (in the main document, we use `Nsamples = 100000`)

#### ABC SMC Model Selection
* Parallel computing toolbox highly recommend
* This function can take a large amount of time to run for `N > 10`
  * We have set `N = 8` and not used the full sequence of thresholds for demonstrative purposes

Run `Example_ABCSMC.m` code (in the main document, we use `Nsamples = 5000`)
