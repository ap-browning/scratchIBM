****************************
CODE
READ ME
****************************



****************************
REQUIREMENTS
***

- MATLAB R2018a or later
- MEX and C compiler
- (optional) Parallel computing toolbox
      If you do not have this toolbox, you must replace 'parfor' with 'for' in the following files:
	1. ABC_SimulateModel.m
	2. ABC_SMC.m


****************************
INSTALLATION
***

Before use, you must compile the MEX code contained in IBM.c using the 'Install.m' script or typing:
	"mex BinnyIBM.c  -R2018a"
into the command window.


****************************
SIMULATING MODEL
***

See Example_IBM.m for an example IBM realisation


****************************
ABC Rejection 
***
( ! Important ! )
	- Parallel computing toolbox highly recommend
	- This function can take a large amount of time to run for N > 100
	- We have set N = 100 and alpha = 0.1 for demonstrative purposes

Run Example_ABCRejection.m code (in the main document, we use N = 100,000)


****************************
ABC SMC Model Selection
***
( ! Important ! )
	- Parallel computing toolbox highly recommend
	- This function can take a large amount of time to run
	- We have set N = 100 and not used the full sequence of thresholds for demonstrative purposes

Run Example_ABCSMC.m code (in the main document, we use Nsamples = 5,000)
