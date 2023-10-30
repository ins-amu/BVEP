B6 BENCHMARK - C IMPLEMENTATION

This is the README file for the C implementation of benchmark B6 in the 
BioPreDyn-Bench collection.
http://www.iim.csic.es/~gingproc/biopredynbench/


***Note: this implementation works under Linux***


INSTALLATION
----------------
(1) install GSL systemwise, or locally. 
(2) go to trunk folder
(3a) if GSL libs are not installed in /usr/lib64, /usr/lib or some similar system folder, the GSL 
variable should be defined in the Makefile.
(3b) if GSL includes are not installed in /usr/local, then you should point to them by adding 
-I/your/include/folder to the INCLUDES variable in Makefile
(4) $make deps (to do only the first time)
(5) $make

CLEANING
---------------
'make clean' -> after cleaning you need to recompile with 'make'
'make veryclean' -> after verycleaning you need to regenerate dependencies with 'make deps' and then recompile with 'make'

RUNNING DEMO
---------
- Demo could be runned just by typing ./ggn in the trunk folder. This will run the model with the standard (hardcoded) set of 
parameters and return the score. 
- The standard set of parameters is already optimized.

USING MODEL DEFINED IN ggn.o 
- To use the model as a library, you must compile the ggn.o, include ggn.h to your code and call the ggn function with two 
parameters: an array of model parameters and an integer representing the size of that array.
- The number of parameters must correspond to the number of "ones" in the mask array (hardcoded in the ggn function). The default 
mask has 37 parameters set to "one".
- The mask is an array of n elements where n is the total number of parameters (fiksed and variable). The mask element is set to 1
if the corresponding parameter is variable (and therefore we can optimize on it), and 0 otherwise (initial value is kept during the 
whole optimization process).
- The model input file dm_hkgn53_wls_5_003 must be in the trunk folder.
- The function will calculate the model and return the score.

INPUT FILE
- various sections are defined in the input file. The ones to take into account when using with different optimisers are
$input - initial set of parameters (fixed and variable). 
$limits - limits for each parameter. Each parameter must be within the corresponding limits, otherwise a VERY_BIG_NUMBER is returned as a score.





