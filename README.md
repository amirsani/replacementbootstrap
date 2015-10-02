# replacementbootstrap
The Replacement Bootstrap for Dependent Data

These files are the minimum viable files to run the algorithm from 
- The Replacement Bootstrap for Dependent Data, with Alessandro Lazaric and Daniil Ryabko, ISIT 2015

# Notes
The existing code does not take advantage of any memory optimizations (preferring cache memory when possible) or parallelization (it's naively parallelizable in the number of bootstraps and alphabet size) because the utility depends on the parameter settings. Auto-detection if memory optimization or multiple-cores are possible would be a great upgrade. Additionally, it would be nice to have detection of CUDA and auto-utilization of CUDA if available. The typical user will have a K40. CUDA would be great to run on the loop iterations for a bootstrap or as the loop iteration for the specific alphabet symbol. 

## Minimum viable run file
The run.sh file will demonstrate the code.

# Parameters
## orig_file_name
- The sequence to bootstrap
- Currently a binary encoded csv in the form of a gzip file. 
- Would be nice to include flags for csv or json input files as well, along with a skip header option.

## output_file_name
- Would be nice to include this as an option.

## bootstraps
- Number of bootstraps
- default = 1,000
- minimum = 1,000
- maximum = 100,000

## A
- Alphabet Size must be a power of 2 for now.
- It would be nice to include the option of optimal binning or other unsupervised quantization/discretization.
- The constraint of power of 2 is in the existing code design and not the algorithm.
- default = 2 (binary)
- minimum power value = 2
- maximum power value = 1,024

## n
- Sequence length
- default = autodetect length based on input sequence
- minimum = 10
- maximum = 100000

## replacement_percentage
- The replacement percentage for the algorithm
- default = 100
- minimum = 0
- maximum = 10,000

## possible option
- It would be interesting to allow for other density estimation procedures like kernel density estimation.
