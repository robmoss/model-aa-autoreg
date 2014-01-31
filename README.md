Introduction
=============================================================================

This repository contains the source code for an implementation of the
spatially distributed mathematical model of the rat renal juxtamedullary
afferent arteriole described in the following article:

R. Feldberg, M. Colding-Jørgensen, N.-H. Holstein-Rathlou.
Analysis of interaction between TGF and the myogenic response in renal blood
flow autoregulation.
Am J Physiol Renal Physiol 269(4): F581--F593, 2013.
PubMed ID: [7485545](http://www.ncbi.nlm.nih.gov/pubmed/7485545)

This model has also been extended with a permeable glomerulus that calculates
the hydrostatic pressure profile along the glomerular capillary bed and the
single-nephron glomerular filtration rate.

Licence
=============================================================================

This work is made available under the BSD 3-Clause license (see `LICENSE`).

Requirements
=============================================================================

Version numbers listed below are **known-good**.
Newer versions should also be fine; older versions may or may not work.

- A Fortran compiler (gfortran 4.4.7).

- The `R` statistical computing environment (version 2.15).

- The `ggplot2` plotting library for `R` (version 0.9.3).

Structure
=============================================================================

- `model_base.f90`: an implementation of the original afferent model.

- `model_glom.f90`: the extended model with an explicit glomerulus.

- `sngfr.f90`: generate lookup tables for the afferent arteriole response
   that can be used to determine SNGFR with no computational cost.

- `tests.f90`: test cases that reproduce figures from the original paper.

- `cmdline.f90`: parses command-line arguments.

- `main.f90`: model entry point.

Usage
=============================================================================

1. Compile the model by running `make` from the root directory.

2. Run the model by executing `./bin/run_model`.
   This will print out command-line documentation.

3. Run test cases that reproduce several figures from the original paper with
   `./bin/run_model test`.
   The test cases will produce data files (`*.ssv`) in the working directory.

4. Plot the results of the test cases by executing `./bin/plot.R`.

5. Generate a lookup table for model SNGFR by executing
   `./bin/run_model sngfr --glom`.
   This can take a **long** time to complete.
