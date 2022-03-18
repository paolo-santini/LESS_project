# The LESS project repository

This repository contains the software artifacts related to the LESS 
signature algorithm.

The *reference_code* folder contains a portable C11 implementation
of the LESS signature algorithm for the three parameter sets specified
in the extended paper reference. The reference implementation provides
a functionality test binary, *less_test*, and a microbenchmark for
x86_64 platforms, *less_benchmark*.

The *attacks* folder contains Sage scripts to simulate codewords
enumerating methods to attack the Code Equivalence Problem; scripts
to estimate the attacks cost are provided as well.


## Paper references
* _Alessandro Barenghi, Jean-François Biasse, Edoardo Persichetti, Paolo Santini_:  
**LESS-FM: Fine-Tuning Signatures from the Code Equivalence Problem**. [PQCrypto 2021](https://dblp.org/db/conf/pqcrypto/pqcrypto2021.html#BarenghiBPS21)
 (Extended version available at [https://eprint.iacr.org/2021/396.pdf])
