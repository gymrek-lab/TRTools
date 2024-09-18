.. overview_directive
.. |SISTR overview| replace:: SISTR infers selection coefficients STRs
.. overview_directive_done


siSTR
========

|SISTR overview|

(Under development)

SISTR consists of the following subcommands:

* :code:`sistr index`: build ABC and LRT lookup tables needed to run other siSTR commands
* :code:`sistr infer`: run siSTR to compute selection coefficients at individual STRs based on population allele frequencies
* :code:`sistr score`: Use computed selection coefficients to score observed alleles in a cohort
* :code:`sistr joint`: infer siSTR models jointly across loci (SISTR2 method)