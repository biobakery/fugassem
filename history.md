
# FUGAsseM History #

## v0.3.2 2022-04-08 ##
* Updated GO-obo files
* Added option for prediction using universal GO set 
* Added option for prediction using global GO set
* Added option for predicting terms from all GO ontologies
* Refined generate_annotation_input module to prepare co-homology evidence 
* Added function for visualization
* Fixed issue on processing fugassem_stratified_abundance module
* Added example data

## v0.3.1 2022-03-02 ##
* Fixed small issue from geneontology.py
* Added function for detailed annotation when '--go-level none'
* Tweaked stratified_abundance.py for accepting compressed input files
* Added format_function module for formatting input file
* Fixed issue when using different zero-filtering approach

## v0.3.0 2022-02-05 ##
* Accepted input matrix evidence with connection weights
* Added options to accept either pair- or unpair-based matrix evidence

## v0.2.1 2022-01-08 ##
* Fixed issues of utilities for preparing input files
* Added option for skipping to integrate mtx into finalized prediction 

## v0.2.0 2021-12-14 ##
* Added multi-tasking module for processing taxa in parallel

## v0.1.1 2021-12-10 ##
* Added parameters for taking vector-based and matrix-based evidence files

## v0.1.0 2021-12-06 ##
* Added initial preprocessing module, prediction module and utilitites of FUGAsseM
