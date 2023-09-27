
# FUGAsseM History #

## v0.3.5 2023-03-17 ##
* Fixed issue when using '--bypass-coexp' option
* Tweaked default options for fugassem_generate_annotation_input utility
* Update demo run
* Tweaked workflow name
* Release v0.3.5

## v0.3.4 2022-11-11 ##
* Renamed options for '--go-mode'
* Added 'bypass-coexp' option to skip calculating coexpression
* Fixed issue when using '--bypass-prediction' option
* Added options for calculating feature vector

## v0.3.3 2022-08-16 ##
* Reorganize output folders
* Added a finalized combined results file
* Added stratified option for 'Terminal' or 'Strain' level
* Fixed issue when filtering based on prevalence
* Fixed issue when merging results

## v0.3.2 2022-04-08 ##
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
