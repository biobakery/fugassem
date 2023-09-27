# FUGAsseM User Manual

**FUGAsseM** (**F**unction predictor of **U**ncharacterized **G**ene products by **Asse**ssing high-dimensional community data in **M**icrobiomes) is a computational tool based on a “guilt by association” approach to predict functions of novel gene products in the context of microbial communities. It uses machine learning methods to predict functions of microbial proteins by integrating multiple types of community-based data, such as co-expression patterns from metatranscriptomics (MTX), co-localization patterns from metagenomic (MGX) assemblies, homologous-based annotation, domain-based annotations, etc.


## Citing FUGAsseM

**If you use FUGAsseM, please cite our manuscript:**

Yancong Zhang, Amrisha Bhosle, Sena Bae, Kelly Eckenrode, Xueying (Sonia) Huang, Jingjing Tang, Danylo Lavrentovich, Lana Awad, Ji Hua, Xochitl C. Morgan, Andy Krueger, Wendy S. Garrett, Eric A. Franzosa\*, Curtis Huttenhower\*. "Predicting functions of uncharacterized gene products in microbiomes" *[In preparation]*.

**Please add the software link in your Methods if you cite FUGAsseM:**

[http://huttenhower.sph.harvard.edu/fugassem](http://huttenhower.sph.harvard.edu/fugassem)

**For additional information, read the [FUGAsseM Tutorial](https://github.com/biobakery/biobakery/wiki/FUGAsseM)**


If you have questions about FUGAsseM, please direct it to [the FUGAsseM channel](https://forum.biobakery.org/c/Microbial-community-profiling/fugassem) of the bioBakery Support Forum.

***

## Contents ##
* [Workflow](#workflow)
	* [Workflow by bypass mode](#workflow-by-bypass-mode) 
* [Install FUGAsseM](#install-fugassem)
    * [Requirements](#requirements)
    * [Installation](#installation)
    	* [Install FUGAsseM](#install-fugassem)
* [How to run](#how-to-run)
    * [Basic usage](#basic-usage)
    * [Demo runs](#demo-runs)
    	* [Demo run of canonical prediction](#demo-run-of-canonical-prediction) 
    		* [Input files](#input-files)
    		* [Running command](#running-command)
    		* [Run FUGAsseM-MTX model](#run-fugassem-mtx-model)
    		* [Output files](#output-files)
    	* [Demo run of integrated prediction](#demo-run-of-integrated-prediction) 
    		* [Input files](#input-files)
    		* [Running command](#running-command)
    		* [Run FUGAsseM-full model](#run-fugassem-full-model)
    		* [Output files](#output-files)
* [Guides to FUGAsseM Utilities](#guides-to-fugassem-utilities)
	* [Preparing stratified MTX-based abundance input](#preparing-stratified-mtx-based-abundance-input)
		* [Input files of preparing MTX abundance utility](#input-files-of-preparing-MTX-abundance-utility)
		* [Demo run of preparing MTX abundance utility](#demo-run-of-preparing-MTX-abundance-utility) 
		* [Output files of preparing MTX abundance utility](#output-files-of-preparing-MTX-abundance-utility)
	* [Preparing evidence input](#preparing-evidence-input)
		* [Input files of preparing evidence utility](#input-files-of-preparing-evidence-utility)
		* [Demo run of preparing evidence utility](#demo-run-of-preparing-evidence-utility) 
		* [Output files of preparing evidence utility](#output-files-of-preparing-evidence-utility)
    	
***


## Workflow
### Description
FUGAsseM predicts functions for uncharacterized protein families in host- and non-host-associated microbial communities. This software starts with protein families from communities and combines expression profiles in MTX, sequence context in MGX, and secondary-structure-based functional annotations or interactions to infer protein functions at the community-wide scale.


### Workflow by bypass mode
There are multiple bypass options that will allow you to adjust the standard workflow.

Bypass options:

* --bypass-preparing-taxa
	* do not run the module for splitting stratified MTX into individual taxa
* --bypass-preprocessing
	* do not run the module for preprocessing evidences and functions for prediction
* --bypass-prediction
	* do not run the module for predicting functions
* --bypass-mtx
	* do not integrate MTX for finalized predicting functions
* --bypass-coexp
	* do not calculate MTX-based coexpression


## Install FUGAsseM
### Requirements
1. [Python](https://www.python.org/) (version >= 3.7, requiring numpy, pandas
multiprocessing, sklearn, matplotlib, scipy, goatools, statistics python packages; *tested 3.7*)
2. [AnADAMA2](https://huttenhower.sph.harvard.edu/anadama2) (version >= 0.8.0; *tested 0.8.0*)


### Installation
You only need to do **any one** of the following options to install the FUGAsseM package. 	

**Option 1: Installing with conda**

* `$ conda install -c biobakery fugassem`

**Option 2: Installing with pip**

* `$ pip install fugassem`
* If you do not have write permissions to `/usr/lib/`, then add the option --user to the install command. This will install the python package into subdirectories of `~/.local/`. Please note when using the --user install option on some platforms, you might need to add `~/.local/bin/` to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message `fugassem: command not found` when trying to run FUGAsseM after installing with the --user option.

**Option 3: Installing from source**

1. Download FUGAsseM

	You can download the latest FUGAsseM release or the development version by **any one** of the following options. The source contains example files.
	
	* Option 1: Latest Release (Recommended)
	
		Download [fugassem.zip](https://github.com/biobakery/fugassem/archive/master.zip) and unpack the latest release of FUGAsseM.
	
	* Option 2: Development Version
	
		Create a clone of the repository (Note: Creating a clone of the repository requires Github to be installed):

		`$ git clone https://github.com/biobakery/fugassem`

2. Move to the FUGAsseM directory
	* `$ cd $FUGAsseM_PATH`

3. Install FUGAsseM package
	* `$ python setup.py install`
	* If you do not have write permissions to '/usr/lib/', then add the option --user to the install command. This will install the python package into subdirectories of '\~/.local'. Please note when using the '--user' install option on some platforms, you might need to add '\~/.local/bin/' to your $PATH as it might not be included by default. You will know if it needs to be added if you see the following message `fugassem: command not found` when trying to run FUGAsseM after installing with the '--user' option.
 

## How to run
### Basic usage
* For a list of command line options, run:

	`$ fugassem --help`

	This command yields:

	```
	usage: fugassem_workflow.py [-h] [--version]
                        [--taxon-level {MSP,Terminal,Species,Genus,Family,Order,Class,Phylum}]
                        [--minimum-prevalence MINIMUM_PREVALENCE]
                        [--minimum-abundance MINIMUM_ABUNDANCE]
                        [--minimum-detected MINIMUM_DETECTED]
                        [--minimum-coverage MINIMUM_COVERAGE]
                        [--minimum-number MINIMUM_NUMBER]
                        [--filtering-zero {lenient,semi-strict,strict,None}]
                        [--covariate-taxon COVARIATE_TAXON]
                        [--correlation-method {Pearson,Spearman,Kendall,Pearson_SE}]
                        [--go-level GO_LEVEL]
                        [--go-mode {taxon-specific,universal,union}]
                        [--func-type {GO,BP,CC,MF}] [--ml-type {RF,NB,DT}]
                        [--vector-list VECTOR_LIST]
                        [--matrix-list MATRIX_LIST] [--matrix-pair]
                        [--basename BASENAME] --input-annotation
                        INPUT_ANNOTATION [--bypass-preparing-taxa]
                        [--bypass-preprocessing] [--bypass-prediction]
                        [--bypass-mtx] [--bypass-coexp] [--threads THREADS] [--memory MEMORY]
                        [--time TIME] [--output OUTPUT] [-i INPUT]
                        [--config CONFIG] [--local-jobs JOBS]
                        [--grid-jobs GRID_JOBS] [--grid GRID]
                        [--grid-partition GRID_PARTITION]
                        [--grid-benchmark {on,off}]
                        [--grid-options GRID_OPTIONS]
                        [--grid-environment GRID_ENVIRONMENT]
                        [--grid-scratch GRID_SCRATCH] [--dry-run]
                        [--skip-nothing] [--quit-early]
                        [--until-task UNTIL_TASK]
                        [--exclude-task EXCLUDE_TASK] [--target TARGET]
                        [--exclude-target EXCLUDE_TARGET]
                        [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]

	```

* **Option: --go-mode**
	
	Mode of Gene Ontology (GO) set used for prediction:
	
	* "taxon-specific": taxon-specific informative terms (where "informative" term is defined as a term that contains >X genes and all of its children terms contain < X genes)
	* "universal": universal informative terms 
	* "union": merged taxon-specific informative terms for overall prediction *[Recommended]*
                       

* **Option: --go-level** 

	GO informative level used for trimming terms that are informative at a given level:
	
	* \<number OR fraction of genes\>: specify numeric level for trimming
	* \<all\>: keep all terms
	* \<none\>: skip trimming               
                        
* **Option: --vector-list**

	\<string\> a comma separated list of vector-based evidence files, such as protein-over-function type of files.
	
* **Option: --matrix-list**

	\<string\> a comma separated list of matrix-based evidence files, such as protein-over-protein network type of file).

* **Option: --minimum-prevalence**

 	\<fraction\> minimum prevalence of each protein family in normalized MTX stratified by taxa scaling to 1.

* **Option: --minimum-coverage**

	\<fraction\> minimum fraction of annotated genes per taxon.

* **Option: --minimum-number**

	\<integral number\> minimum number of total genes per taxon.

* **Parallelization Options**

	When running any workflow you can add the following command line options to make use of existing computing resources:
	* --threads <1>  number of threads/cores for each task to use
	* --local-jobs <1> : Run multiple tasks locally in parallel. Provide the max number of tasks to run at once. The default is one task running at a time.
	* --grid-jobs <0> : Run multiple tasks on a grid in parallel. Provide the max number of grid jobs to run at once. The default is zero tasks are submitted to a grid resulting in all tasks running locally.
	* --grid \<slurm> : Set the grid available on your machine. This will default to the grid found on the machine with options of slurm and sge.
	* --grid-partition \<serial_requeue> : Jobs will be submitted to the partition selected. The default partition selected is based on the default grid.

	For additional workflow options, see the [AnADAMA2](https://github.com/biobakery/anadama2) user manual.


### Demo runs
#### Demo run of canonical prediction

A standard workflow of FUGAsseM using MTX model by taking MTX coexpression profiles and raw GO annotations as inputs, which (1) preparing protein families and annotations, (2) building coexpression profiles of proteins within each taxon, and (3) building a machine learning classifier for function prediction.

##### Input files
* normalized MTX-based abundance table stratified by taxa (TSV format file), e.g. [demo\_proteinfamilies\_rna\_CPM.stratified\_Species_mtx.tsv](https://github.com/biobakery/fugassem/blob/master/examples/input/demo_proteinfamilies_rna_CPM.stratified_Species_mtx.tsv)
* raw annotations for protein families (TSV format file), e.g. [demo_proteinfamilies.GO.simple.tsv](https://github.com/biobakery/fugassem/tree/master/examples/input/demo_proteinfamilies.GO.simple.tsv)
	
##### Running command
	
```
$ fugassem --basename $BASENAME \ 
--input $INPUT_MTX \
--input-annotation $INPUT_annotation \ 
--output $OUTPUT_DIR
```

* The command replaces `$INPUT_MTX `, `$INPUT_annotation` with two input files, `$OUTPUT_DIR` with the path to the folder to write output files. See the section on **parallelization options** to optimize the run based on your computing resources. 
* This runs with the default settings to run all modules. These settings will work for most data sets. However, if you need to customize your settings to modify the default settings. You can customize which modules you want to run in your own local configuration file.
	
##### Run FUGAsseM-MTX model

```
$ fugassem --basename $BASENAME \ 
--input input/demo_proteinfamilies_rna_CPM.stratified_Species_mtx.tsv \ 
--input-annotation input/demo_proteinfamilies.GO.simple.tsv \ 
--output $OUTPUT_DIR
```

##### Output files
When FUGAsseM is run, the merged prediction files of all taxa will be created at `$OUTPUT_DIR/merged`:

**1. Merged finalized prediction file**
		
```
taxon   feature func    category        score   raw_ann
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0000155      GO      0.75    1
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0003700      GO      0.73    1
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0003824      GO      0.22    0
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0004673      GO      0.8     1
...
Escherichia_coli        Cluster_100559  GO:0000271      GO      0.45    0
Escherichia_coli        Cluster_100559  GO:0003677      GO      0.35    0
Escherichia_coli        Cluster_100559  GO:0003700      GO      0.61    0
Escherichia_coli        Cluster_100559  GO:0003979      GO      0.34    0
...
```
	
* File name: `$OUTPUT_DIR/merged/$BASENAME.finalized_ML.prediction.tsv`
* This file includes the finalized predictions by MTX coexpression (TSV format file).
* `$OUTPUT_DIR` = the output folder
* `$BASENAME` = the basename of output files
* This file details the prediction of each protein family per each function per taxon.
* The predictions for each protein family combined by coexpression from MTX

**2. Finalized Prediction file of each taxon**

```	
taxon   feature func    category        score   raw_ann
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0000155      GO      0.75    1
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0003700      GO      0.73    1
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0003824      GO      0.22    0
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0004673      GO      0.8     1
...
```
				
* The finalized prediction resDomainults using integrated evidence per taxon are in the file: `$OUTPUT_DIR/main/$TAXON_NAME/$BASENAME.$TAXON_NAME.finalized_ML.prediction.tsv`.
* `$OUTPUT_DIR` = the output folder
* `$TAXON_NAME ` = taxon name
* `$BASENAME` = the basename of output files


#### Demo run of integrated prediction

When other community-wide data are available, FUGAsseM can be run using an integrated model. This includes: (1) preparing protein families and annotations, (2) building co-expression profiles of proteins within each taxon, (3) building individual machine learning classifiers for function prediction per type of evidence data, and (4) integration to generate an ensemble classifier for final function prediction. Evidence such as homology between protein families, gene neighborhood, and domain-domain interactions may be included.

##### Input files
* normalized MTX-based abundance table stratified by taxa (TSV format file), e.g. [demo\_proteinfamilies\_rna\_CPM.stratified\_Species_mtx.tsv](https://github.com/biobakery/fugassem/tree/master/examples/input/demo_proteinfamilies_rna_CPM.stratified_Species_mtx.tsv)
* raw annotations for protein families (TSV format file), e.g. [demo_proteinfamilies.GO.simple.tsv](https://github.com/biobakery/fugassem/tree/master/examples/input/demo_proteinfamilies.GO.simple.tsv)
* vector-based evidence file (TSV format file), e.g. [demo_proteinfamilies.GO.homology.tsv](https://github.com/biobakery/fugassem/tree/master/examples/input/demo_proteinfamilies.GO.homology.tsv)
* matrix-based evidence files (TSV format file), e.g.
	* [demo_proteinfamilies.DDI.simple.tsv](https://github.com/biobakery/fugassem/tree/master/examples/demo_proteinfamilies.DDI.simple.tsv): Domain-domain interactions for building DDI network for prediction
	* [demo_proteinfamilies.contig.simple.tsv](https://github.com/biobakery/fugassem/tree/master/examples/demo_proteinfamilies.contig.simple.tsv): MGX-based contigs of protein families for building co-contig network for prediction

##### Running command
	
```
$ fugassem --basename $BASENAME --input $INPUT_MTX \ 
--input-annotation $INPUT_annotation \ 
--vector-list $VECTOR_list --matrix-list $METRIX_list \ 
--output $OUTPUT_DIR
```

* The command replaces `$INPUT_MTX `, `$INPUT_annotation` with two input files, `$OUTPUT_DIR` with the path to the folder to write output files. See the section on **parallelization options** to optimize the run based on your computing resources. 
* This runs with the default settings to run all modules. These settings will work for most data sets. However, if you need to customize your settings to modify the default settings. You can customize which modules you want to run in your own local configuration file.
	
##### Run FUGAsseM-full model

```
$ fugassem --basename $BASENAME --input input/demo_proteinfamilies_rna_CPM.stratified_Species_mtx.tsv \
--vector-list input/demo_proteinfamilies.GO.homology.tsv \
--matrix-list input/demo_proteinfamilies.DDI.simple.tsv,input/demo_proteinfamilies.contig.simple.tsv \
--input-annotation input/demo_proteinfamilies.GO.simple.tsv \
--output $OUTPUT_DIR
```

##### Output files
When FUGAsseM is run, the merged prediction files of all taxa will be created at `$OUTPUT_DIR/merged`:

**1. Finalized prediction file**
		
```
taxon   feature func    category        score   raw_ann
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0000155      GO      0.97    1
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0003700      GO      0.97    1
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0003824      GO      0.34    0
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0004673      GO      0.73    1
...
Escherichia_coli        Cluster_100559  GO:0000271      GO      0.16    0
Escherichia_coli        Cluster_100559  GO:0003677      GO      0.17    0
Escherichia_coli        Cluster_100559  GO:0003700      GO      0.29    0
Escherichia_coli        Cluster_100559  GO:0003979      GO      0.15    0
...
```
	
* File name: `$OUTPUT_DIR/merged/$BASENAME.finalized_ML.prediction.tsv`
* This file includes the finalized predictions by integrating multiple machine learning (ML) classifiers (TSV format file).
* `$OUTPUT_DIR` = the output folder
* `$BASENAME` = the basename of output files
* This file details the prediction of each protein family per each function per taxon.
* The predictions for each protein family are combined by multiple information sources, e.g. coexpression from MTX, co-localization from MGX, sequence similarity, interactions, etc.
			
**2. Prediction file for each type of evidence**
		
```
taxon   feature func    category        score   raw_ann
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0000155      GO      0.75    1
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0003700      GO      0.73    1
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0003824      GO      0.22    0
Bacteroides_thetaiotaomicron    Cluster_1024034 GO:0004673      GO      0.8     1
...
Escherichia_coli        Cluster_100559  GO:0000271      GO      0.45    0
Escherichia_coli        Cluster_100559  GO:0003677      GO      0.35    0
Escherichia_coli        Cluster_100559  GO:0003700      GO      0.61    0
Escherichia_coli        Cluster_100559  GO:0003979      GO      0.34    0
...
```
	
* File name: `$OUTPUT_DIR/merged/$BASENAME.$EVIDENCE_TYPE_ML.prediction.tsv` (where `$EVIDENCE_TYPE` = the basename of each evidence).
* This file includes the predictions based on individual type of evidence (TSV format file).
* `$OUTPUT_DIR` = the output folder
* `$BASENAME` = the basename of output files
* This file details the prediction of each protein family per each function per taxon.
* The predictions for each protein family using one type of evidence data to build the ML classifier.

**3. Prediction file for each taxon**

* Predictions of each taxon
	* FUGAsseM predicts functions based on input evidence data.
	* The finalized prediction results using integrated evidence per taxon are in the file: `$OUTPUT_DIR/main/$TAXON_NAME/$BASENAME.$TAXON_NAME.finalized_ML.prediction.tsv`.
	* The prediction results by using individual evidence per taxon are in the file: `$OUTPUT_DIR/$TAXON_NAME/$BASENAME.$TAXON_NAME.$EVIDENCE_TYPE_ML.prediction.tsv` (where `$EVIDENCE_TYPE` = the basename of each evidence).

**4. Intermediate output files**
	
* Preprocessing features of each taxon
	* FUGAsseM preprocesses input evidence data and prepares feature tables for machine learning per taxon. Each type of feature will be used to build a ML classifier. These input data are in the folder: `$OUTPUT_DIR/main/$TAXON_NAME/data/`.
	* All intermediate preprocessing results are in the folder: `$OUTPUT_DIR/main/$TAXON_NAME/preprocessing/`.
	* All intermediate prediction results are in the folder per taxon: `$OUTPUT_DIR/main/$TAXON_NAME/prediction/`.



## Guides to FUGAsseM Utilities

### Preparing stratified MTX-based abundance input

FUGAsseM takes a MTX-based abundance (that is normalized within each taxon) table of protein families as the input. Users may provide this table using their own analysis or use one of two options as below:

* **Option 1: Reference-based approach**
	
	The first option uses HUMAnN to map MTX reads against reference proteins and get quantified MTX abundance stratified into species. This stratified MTX abundance should be normalized within each taxon by either [HUMAnN](https://github.com/biobakery/humann)’s utility [humann\_renorm\_table](https://github.com/biobakery/humann#humann_renorm_table) or FUGAsseM’s utility fugassem\_abundance\_normalization. FUGAsseM accepts this normalized MTX abundance stratified into taxon as the input.
	
* **Option 2: Assembled-based approach**

	The second option relies on assembled protein families from MGX. We provide a utility in the FUGAsseM package called fugassem\_generate\_stratified\_mtx\_input to help generate the MTX abundance table. This utility (1) maps MTX shotgun reads against MGX-assembled gene catalogs, (2) sums up the quantified abundance of gene catalogs to protein families level, (3) normalizes the protein-family-based MTX abundance within each stratified taxon. The output of this utility can be taken as FUGAsseM's input.

	#### Generating workflow
	`$ fugassem_generate_stratified_mtx_input --help`
	
	This command yields:
	
	```
	usage: fugassem_generate_stratified_mtx_input [-h] [--version]
                                              [--extension-paired EXTENSION_PAIRED]
                                              [--extension {.fastq.gz,.fastq}]
                                              [--taxon-level {MSP,Species}]
                                              [--taxon-prevalence TAXON_PREVALENCE]
                                              [--normalized-method NORMALIZED_METHOD]
                                              --gene-catalog GENE_CATALOG
                                              --gene-catalog-seq
                                              GENE_CATALOG_SEQ
                                              --protein-family PROTEIN_FAMILY
                                              --family-taxonomy
                                              FAMILY_TAXONOMY
                                              [--basename BASENAME]
                                              [--threads THREADS]
                                              [--memory MEMORY] [--time TIME]
                                              [--output OUTPUT] [-i INPUT]
                                              [--config CONFIG]
                                              [--local-jobs JOBS]
                                              [--grid-jobs GRID_JOBS]
                                              [--grid GRID]
                                              [--grid-partition GRID_PARTITION]
                                              [--grid-benchmark {on,off}]
                                              [--grid-options GRID_OPTIONS]
                                              [--grid-environment GRID_ENVIRONMENT]
                                              [--grid-scratch GRID_SCRATCH]
                                              [--dry-run] [--skip-nothing]
                                              [--quit-early]
                                              [--until-task UNTIL_TASK]
                                              [--exclude-task EXCLUDE_TASK]
                                              [--target TARGET]
                                              [--exclude-target EXCLUDE_TARGET]
                                              [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
	                        
    ```
	
	* `--input`: the input directory where a set of fastq (or fastq.gz) files (single-end or paired-end) passing through QC are stored. The files are expected to be named `$SAMPLE.paired_R1.gz`, `$SAMPLE.paired_R2.gz`, `$SAMPLE.orphan_R1.gz` and `$SAMPLE.orphan_R2.gz` where `$SAMPLE` is the sample name or identifier corresponding to the sequences. `$SAMPLE` can contain any characters except spaces or periods.
	* `--extension-paired` indicates the extension for paired fastq files using comma to separate. It should be specified as ".R1.fastq.gz,.R2.fastq.gz" if the paired fastq files are `$SAMPLE.R1.fastq.gz` and `$SAMPLE.R2.fastq.gz`  
	* `--extension` indicates the extension for all fastq files. It should be specified as ".fastq.gz" if the fastq files are `$SAMPLE.fastq.gz` 
	* `--taxon-level`: taxonomic level used for stratification [ Default: Species ]
	* `--gene-catalog`: input clustering file in extended-fasta format for non-redundant gene catalogs, which can be generated by [MetaWIBELE](https://github.com/biobakery/metawibele)
	* `--gene-catalog-seq`: input fasta file of nucleotide sequences of representatives for non-redundant gene catalogs, which can be generated by [MetaWIBELE](https://github.com/biobakery/metawibele)
	* `--protein-family`: input clustering file in extended-fasta format for protein families clustered by non-redundant gene catalogs, which can be generated by [MetaWIBELE](https://github.com/biobakery/metawibele)
	* `--family-taxonomy`: input taxonomy file for protein families, which can be generated by [MetaWIBELE](https://github.com/biobakery/metawibele)
	* `--output`: the output directory. 
	
	#### Input files of preparing MTX abundance utility
	* QC'ed shotgun sequencing metatranscriptome file (fastq, fastq.gz, fasta, or fasta.gz format), e.g. "raw_input" folder including:
		- [sample1.R1.fastq.gz](https://github.com/biobakery/fugassem/raw/master/examples/raw_input/sample1.R1.fastq.gz)
		- [sample1.R2.fastq.gz](https://github.com/biobakery/fugassem/raw/master/examples/raw_input/sample1.R2.fastq.gz)
		- [sample2.R1.fastq.gz](https://github.com/biobakery/fugassem/raw/master/examples/raw_input/sample2.R1.fastq.gz)
		- [sample2.R2.fastq.gz](https://github.com/biobakery/fugassem/raw/master/examples/raw_input/sample2.R2.fastq.gz)
	* clustering file of non-redundant gene catalogs (extended-fasta format): [demo_genecatalogs.clstr](https://github.com/biobakery/fugassem/raw/master/examples/raw_input/demo_genecatalogs.clstr)
	* nucleotide sequences of representatives for non-redundant gene catalogs (fasta format): [demo_genecatalogs.centroid.fna](https://github.com/biobakery/fugassem/raw/master/examples/raw_input/demo_genecatalogs.centroid.fna)
	* clustering file of protein families clustered by non-redundant gene catalogs (extended-fasta format): [demo_proteinfamilies.clstr](https://github.com/biobakery/fugassem/raw/master/examples/raw_input/demo_proteinfamilies.clstr)
	* taxonomy file of protein families clustered by non-redundant gene catalogs (tsv format): [demo\_proteinfamilies\_annotation.taxonomy.tsv](https://github.com/biobakery/fugassem/raw/master/examples/raw_input/demo_proteinfamilies_annotation.taxonomy.tsv)
	* See the section on parallelization options to optimize the workflow run based on your computing resources. 
	* The workflow runs with the default settings for all main tool subtasks. If you need to customize your workflow settings for the preprocessing workflow to modify the default settings, you can change the parameter settings.
		* For example, `--extension-paired "$R1_suffix,$R2_suffix"`, `--extension "$fastq_suffix"` (what are the following part after `$SAMPLE` in the input file names) will modify the default settings when running the assembly task.
	
	#### Demo run of preparing MTX abundance utility
	
	```
	$ fugassem_generate_stratified_mtx_input --taxon-level Species \
	--gene-catalog input/demo_genecatalogs.clstr \ 
	--gene-catalog-seq raw_input/demo_genecatalogs.centroid.fna \ 
	--protein-family raw_input/demo_prshotgunoteinfamilies.clstr \ 
	--family-taxonomy raw_input/demo_proteinfamilies_annotation.taxonomy.tsv \ 
	--basename $BASENMAE --input raw_reads \
	--output $OUTPUT_DIR
	```
	
	#### Output files of preparing MTX abundance utility
	
	The main output file is a normalized MTX abundance tale:
	
	```
	ID  sample1 sample2
	Cluster_100569|k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Enterobacteriaceae.g__Escherichia.s__Escherichia_coli.t__Escherichia_coli_O157_H7.msp__msp_008   31585.1 23428.9 
	Cluster_1022788|k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Enterobacteriaceae.g__Escherichia.s__Escherichia_coli.msp__msp_168  0   440932  
	Cluster_1022791|k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Enterobacteriaceae.g__Escherichia.s__Escherichia_coli.msp__msp_168  0   434418  
	Cluster_1022793|k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Enterobacteriaceae.g__Escherichia.s__Escherichia_coli.msp__msp_168  0   124650  
	...

	```
	
	* File name: `$OUTPUT_DIR/$BASENAME.proteinfamilies.nrm.tsv`
	* This file includes the normalized MTX abundance stratified in each taxon (TSV format file).
	* `$OUTPUT_DIR` = the output folder
	* `$BASENAME` = the basename of output file
	* This file provides the normalized MTX abundance stratified in each taxon for each gene family (rows) in each sample (columns) (TSV format file). 
	

### Preparing evidence input

FUGAsseM provides another type of utilities to prepare evidence inputs used by FUGAsseM based on the outputs of [MetaWIBELE](https://github.com/biobakery/metawibele) and homology-based annotation, including the annotation file of MetaWIBELE, assembled information file, homologous information etc. This utility generates evidence files such as raw GO annotations, domain-based annotation evidence, assembly-based annotation evidence.

#### Generating workflow

`$ fugassem_generate_annotation_input --help`
	
This command yields:
	
```
usage: fugassem_generate_annotation_input [-h] [--version]
		                                          [--func-type {GO,BP,MF,CC,UniRef90_GO,UniRef90_GO_BP,UniRef90_GO_CC,UniRef90_GO_MF,UniRef90_COG,UniRef90_eggNOG,UniRef90_KEGG-KOs,InterProScan_PfamDomain,Denovo_transmembrane,Denovo_signaling,DOMINE_interaction}]
	                                          [--pfam PFAM] [--ddi DDI] [--contig]
	                                          [--clust-file CLUST_FILE]
	                                          [--gene-info GENE_INFO]
	                                          [--coann-list COANN_LIST]
	                                          [--homology]
	                                          [--homology-ann HOMOLOGY_ANN]
	                                          [--basename BASENAME]
	                                          [--threads THREADS]
	                                          [--memory MEMORY] [--time TIME]
	                                          [--output OUTPUT] [-i INPUT]
	                                          [--config CONFIG]
	                                          [--local-jobs JOBS]
	                                          [--grid-jobs GRID_JOBS]
	                                          [--grid GRID]
	                                          [--grid-partition GRID_PARTITION]
	                                          [--grid-benchmark {on,off}]
	                                          [--grid-options GRID_OPTIONS]
	                                          [--grid-environment GRID_ENVIRONMENT]
	                                          [--grid-scratch GRID_SCRATCH]
	                                          [--dry-run] [--skip-nothing]
	                                          [--quit-early]
	                                          [--until-task UNTIL_TASK]
	                                          [--exclude-task EXCLUDE_TASK]
	                                          [--target TARGET]
	                                          [--exclude-target EXCLUDE_TARGET]
	                                          [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}]
```

* `--input`: annotation file generated by MetaWIBELE](https://github.com/biobakery/metawibele)
* `--contig`: if specified, will extract contig annotation
* `--clust-file`: clustering file of protein families (extended-fasta format), which can be generated by [MetaWIBELE](https://github.com/biobakery/metawibele)
* `--gene-info`: gene info file including the contig information of each member in protein families
* `--homology`: If specified, will extract co-homology annotation
* `--homology-ann`: homology-based annotation file, e.g. co-UniRef50 clustering annotation
* `--output`: the output directory.

#### Input files of preparing evidence utility
* MetaWIBELE's annotationfile (TSV format): [demo\_proteinfamilies\_annotation.tsv](https://github.com/biobakery/fugassem/raw/master/examples/raw_input/demo_proteinfamilies_annotation.tsv)
* clustering file of protein families clustered by non-redundant gene catalogs (extended-fasta format): [demo_proteinfamilies.clstr](https://github.com/biobakery/fugassem/raw/master/examples/raw_input/demo_proteinfamilies.clstr)
* source assembled contigs of assembled genes file (TSV format): [demo\_gene\_info.tsv](https://github.com/biobakery/fugassem/raw/master/examples/raw_input/demo_gene_info.tsv)
* UniRef50-like clustering annotation file (TSV format): [demo\_map\_proteinfamilies.ident50.tsv] (https://github.com/biobakery/fugassem/raw/master/examples/raw_input/demo_map_proteinfamilies.ident50.tsv) 
* See the section on parallelization options to optimize the workflow run based on your computing resources. 
* The workflow runs with the default settings for all main tool subtasks. If you need to customize your workflow settings for the preprocessing workflow to modify the default settings, you can change the parameter settings.
	
#### Demo run of preparing evidence utility
	
```
$ fugassem_generate_annotation_input --clust-file raw_input/demo_proteinfamilies.clstr \ 
--contig --gene-info raw_input/demo_gene_info.tsv \ 
--homology --homology-ann raw_input/demo_map_proteinfamilies.ident50.tsv \ 
--basename demo_proteinfamilies \ 
--input raw_input/demo_proteinfamilies_annotation.tsv \ 
--output $OUTPUT_DIR
```
	
#### Output files of preparing evidence utility
	
The main output file is a normalized MTX abundance tale:
	
The following are the main output files of preparing evidence utility that are followed the same used for FUGAsseM:
	
```
$OUTPUT_DIR/$BASENMAE_proteinfamilies.GO.simple.tsv
$OUTPUT_DIR/$BASENMAE_proteinfamilies.DDI.simple.tsv
$OUTPUT_DIR/$BASENMAE_proteinfamilies.contig.simple.tsv
$OUTPUT_DIR/$BASENMAE_proteinfamilies.GO.homology.tsv
```
	
**1. demo\_proteinfamilies.GO.simple.tsv**
			
```
Cluster_100559  GO:0009058
Cluster_100559  GO:0016788
Cluster_1008935 GO:0005985
Cluster_1008935 GO:0005737
...
```
		
* This file provides original function (e.g. GO) annotations for each protein family that are used for training models (TSV format).
	
**2. demo\_proteinfamilies.DDI.simple.tsv**

```
Cluster_100559  PF00501:PF00109
Cluster_100559  PF00501:PF00975
Cluster_100559  PF00501:PF02801
Cluster_100559  PF00668:PF00109
...
```
	
* This file provides the domain-domain annotations of protein families that will be used to build DDI networks as one type of matrix evidence (TSV format). 

**3. demo\_proteinfamilies.contig.simple.tsv**

```
Cluster_7792    CSM7KORG_contig_k105_16193
Cluster_10199   CSM79HI7_contig_k105_24044
Cluster_10250   HSM7J4MW_contig_k105_17216
Cluster_12531   CSM79HI7_contig_k105_3486
...
```
	
* This file provides the contig annotations of protein families that will be used to build co-contig networks as one type of matrix evidence (TSV format). 

**5. demo\_proteinfamilies.GO.homology.tsv**

```
ID      GO:0000155__seqSimilarity       GO:0000166__seqSimilarity       GO:0000271__seqSimilarity       GO:0003674__seqSimilarity
Cluster_100559  0       1.0     0       1.0     0       0       0       1.0     0       0       0       0       0       0
Cluster_100569  0       0       0       0       0       0       0       0       0       0       0       0       0       0
Cluster_1008935 0       0       0       1.0     0       0       0       1.0     0       1.0     1.0     0       0       0
Cluster_101048  0       1.0     1.0     1.0     0       0       0       1.0     1.0     0       0       0       0       0
...
```
	
* This file provides the co-homology (e.g. co-uniref50) based annotations of protein families that will be used as one type of vector evidence (TSV format).

----



