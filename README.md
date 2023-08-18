# In house scripts used in the paper "Multi-Omics Integration Can be Used to Rescue Metabolic Information for Some of the Dark Region of the *Pseudomonas putida* Proteome"

## Within species Guilt by Association model

### GBA_integration
**SLT01_baselineDataSetup.py:** Parses the go.obo into a list of python objects and a networkx graph, collects initial annotation information and makes a list of protein objects.
**SLT01_combineEdgelist.py:** combines individual edgelist files for use by the first within-species model
**SLT01_makeEdgelist.py:** formats raw similarity data into edgelists
**SLT01_randomForestModel-final.py:** Trains and runs inference for the initial protein-protein similarity model
**SLT01_termCentric-inputData.py:** Formats the data for the term transfer model
**SLT01_termCentric-semiSupervisedModel.py:** Fits the semisupervised term transfer model
**test_set.txt:** The proteins held out as a test set for the prediction model

### coexpression
**names.txt:** A mapping between the internally used mass spec data filenames and the filenames from the online repositories
**valid_defects.txt:** A list of mass defects corresponding to plausably biological PTMs for filtering the ANNSoLo output
**SLT01_spectraST-commands.txt:** The specific commands used when processing spectral libraries with SpectraST
**SLT01_coEx-annsoloToFlashLFQ.py:** Processes the ANNSoLo output, does FDR control, and filters implausable PTMs then formats the output for quantificaiton by FlashLFQ
**SLT01_coEx-processDecoyLibrary.py:** Filters the NIST mouse spectral library of any peptides that are also potentially generated by *P. putida* 

### evolutionary_correlation
**SLT01_corrEvoFixedTopo-scheduler.py:** Run RAxML-ng on all orthogroups with the species tree as a topological constraint
**SLT01_corrEvoPargenesInput.py:** Cleans and formats the MSAs for Pargenes
**SLT01_corrEvoRenameMSAgenes.py:** Renames genes in the MSAs to their species for use in Pargenes
**SLT01_corrEvoTreeCMP-cleanNewick2.py:** Processes the Newick outputs of the topologically constrained RAxML-ng runs to work with TreeCMP
**SLT01_corrEvoTreeCMP-fullComparison.py:** Breaks the TreeCMP inputs into blocks and runs TreeCMP on those blocks in parallell 
**SLT01_MAFFTalign.py:** Runs MAFFT on all orthogroups

### structural_similarity
**SLT01_TMalignScheduler.py:** Runs TM-align on all *P. putida* proteins
**SLT01_trimAF.py:** Trimms low confidence regions at the termini of proteins as a preprocessing step for alignment

## Between Species Structural Similarity Model
### structural_similarity_model
**SLT01_RUPEEmapNamesAndcollectInputData.py:** Maps names from PDB hits to Uniprot IDs
**SLT01_RUPEEpairwiseAlignments.py:** Runs NWalign on the hits and collects pairwise similarity information
**SLT01_RUPEEsemisupervisedInput.py:** Collects the input information and formats it for the term transfer model
**SLT01_RUPEEsemisupervisedModel.py:** Fits the semi-supervised random forest model
**SLT01_RUPEEwebdriver.py:** automates searching *P. putida* proteins with the online RUPEE search tool

## Other scripts

**SLT01_paperFigures.py:** generates the underlying plots for all figures used in the paper. 

### GO_enrichment
**SLT01_GOenrichmentModel.stan:** The statistical model for the GO enrichment analysis in Figure 5
**SLT01_GOenrichmentModel-driver.py:** Formats data and runs the above model
**SLT01_GOexpectationAnalysis.stan:** The statistical model for the predicted term count model in Figure 6
**SLT01_GOexpectationAnalysis-driver.py:** Formats data and runs the above model

### interpro_analysis
**SLT01_InterproEnrichmentModel.stan:** The statistical model for assessing interpro term enrichments among PUFs. Used for figure S3.
**SLT01_InterproEnrichmentModel-driver.py:** Formats the input data and runs SLT01_InterproEnrichmentModel.stan