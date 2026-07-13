# IgCaller

<img src="IgCaller_logo.png" align="right" width="220">

#### Reconstructing immunoglobulin (IG) B-cell receptor and T-cell receptor (TCR) gene rearrangements and oncogenic translocations from genomic data in lymphoid neoplasms

IgCaller is a python program designed to fully characterize the immunoglobulin (IG) B-cell receptor and T-cell receptor (TCR) gene rearrangements and oncogenic translocations in lymphoid neoplasms. It was originally developed to reconstruct the IG gene from WGS data. More recently, it has been extended to characterize also the TCR and to work with WES, capture-based NGS, and amplicon-based NGS data. For more information, read the original [publication](https://rdcu.be/b5tsw) and the manual below.


### Requirements

* IgCaller is based on python3 and requires the following modules: statistics, regex (v2.5.33 and 2.5.109), argparse (v1.1), numpy (1.18.3 and 1.24.4), scipy (v1.5.2 and v1.8.0), biopython (v1.81 and v1.85), pandas (v1.0.3 and v2.3.1), seaborn (v0.11.2 and v0.13.2), and matplotlib (v3.0.3 and v3.10.5). Although providing the versions of the modules tested, we are not aware about any specific version requirement for running IgCaller. Other modules used by IgCaller but already included in base python are: subprocess, sys, os, itertools, operator, collections, gzip, pickle, and difflib.
* The only required non-python program is [samtools](http://www.htslib.org) (v1.16.1 and v1.22 have been tested).
* IgCaller has only been validated with BAM files obtained using [BWA-MEM](http://bio-bwa.sourceforge.net/) (v0.7.15 and v0.7.17 have been tested).

### Installation

Download and uncompress the ZIP file or clone the IgCaller repository:

```
git clone https://github.com/ferrannadeu/IgCaller
```

Starting at version 1.4-beta, run the following command to allow executable permissions:

```
chmod +x /path/to/IgCaller/IgCaller
```

...and consider adding it to your PATH:

```
export PATH=/path/to/IgCaller/:$PATH
```

### Running IgCaller

#### Basic command line:

Basic command line to analyze the IG of a tumor sample processed by WGS at 30x with paired normal sample available:
```
/path/to/IgCaller/IgCaller -I /path/to/IgCaller/IgCaller_reference_files/ -V hg19 -C ensembl -T /path/to/bams/tumor.bam -seq wgs -seqDepth low -N /path/to/bams/normal.bam -pm paired -R /path/to/reference/genome_hg19.fa -o /path/to/IgCaller/outputs/
```

Adjust the parameters for sequencing technique, sequencing coverage, gene to be analyzed, etc. as explained below:

#### Mandatory arguments:
*	inputsFolder (-I): Path to folder containing IgCaller reference files.
*	genomeVersion (-V): Reference genome version [hg19, hg38].
*	chromosomeAnnotation (-C): Chromosome annotation [ensembl = without 'chr' (i.e. 1); ucsc = with 'chr' (i.e. chr1)].
*	bamT (-T): Path to tumor BAM file.
* sequencing (-seq): Sequencing technique [whole-genome sequencing (wgs), whole-exome sequencing (wes), high-coverage capture NGS (capture), or amplicon-based NGS (amplicon)].
* sequencingDepth (-seqDepth): Sequencing depth [low (~30x), int (50-200x), high (>250x)].
*	bamN (-N): Path to normal BAM file, if available.
*	pairedMode (-pm): Tumor and normal paired status [paired/unpaired, default=None]. Need to be specified when bamN is specified. Paired = Normal BAM used in all analyses (V(D)J sequence reconstruction to handle SNPs, analysis of CSR, oncogenic rearrangements, and estimation of tumor purity). Unpaired = Normal BAM only used for the analysis of CSR, oncogenic rearrangements, and estimation of tumor purity.
*	refGenome (-R): Path to reference genome FASTA file. Not mandatory, but recommended, when specifying a normal BAM file using the argument 'bamN'. Mandatory when '-bamN' is not specified and when '-seq' is set to 'amplicon'.

#### Optional arguments:
###### Gene/receptor to be analyzed:
*	geneToAnalyze (-g): Gene/receptor to be analyzed [ig/tcr/both, default=ig].

###### Amplicon-based sequencing:
*	primer (-pr): Location of the primer in the V gene for amplicon-based data [leader/fr1/cdr1/fr2/cdr2/fr3/Vseq, default=None]. It is used to cut the sequence after the primer location (if leader, fr1, cdr1, fr2, cdr2, fr3) or sequence (if Vseq) prior annotation and identity calculation. Only applicable when IgBLAST is used for the annotation using the '-a' argument. Required if '-seq amplicon'. See parameter '--primerStringency' and 'primerFasta' for further tunning.
* primerFasta (-prf): Fasta file with primer sequences annealing to V genes that need to be removed from the sequence. Mandatory when '-pr/--primer' is 'Vseq'. Sequence name (or header line) in fasta file must start with '>locus', where locus is IGH, IGK, IGL, TRA, TRAB, TRAG, or TRAD to allow locus-specific mapping and elimination of primers. Only applicable when IgBLAST is used for the annotation using the '-a' argument (recommended).
*	primerStringency (-prs): Primer stringency [strict/permissive, default=strict]. Mandatory when '-pr/--primer' is 'leader', 'fr1', 'cdr1', 'fr2', 'cdr2', or 'fr3'. 'strict' = only annotate sequences with nucleotides found after the expected location of the primer. 'permissive' = annotate sequences even if they start later than the expected location of the primer. Only applicable when IgBLAST is used for the annotation using the '-a' argument (recommended).

###### Purity of the tumor sample:
*	tumorPurity (-p): Purity (or tumor cell contect) of the tumor sample, if known [0-1, default=1]. It is used to adjust the scores and some internal cutoffs during the analysis. If unknown, assumes 1.

###### Output path and options:
*	outputPath (-o): A folder inside this directory will be created with the output [default = current directory].
* keepMiniIgBams (-kmb): Should IgCaller keep (ie no remove) mini IG BAM files used in the analysis? [no/yes, default=no].
* reportReadNames (-rrn): Report read names associated with each rearrangement found [no/yes, default=no].

###### Annotation tool and database:
* annotateSeq (-a): Annotate sequence using IgCaller built-in annotation workflow or using IgBLAST [builtin/igblast, default=igblast].
* annotateSeqDB (-aa): Database of sequences to be used by IgBLAST [imgt/ogrdb, default=imgt]. IMGT = release 202430-2 (23 July 2024); OGRDB = release 2024-10-12.

###### Chronic lymphocytic leukemia (CLL)-specific annotations:
* subsetsAnnotation (-subsets): Should CLL stereotype subsets be annotated? Only applicable to productive IGH gene rearrangements. 'imgt' to annotate #2 and #8 using IMGT criteria; 'major' to annotate all major subsets based on the definitions reported in Agathangelidis et al Blood 2021 [no/imgt/major, default = no].
* R110annotation (-R110): Should the R110 mutation in IGLV3-21 be annotated? Annotation based on the definitions reported in Nadeu et al Blood 2021 [no/yes, default = no].

###### Acquired N-glycosylation sites (AGS) annotation:
* agsAnnotation (-ags): Should acquired N-glycosylation sites (AGS) be annotated? Only applicable to productive IGH gene rearrangements and requires '-a igblast'. AGS annotation based on the definitions reported in Tatterton et al Blood 2025 [no/yes, default = no].

###### IG/TCR reconstruction:
* highSensitivity (-hs): highSensitivity = no: runs faster by skipping low-confidence rearrangements. highSensitivity = yes: may run significantly slower in some samples since it tries to recover low-confidence rearrangements [no/yes, default = yes].
* minimumNumberOfNucleotidesSoft (-mnns): Minimum number of soft-clipped nucleotides to consider a read as split-read and number of nucleotides to be aligned on a J-V, J-D or D-V break in order to recover the read as a split-read covering the rearrangement [integer number higher than 0, default = 5].
* phaseReadsBasedOnMutations (-prbm): Phase reads based on mutations identified in the V gene. They will be used to reconstruct the sequence on the V gene but not considered to calculate the score and mapping quality [no/yes, default = yes].
* keepInsertSizeOnlyRearrangements (-kisor): Specify if IG/TCR rearrangements supported only by insert-size reads (i.e., no split-reads) should be kept [no/yes, default = no]. If kept, they are labeled with the tag 'Rearrangement without junction coverage' and the nucleotide sequences, functionality, identity to germline, and other features are not computed. If specified, these rearrangements are reported in the filtered output independently of the value specified in the '-rop' argument.
* collapseSequences (-cs): Collapse IG/TCR sequences if they are similar at nucleotide level. If so, the sequence of the rearrangement supported with the highest score is kept, but number of reads and read names are merged in the output [no/yes, default = yes].
* collapseSequencesSimilarity (-css): Cutoff of similarity between two IG/TCR sequences in order to be collapsed [0-1, default=0.98].
* shortReportedVseq (-shortV): Make V sequence start at (approx) FR1 [yes/no, default=yes].
*	mappingQuality (-mq): Mapping quality cut off to filter out reads for IG/TCR V(D)J reconstruction [default=0].
*	baseQuality (-bq): Base quality cut off for samtools mpileup for mutation analysis [default=13].
*	minDepth (-d): Depth cut off to consider a position [0-inf, default=1].
*	minAltDepth (-ad): Alt depth cut off to consider a nucleotide [0-inf, default=1].
*	vafCutoff (-vaf): VAF cut off to consider a mutation when working with phased reads and without phased reads, respectively [two numbers between 0-1 separated by comma, default=0.5,0.1].
*	vafCutoffNormal (-vafN): VAF cut off to consider a variant in the normal sample [0-1, default=0.20].
* scoreCutoff (-s): Minimum score supporting a gene rearrangement in order to be considered as high confidence [default='empty'; will consider 5 for seqDepth = low, 10 for int, and 15 for high].
* scoreCutoffFilter (-sf): Minimum score supporting a gene rearrangement in order to be kept during the analysis (intermediate filtering step to speed up the analysis) [default='empty'; will consider 2 for seqDepth = low, 5 for int, and 10 for high].
* scoreCutoffCSR (-scsr): Minimum score supporting a CSR rearrangement in order to be considered as high confidence [default='empty'; will consider 5 for seqDepth = low, 10 for int, and 15 for high].
*	reportOnlyProductive (-rop): Report only productive rearrangements in the 'output_filtered' file [no/oof/yes, defaul=no]. 'no' = report all rearrangements identified. 'oof' means out-of-frame and is used to filter out as unproductive only those rearrangements with an out-of-frame junction or, in other words, to report both productive rearrangements and in-frame rearrangements that are unproductive due to stop codons. 'yes' = report only productive rearrangements.

###### Oncogenic rearrangements (i.e., translocations, etc.):
*	runOncogenicRearrangements (-ror): Run the analysis of oncogenic IG/TCR rearrangements [yes/no, default=yes].
*	runOnlyOncogenicRearrangements (-roor): Run only the analysis of oncogenic IG/TCR rearrangements [yes/no, default=no].
*	mappingQualityOncoIg (-mqOnco): Mapping quality cut off to filter out reads when analyzing oncogenic IG/TCR rearrangements [default=0].
*	minNumberReadsTumorOncoIg (-mntonco): Minimum score supporting an oncogenic IG/TCR rearrangement in order to be annotated [default=5].
*	minNumberReadsTumorOncoIgPass (-mntoncoPass): Minimum score supporting an oncogenic IG/TCR rearrangement in order to be considered as high confidence [default='empty'; will consider 6 for seqDepth = low, 10 for int, and 15 for high].
*	vafOncoIgPass (-vafOnco): Minimum VAF of an oncogenic IG/TCR rearrangement in order to be considered as high confidence [range: 0-1; default=0.05].
*	maxNumberReadsNormalOncoIg (-mnnonco): Maximum number of reads supporting an oncogenic IG/TCR rearrangement in the normal sample in order to be considered as high confidence [default=2].
* maxNumberCountInPoN (-mncPoN): Maximum number of count in panel of normals (PoN) in order to be considered as high confidence [default=2].
* genesOncoIg (-gOnco): Genes to be annotated based on proximity to the oncogenic IG/TCR rearrangement breakpoints. All genes, only protein coding genes, or only canonical portein coding transcripts. For hg19 reference, protein_coding_canonical is not supported and it is treated as protein_coding. [all/protein_coding/protein_coding_canonical, default=protein_coding_canonical].
* customGenesOncoIg (-cgOnco): Comma-separated list of genes to be annotated with higher priority on the non-IG/TCR breakpoint if they are found within '-gOncoDist' of the breakpoint [default=''].
* genesOncoIgDistance (-gOncoDist): Maximum distance in base pairs from the non-IG/TCR breakpoint to the closest gene to annotate it [default=250000].
* customGenesOncoIgDistance (-cgOncoDist): For genes provided in -cgOnco, maximum distance in base pairs from the non-IG/TCR breakpoint to the closest gene to annotate it [default=500000].

###### Purity calculation:
* estimatePurity (-ep): Estimate purity [no/yes, default=yes]. If specified, the purity is estimated based on the IG/TCR rearrangements and reported in the output. It does not effect the scores and internal cutoffs.
* estimatePurityCoverage (-epc): Count gene rearrangements when estimating purity based on drop of coverage even if IgCaller has not called the rearrangement [no/yes/igh, default=no]. 'igh' means only applied to IGH locus.
* scoreCutoffPurity (-scp): Minimum score supporting a rearrangement in order to be considered during purity calculation [default='empty'; same as 'scoreCutoff'].
* plotPurityCoverage (-ppc): Plot coverage along J genes to visually evaluate coverage estimations [no/yes, default='yes'].

###### Samtools-related arguments:
*	pathToSamtools (-ptsam): Path to the directory where samtools is installed. No need to specify it if samtools is found in 'PATH' [default = 'empty', assuming it is in 'PATH'].
*	numThreads (-@): Maximum number of threads used for samtools [default=1].


#### Running time:
IgCaller only requires 1 CPU, and it usually takes around 5-10 minutes to characterize the complete IG/TCR of one tumor sample, including the analysis of oncogenic rearrangements and purity estimations. Execution time increases with sequencing coverage and number of potential rearrangements present in the sample.

#### Tested on:
IgCaller was tested on a MacBook Pro (macOS Mojave, Big Sur, Ventura), Ubuntu (16.04, 18.04, 22.04), and MareNostrum 4 (Barcelona Supercomputing Center, SUSE Linux Enterpirse Server 12 SP2 with python/3.6.1).

#### Demo dataset:
A demo dataset to test IgCaller is provided under the "Demo" folder.

### Outputs

IgCaller returns a set of tab-separated files:

*	tumor_sample_output_filtered.tsv: High confidence rearrangements passing the defined filters.
*	tumor_sample_output_IGH.tsv: File containing all IGH rearrangements [for IG analyses].
*	tumor_sample_output_IGK.tsv: File containing all IGK rearrangements [for IG analyses].
*	tumor_sample_output_IGL.tsv: File containing all IGL rearrangements [for IG analyses].
*	tumor_sample_output_class_switch.tsv: File containing all CSR rearrangements [for IG analyses].
*	tumor_sample_output_TRA.tsv: File containing all TRA rearrangements [for TCR analyses].
*	tumor_sample_output_TRB.tsv: File containing all TRB rearrangements [for TCR analyses].
*	tumor_sample_output_TRG.tsv: File containing all TRG rearrangements [for TCR analyses].
*	tumor_sample_output_TRD.tsv: File containing all TRD rearrangements [for TCR analyses].
*	tumor_sample_output_oncogenic_rearrangements.tsv: File containing all oncogenic IG/TCR rearrangements (translocations, deletions, inversions, and gains) identified genome-wide.
*	tumor_sample_output_purity.tsv: File containing the estimated purity of each gene considered during purity calculation.
* tumor_sample_output_purity_RawTumorCoverage.pdf: File containing the plots of the raw tumor coverage along the genes used to estimate the purity.
* tumor_sample_output_purity_NormalizedTumorCoverage.pdf: File containing the plots of the normalized tumor coverage along the genes used to estimate the purity. Coverage is normalized by mean coverage of the tumor sample in the region considered and by the coverage of the normal sample (i.e., baseline coverage distribution), if available. 

### Other notes

An R script to help the study of mutational signatures in CLL is available under the "Mutational_signature_analysis_in_CLL" folder. This script aims to determine the presence/absence of non-canonical AID mutations (signature 9) in CLL patients using an already defined catalogue of single nucleotide variants.

### Citation

If you use IgCaller, please cite:

Nadeu, F., Mas-de-les-Valls, R., Navarro, A. et al. IgCaller for reconstructing immunoglobulin gene rearrangements and oncogenic translocations from whole-genome sequencing in lymphoid neoplasms. Nature Communications 11, 3390 (2020). https://doi.org/10.1038/s41467-020-17095-7.

### Contact

Bugs, comments and improvements can be submitted as GitHub [issues](https://github.com/ferrannadeu/IgCaller/issues) or directly to *nadeu@recerca.clinic.cat*. If running into any bugs or issues, please share a reproducible example.

### Releases
* v2.0:
  * Added functionality to reconstruct the T-cell receptor (TCR). See [issue #10](https://github.com/ferrannadeu/IgCaller/issues/10) and argument -g for details.
  * Added the possibility to annotate the reconstructed IG/TCR sequences using either the built-in annotation scheme of IgCaller or IgBLAST using IMGT or OGRDB databases (see -a and -aa). We recommend using IgBLAST (default), especially when paired normal (i.e., germline) BAM file is not available.
  * Added the reconstruction of partial (J-D or D-V only) rearrangements.
  * Added compatibility with data generated using amplicon-based NGS approaches (i.e., primer-based PCR amplification of IG/TCR rearrangements). See arguments -seq, -pr, -prf, and -prs for further details.
  * Added a module to calculate tumor purity based on the IG/TCR gene rearrangements (see arguments -ep, -epc, -scp, -ppc, and -pm, as well as the new output files *output_purity*).
  * Added chronic lymphocytic leukemia (CLL)-specific annotations: annotation of CLL stereotyped subsets (see -subsets) and IGLV3-21 R110 mutation (see -R110). If annotated, a tag is added next to the V(D)J genes annotation both for subsets (i.e. [CLL#2]) and R110 (i.e. [R110]).
  * Added acquired N-glycosylation sites (AGS) annotation: annotation of AGS (see -ags). AGS are defined as N-X-T/S, where X is any amino acid except proline. If annotated, the motif(s) found and its location is reported as a tag next to the V(D)J genes annotation (i.e. [CDR-AGS (NTT:CDR3)]). A rearrangement is defined as CDR-located AGS (CDR-AGS) if any of the AGS found are located in any of the CDR regions. AGS that cross FR and CDR borders are classified as CDR. Contraily, the rearrangement is labeled as FR-located AGS (FR-AGS) if none of the AGS are located in a CDR region. The rearrangement is labeled as 'No-AGS' if no AGS are found.
  * Significant improvements on sensitivity and specificity for both IG/TCR gene rearrangements and oncogenic alterations (see -hs, -sf, and -mnns).
  * Improved phasing of reads along the V gene (see -prbm).
  * Added some optional arguments for flexibility (see -rop, -kisor, -cs, -css, -ror, and -roor)
  * Improved annotation of oncogenic alterations, including the identification of N-nucleotides and other improvements (see -vafOnco, -gOnco, -cgOnco, -gOncoDist, and -cgOncoDist for further details).
  * Extended panel of normals for the analysis of oncogenic rearrangements.
  * Added the possibility to report the read names of the reads associated with each specific rearrangement identified (see --rrn).
  * Other minor improvements and edits (virtually) in all functions, including some default arguments.

For information about previous releases see the [releases notes](Releases_notes.md).
