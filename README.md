# IgCaller

<img src="IgCaller_logo.png" align="right" width="220">

#### Reconstructing immunoglobulin (IG) B-cell receptor and T-cell receptor (TCR) gene rearrangements and oncogenic translocations from genomic data in lymphoid neoplasms

IgCaller is a python program designed to fully characterize the immunoglobulin (IG) B-cell receptor and T-cell receptor (TCR) gene rearrangements and oncogenic translocations in lymphoid neoplasms. It was originally developed to reconstruct the IG gene from WGS data. More recently, it has been extended to characterize also the TCR and to work with WES, capture-based NGS, and amplicon-based NGS data. For more information, read the original [publication](https://rdcu.be/b5tsw) and the manual below.


### Requirements

* IgCaller is based on python3 and requires the following modules: statistics, regex (v2.5.29 and v2.5.30), argparse (v1.1), numpy (1.16.2, v1.16.3, and 1.21.6), scipy (v1.2.1, v1.3.0 and 1.7.3), and biopython (v1.81). Although providing the versions of the modules tested, we are not aware about any specific version requirement for running IgCaller. Other modules used by IgCaller but already included in base python are: subprocess, sys, os, itertools, operator, collections, gzip, pickle, difflib.

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
*	pairedMode (-pm): Tumor and normal paired status [paired/unpaired, default=None]. Need to be specified when bamN is specified. Paired = Normal BAM used for all analyses. Unpaired = Normal BAM only used to calculate tumor purity estimates based on coverage.
*	refGenome (-R): Path to reference genome FASTA file. Not mandatory, but recommended, when specifying a normal BAM file using the argument 'bamN'. Mandatory when '-bamN' is not specified and when '-seq' is set to 'amplicon'.

#### Optional arguments:
###### Gene/receptor to be analyzed:
*	geneToAnalyze (-g): Gene/receptor to be analyzed [ig/tcr/both, default=ig].

###### Amplicon-based sequencing:
*	primer (-pr): Location of the primer in the V gene for amplicon-based data [leader/fr1/cdr1/fr2/cdr2/fr3, default=None]. It is used to calculate the percentage of identity starting with the sequence after the primer location. Only applicable when IgBLAST is used for the annotation using the '-a' argument (recommended). Required if '-seq amplicon'. See parameter '--primerStringency' for further tunning.
*	primerStringency (-prs): Primer stringency [strict/permissive, default=strict]: strict = only annotate sequences with nucleotides found after the theoric location of the primer. permissive = annotate sequences even if they start later than the theoric location of the primer. Only applicable when IgBLAST is used for the annotation using the '-a' argument (recommended).

###### Purity of the tumor sample:
*	tumorPurity (-p): Purity (or tumor cell contect) of the tumor sample, if known [0-1, default=1]. It is used to adjust the scores and some internal cutoffs during the analysis. If unknown, use 1.

###### Output path and options:
*	outputPath (-o): A folder inside this directory will be created with the output [default = current directory].
* keepMiniIgBams (-kmb): Should IgCaller keep (ie no remove) mini IG BAM files used in the analysis? [no/yes, default=no].
* reportReadNames (-rrn): Report read names associated with each rearrangement found [no/yes, default=no].

###### Annotation tool and database:
* annotateSeq (-a): Annotate sequence using IgCaller built-in annotation workflow or using IgBLAST [builtin/igblast, default=igblast].
* annotateSeqDB (-aa): Database of sequences to be used by IgBLAST [imgt/ogrdb, default=imgt]. IMGT = release 202430-2 (23 July 2024); OGRDB = release 2024-10-12.

###### CLL-specific annotations:
* subsetsAnnotation (-subsets): Should CLL stereotype subset be annotated? 'imgt' to annotate #2 and #8 using IMGT criteria; 'major' to annotate all major subsets using Agathangelidis et al Blood 2021 definitions [no/imgt/major, default = no].
* R110annotation (-R110): Should the R110 mutation in IGLV3-21 be annotated? [no/yes, default = no].

###### IG/TCR reconstruction:
* highSensitivity (-hs): highSensitivity = no: runs faster by skipping low-confidence rearrangements. highSensitivity = yes: may run significantly slower in some samples since it tries to recover low-confidence rearrangements [no/yes, default = yes].
* keepInsertSizeOnlyRearrangements (-kisor): Specify if rearrangements supported only by insert-size reads (i.e., no split-reads) should be kept. keepInsertSizeOnlyRearrangements = yes: may run significantly slower [no/yes, default = no].
* phaseReadsBasedOnMutations (-prbm): Phase reads based on mutations identified in the V gene. They will be used to reconstruct the sequence on the V gene but not considered to calculate the score and mapping quality [no/yes, default = yes].
* minimumNumberOfNucleotidesSoft (-mnns): Minimum number of soft-cliped nucleotides to consider a read as split-read and needed to be align on a J/V break in order to recover the read as split [integer number higher than 0, default = 5].
* shortReportedVseq (-shortV): Make V sequence start at (approx) FR1 [yes/no, default=yes].
*	mappingQuality (-mq): Mapping quality cut off to filter out reads for IG/TCR V(D)J reconstruction [default=0].
*	baseQuality (-bq): Base quality cut off for samtools mpileup for mutation analysis [default=13].
*	minDepth (-d): Depth cut off to consider a position [0-inf, default=2].
*	minAltDepth (-ad): Alt depth cut off to consider a nucleotide [0-inf, default='empty'; will consider 1 for seqDepth = low/int and 2 for high].
*	vafCutoff (-vaf): VAF cut off to consider a mutation when working with phased reads and without phased reads, respectively [two numbers between 0-1 separated by a comma, default=0.66,0.1].
*	vafCutoffNormal (-vafN): VAF cut off to consider a variant in the normal sample [0-1, default=0.20].
*	reportOnlyProductive (-rop): Report only productive rearrangements in the 'output_filtered' file [yes/no, defaul=no].
* scoreCutoff (-s): Minimum score supporting a gene rearrangement in order to be considered as high confidence [default='empty'; will consider 5 for seqDepth = low, 10 for int, and 15 for high].
* scoreCutoffFilter (-sf): Minimum score supporting a gene rearrangement in order to be kept during the analysis [default='empty'; will consider 2 for seqDepth = low, 4 for int, and 6 for high].
* scoreCutoffCSR (-scsr): Minimum score supporting a CSR rearrangement in order to be considered as high confidence [default='empty'; will consider 5 for seqDepth = low, 10 for int, and 15 for high].

###### Oncogenic rearrangements (i.e., translocations, etc.):
*	runOncogenicRearrangements (-ror): Run the analysis of oncogenic IG/TCR rearrangements [yes/no, default=yes].
*	runOnlyOncogenicRearrangements (-roor): Run only the analysis of oncogenic IG/TCR rearrangements [yes/no, default=no].
*	mappingQualityOncoIg (-mqOnco): Mapping quality cut off to filter out reads when analyzing oncogenic IG/TCR rearrangements [default=15].
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

###### Samtools-related arguments:
*	pathToSamtools (-ptsam): Path to the directory where samtools is installed. No need to specify it if samtools is found in 'PATH' [default = 'empty', assuming it is in 'PATH'].
*	numThreads (-@): Maximum number of threads used for samtools [default=1].


#### Running time:
IgCaller only requires 1 CPU, and it usually takes around 5 minutes to characterize the complete IG/TCR of one tumor sample. Execution time increases with sequencing coverage and number of potential rearrangements present in the sample.

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
  * Added compatibility with data generated using amplicon-based NGS approaches (i.e., primer-based PCR amplification of rearrangements). See arguments -seq, -pr, and -prs for further details.
  * Added the reconstruction of partial (J-D only) rearrangements.
  * Added a module to calculate tumor purity based on the IG/TCR gene rearrangements (see arguments -pp, -ppc, and -sp).
  * Added the possibility to annotate the reconstructed IG/TCR sequences using either the built-in annotation scheme of IgCaller or IgBLAST using IMGT or OGRDB databases (see arguments -a and -aa). We recommend using IgBLAST, especially when the normal (i.e., germline) BAM file is not available.
  * Added chronic lymphocytic leukemia (CLL)-specific annotations: annotation of CLL stereotyped subsets (see -subsets) and IGLV3-21 R110 mutation (see -R110)
  * Improved annotation of oncogenic alterations, including the identification of N-nucleotides. See arguments -vafOnco, -gOnco, -cgOnco, -gOncoDist, and -cgOncoDist for further details.
  * Significant improvements on sensitivity and specificity for both IG/TCR gene rearrangements and oncogenic alterations.
  * Improved phasing of reads along the V gene.
  * Added the possibility to report the read names of the reads associated with each specific rearrangement identified.
  * Extended panel of normals for the analysis of oncogenic rearrangements.
  * Other minor improvements and edits.
* v1.4-beta (not stored as a formal release):
  * Dockerfile and shebang included following [pull #6](https://github.com/ferrannadeu/IgCaller/pull/6).
  * Small edit to correct a warning message when trying to remove non-created temporary files in some situations.
* v1.3:
  * Improved annotation of the IG genes involved in the translocations identified.
  * Improved selection of high-confidence rearrangements.
  * From v1.2.2:
    * Added phasing of reads based on exact breakpoints/split reads prior to consider phasing reads just if they span the given J-V pair. It allows the identification of multiple rearrangements within the same J and V genes if their breakpoints are different. 
    * Improved selection of high-confidence rearrangements after additional testing.
  * From v1.2.1:
    * Improved phasing of mutations within the rearranged allele. When reads spanning the J-V rearrangement are not covering a fraction of the V gene, IgCaller now tries to phase reads on the fly based on the previous mutations identified in the V gene. This new functionality improves the capacity to phase mutations and SNPs while increases the sensitivity to detect mutations specially in heavily mutated V genes. If this phasing of reads based on mutations does not recover any phased read covering a given position, all reads covering this specific position are then considered.
    * Phasing information is now considered when selecting high-confidence rearrangements to be included in 'tumor_sample_output_filtered.tsv'.
    * Corrected internal split-read annotation that made IgCaller crash in a subset of samples aligned to a reference genome that included "alt" contigs. See [issue #5](https://github.com/ferrannadeu/IgCaller/issues/5) for details. 
* v1.2:
  * General improvements (major):
    * Added compatibility for reference genomes containing lowercase nucleotides (lowercase nucleotides are automatically converted to uppercase in the fly).
    * Better tolerance of ambiguous nucleotides in the sequences.
    * Small imperfections when counting the nucleotides found in each genomic position during the reconstruction of the consensus V(D)J sequence have been solved/improved.
    * Increased sensitivity to detect V(D)J rearrangements by expanding the J gene region used to consider ‘insert-size’ reads spanning a given J-V rearrangement. Due to the small size of the J gene and insert size of the NGS libraries, some reads spanning a J-V rearrangement might align 5’ of the J gene. These reads are now recovered. Besides, shorter split reads are kept and tried to ‘re-align’ than in the previous version.
    * Improved (and simplified) pre-defined filtering of high-quality gene rearrangements. High-quality rearrangements found in the “output_filtered.tsv” file are now considered based on the score of the rearrangement, its mapping quality (new), mechanisms [deletion, inversion] (new), CDR3 sequence (new), and V(D)J genes annotated. As in the original version, it assumes a clonal/oligoclonal situation where the same CDR3 sequence and/or combination of V(D)J genes is not possible/real but caused by an inaccurate mapping of the reads leading to multiple solutions.
  * Phasing of mutations/SNPs within the V(D)J rearranged sequence:
    * IgCaller now considers only those reads spanning the J-V rearrangement of interest to better define mutations and SNPs present in the rearranged sequence/allele. If less than ‘minDepth’ reads spanning the J-V rearrangement are found in a given position, this position is analyzed considering all reads covering it, as performed in the previous versions. Note that a column with the details of the “phasing” step has been added (column labeled “V-PhaseInfo”). This column looks like “a/b - c”, where a,b,c are numbers: a = number of mutations phased, b = total number of mutations found in V gene, c = number of SNPs not phased.
    * 'vafCutoff / -vaf' default set to 0.66,0.1 for phased and unphased mutations/SNPs, respectively. This means that a minimum VAF of 0.66 is needed to consider a mutation when a position is analyzed using only “phased reads” (i.e. those spanning the J-V rearrangement), while a VAF of 0.1 is needed if a position is analyzed using all reads covering the position (i.e. “unphased reads”). 
    * Added ‘--shortReportedVseq / -shortV’ to make the V sequence start approx. at FR1 (active by default). Useful when interested only on the FR1-FR3 sequence so the reported number of phased mutations and SNPs is related to this shorter V gene sequence of interest.
    * The coordiantes of some V genes have been manually curated for a proper functionality of the “-shortV” method.
  * Genome-wide oncogenic IG rearrangements:
    * Incoherent filtering of split-read support when analyzing genome-wide IG translocations has been corrected (see issue #4: https://github.com/ferrannadeu/IgCaller/issues/4). It increases the sensitivity to detect IG oncogenic events (translocations, deletions, etc.).
    * The IGL locus considered when analyzing IG translocations has been slightly expanded.
  * General improvements (minor):
    * Added ‘-seq capture’ mode: added to improve compatibility of IgCaller with high coverage, capture-based NGS data (depth of >250/500x). By default, a minimum score of 15 is needed to keep a gene rearrangement as “high quality” in the “output_filtered.tsv”. Besides, to speed up the execution, under this ‘-seq capture’ mode, IgCaller removes in the first steps of the analysis those J-V rearrangements supported by less than two split reads (note that these rearrangements will not be included neither in the locus-specific output file). If the latter functionality is not desired, IgCaller could be run as “-seq wes” to analyze high coverage data, although it might remarkably increase the execution time and number of low-evidence (likely artefactual) rearrangements passing the default parameters.
    * Added '--scoreCutoff / -s': minimum score supporting a gene rearrangement in order to be considered as high confidence (by default will consider 3 for wgs/wes and 15 for capture). It might be useful when the default parameters are to permissive/restrictive. It should be adjusted based on the quality of the data and sequencing coverage.
    * ‘--minDepth / -d’ default set to 2 (it was set to 1 in the previous versions). A minDepth of 2 is necessary for a proper analysis considering “phased reads”.
    * Removed version number from the name of the scripts. It has been added within the scripts and accessible using ‘python3 IgCaller.py -v’.
    * Minor bug that made IgCaller crash in a specific scenario when annotating the gene near the breakpoints of the translocations is now fixed.
    * Other minor improvements related to strand orientation, filter of artefactual/really short V gene rearrangements, etc. have been implemented.
* v1.1: 
  * Added a panel of normals (PoN) with the oncogenic IG rearrangements found by IgCaller in normal WGS data, which can be considered as sequencing artifacts and filtered out from the tumoral samples. When filtering using these PoN, the breakpoint in the IG locus is not considered based on its position but considering the IG locus as a whole (i.e. if the break occurs in the IGH, IGK or IGL locus irrespectively of the exact position), while the exact position of the non-IG breakpoint is considered (with a +/- 1,000 bp window). The normal samples analyzed in the IgCaller manuscript were used to build a PoN for hg19 (n=243 samples) and hg38 (n=161). See optional argument --maxNumberCountInPoN (-mncPoN) to adjust this filter.
  * Added the annotation whether the non-IG breakpoint of the IG rearrangements identified map within interspersed repeats and low complexity sequences according to the [RepeatMasker UCSC track](https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=rmsk). repName, repClass and	repFamily are annotated.
  * Added the annotation of the gene closest to the non-IG breakpoint of the oncogenic IG rearrangements identified based on RefSeq annotations (maximum upstream and downstream distance considered: 250 kb).
  * Minor bug that made IgCaller crash in 1 out of the >1,000 samples tested is now fixed.
  * Added compatibility for BAM files obtained from whole-exome sequencing (see optional argument --sequencing/-seq). The functionality of IgCaller is the same for WGS and WES samples with only a small difference in the pre-defined filtering step.  
* v1.0:
  * First version of IgCaller as described in the original manuscript.
