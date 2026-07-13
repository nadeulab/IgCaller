# Releases notes

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
