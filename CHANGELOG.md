
# IRFinder Changelogs

**2.0.0**
 1. **Novelties**
    1. New **Long** RunMode to process fast[q|a] files from long reads
    using Minimap2 as aligner.
    2. New **-l** argument in **BAM** RunMode, to process long reads using an alternative algorithm. More information in the paper.
     3. New **AI** RunMode that uses a CNN model to detect false IR events on introns without warning in the last column of the result `IRFinder-IR-[non]dir.txt` file. It will generate a file containing only validated introns ( `IRFinder-IR-[non]dir-val.txt` )
    4. New **Diff** RunMode that uses SUPPA2 ( https://github.com/comprna/SUPPA ) or DESeq2 algorithm to identify differential IR events. 
    5. New **CLI** with dedicated helps for each RunMode and a verbose mode.
    6. New **installation script**, to check the dependencies and install or uninstall IRFinder globally and locally.

    7. **Docker** and **Singularity** images available, based on Ubuntu 18 LTS ( bionic ) and containing IRFinder and all his dependencies ( latest versions of STAR, Minimap2 and SUPPA2).

2. **Major changes** ( can impact the results between different versions ) 
    1. **NonUniformIntronCover** warning threshold simplified: now it uses the  25th/50th/75th percentile of intronic depth. Changed from:
``
(max(Column13, Column14) > 2 + Column9 && max(Column13, Column14) > Column9 * 1.5 ) || (min(Column13, Column14) + 2 < Column9 && min(Column13, Column14)*1.5 < Column9 ) 
``
to 
`` Column12-Column10 > Column11 ``
    2. **Default  Mapability read length** is now 100 instead of 70. It's not anymore hard coded and can be changed with the argument **-n** in the RunModes *BuildRef[Process|FromSTARRef]*
    3. **Paired reads with one pair unmapped**  are now processed as single reads instead of being removed.

3. **Minor changes** ( no impact on the results but improves the usability)
   1. The mapability file can be given as argument **-M** in the RunModes *BuildRef[Process|FromSTARRef]*. Precomputed mapabilities for hg38 are available under the git subdirectory */REF/Mapabilities/hg38/* for different read lengths ( 70, 100 and 150). This will reduce drastically the time to build the IRFinder reference.
   2. New argument **-l** in *BuildRefFromSTARRef* to create a **l**ink to the existing reference STAR folder, the genome file and the annotation file, instead of copy them. This will save disk space in case of multiple IRFinder reference directories using the same STAR reference.

**1.3.1**    
1. IRFinder now exits immediately after error, instead of trying to complete the remaining processes.    
2. Improved Perl version judgement during Phase 3 of reference preparation.    
    
**1.3.0**    
New features:    
1. New `BuildRefFromSTARRef` mode. This allows users to use an existing STAR reference to build IRFinder reference, which significantly reduces the total preparation time. This new mode also tries to automatically figure out the original FASTA and GTF files used to generate the existing STAR reference. Call `IRFinder -h` for more details.    
2. `BuildRef` and `BuildRefProcess` mode now support `-j` option to parse an integer that changes the default value of `--sjdbOverhang` argument in STAR.    
3. `FASTQ` mode now supports `-y` option to feed extra STAR arguments to control alignment behaviors.    
    
Improvements:    
1. `FASTQ` mode now outputs a full BAM file in "Unsorted.bam", instead of a BAM file with a trimmed QS column.   
2. IRFinder does not automatically generate "unsorted.frag.bam" to save disk space and to avoid redundancy to "Unsorted.bam". Instead, IRFinder now provides a tool at `bin/TrimBAM4IGV` to generate this kind of trimmed BAM file to facilitate visualization purpose in IGV.     
3. Re-design of standard output information during IRFinder reference preparation. It is easier to recognize occured errors now.    
4. Usage information now can be viewed by `-h` option.     
     
Bug fixes:    
1. The mapability calculation during the IRFinder reference preparation stage has been re-designed. The previous algorithm encountered buffer size issues when dealing with genomes with a huge amount of chromosomes/scaffolds. This has been fixed. Please note, the new algorithm requires `samtools` (>=1.4) executable binary ready in $PATH.    
2. Since Perl 5.28.0, `sort '_mergesort'` is no longer supported. IRFinder now checks the Perl version and uses `sort` functions correspondingly.    
    
**1.2.6**
1. IRFinder now keeps introns with the same effective regions as separate entries in the reference.    
2. IRFinder now automatically checks if the reference preparation stage generates empty reference files, which indicates process failure.    
3. The R object genreated by Differential IR Analysis script now includes an additional slot named "MaxSplice", which represents the maximum splice reads at either end of introns. Each value is the maximum value between Column 17 and 18 in the IR quantification output.    
4. During differential IR analysis, values in "MaxSplice" are now used as the denominators in the GLM, instead of using the values of Column 19 in the IR quantification output. This makes the IR ratio in the differential IR analysis more consistent with the values of Column 20 in the IR quantification output.    
5. User manual has been updated.    
    
**1.2.5**
1. Headers are now correctly added to output files `IRFinder-IR-dir.txt` and `IRFinder-IR-nondir.txt`.

**1.2.4**
1. In the GLM-based method for differential IR comparison, now the orginal matrix for DESeq2 is now made up by IR depth and correct splicing depth. In the previous versions, the latter one is the sum of splicing depth and IR depth. This change is supposed to give a smoother dispersion estimation across all introns.

**1.2.3:**
1. IRFinder now supports GTF attribution tags `gene_type` and `transcript_type` upon the original requirement for typical Ensembl tags `gene_biotype` and `transcript_biotype`. Either of these two pairs is required to correctly build IRFinder reference.    
    
**1.2.2:**
1. In GLM-based differential IR comparison, fixed an error caused by duplicated row names when creating DESeq2 object with a version of DESeq2 later than 1.10.

**1.2.1:**
1. Improved the performance of DESeq2-based GLM analysis for differential IR. This new approach should improve the estimation of dispersion. Normal splicing from IRFinder result is now used as a variable in the GLM, instead of using the value of normal splicing as an offset. This approach is adapted from [detection of allele-specific expression](http://rpubs.com/mikelove/ase) from Michael Love. See Wiki page for details.
2. Updated some out-of-date usage information

**1.2.0:**
1. IRFinder is now compatible with GLM-based analysis. This is achieved by passing IRFinder result to DESeq2 using the function in bin/DESeq2Constructor.R. See Wiki page for details  
2. Fixed the conflict with latest version "bedtools complement" that used to cause failure in preparing IRFinder reference  
3. Improved memory usage when passing lines to bedtools genomecov. This is also supposed to benefit reference preparation of those genomes with a lot of chromosomes contigs. Thanks for the smart solution from Andreas @andpet0101.  
4. Specified the gtf file to be downloaded during reference preparation via automatic downloading. Ensembl currently holds several versions of gtf files for the same genome release. This confused IRFinder BuildRefDownload function in the previous version.
5. Added -v option to print out version number.
