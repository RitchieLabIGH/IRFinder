# IRFinder
Detecting intron retention from RNA-Seq experiments

[User Manual](https://github.com/williamritchie/IRFinder/wiki)

[CHANGELOG](https://github.com/williamritchie/IRFinder/CHANGELOG.md)

## About IRFinder

IRFinder, developed at the [Center for Genomic Medicine of Massachusetts General Hospital](https://cgm.massgeneral.org/), the [CNRS](http://www.cnrs.com) and the [Centenary Institute](https://www.centenary.org.au), implements an end-to-end analysis of intron retention (IR) from mRNA sequencing data in multiple species.    
IRFinder includes alignment via the STAR algorithm, quality controls on the sample analyzed, IR detection, quantification and statistical comparison between multiple samples. IRFinder was capable of estimating IR events with low coverage or low mappability as confirmed by RT-qPCR.  
    
## Before Start: Intron Retention Database - [IRBase](http://mimirna.centenary.org.au/irfinder/database/)
Before diving into IRFinder package, users might also consider [IRBase](http://mimirna.centenary.org.au/irfinder/database/). It is a database for human IR inquiry and visualization, based upon pre-calculated IRFinder results from **over 2000** public available human RNA-Seq samples.     
[IRBase](http://mimirna.centenary.org.au/irfinder/database/) allows users to either 1) enquire, visualize and download single-gene IR results in a tissue/cell-type of interest, or 2) download transcriptome-wide IR results of a sample of interest.     
        
## Cite IRFinder    
Middleton R*, Gao D*, Thomas A, Singh B, Au A, Wong JJ, Bomane A, Cosson B, Eyras E, Rasko JE, Ritchie W. **IRFinder: assessing the impact of intron retention on mammalian gene expression**. [Genome Biol. 2017 Mar 15;18(1):51](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1184-4). doi: 10.1186/s13059-017-1184-4. [PubMed PMID: 28298237](https://www-ncbi-nlm-nih-gov.ezp-prod1.hul.harvard.edu/pubmed/28298237).

