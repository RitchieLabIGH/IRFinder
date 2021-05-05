
# IRFinder
Detecting intron retention from RNA-Seq experiments
To start using IRFinder, read our [wiki user manual.](https://github.com/RitchieLabIGH/IRFinder/wiki)

[CHANGELOG](https://github.com/RitchieLabIGH/IRFinder/CHANGELOG.md)
IRFinder Version 1 is still available at https://github.com/williamritchie/IRFinder but is not anymore maintained.
## About IRFinder

IRFinder, developed at the [Center for Genomic Medicine of Massachusetts General Hospital](https://cgm.massgeneral.org/), the [CNRS](http://www.cnrs.com) and the [Centenary Institute](https://www.centenary.org.au), implements an end-to-end analysis of intron retention (IR) from mRNA sequencing data in multiple species.    
IRFinder includes alignment via the STAR (for short reads) and minimap2 (for long read) algorithm, quality controls on the sample analyzed, IR detection, quantification, convolutional neural network based validation and statistical comparison between multiple samples. 
IRFinder was capable of estimating IR events with low coverage or low mappability as confirmed by RT-qPCR.


    
## Before Start: Intron Retention Database - [IRBase](http://irbase.igh.cnrs.fr/)
Before diving into IRFinder package, users might also consider [IRBase](http://irbase.igh.cnrs.fr/). It is a database for human IR inquiry and visualization, based upon pre-calculated IRFinder results from **over 935** public available human cell lines RNA-Seq sample.     
[IRBase](http://irbase.igh.cnrs.fr/) allows users to enquire, visualize and download single-gene IR results in a tissue/cell-type of interest, download transcriptome-wide IR results of a sample of interest, upload your results to compare with the public ones and share them with the community.

        
## Cite IRFinder    
Middleton R*, Gao D*, Thomas A, Singh B, Au A, Wong JJ, Bomane A, Cosson B, Eyras E, Rasko JE, Ritchie W. **IRFinder: assessing the impact of intron retention on mammalian gene expression**. [Genome Biol. 2017 Mar 15;18(1):51](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1184-4). doi: 10.1186/s13059-017-1184-4. [PubMed PMID: 28298237](https://www-ncbi-nlm-nih-gov.ezp-prod1.hul.harvard.edu/pubmed/28298237).

