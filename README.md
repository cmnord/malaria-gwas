# 6.047 Final Project

[Proposal](https://docs.google.com/document/d/1F0Ke9Pjggio1-GSsk4dtYaaRajI1zjJQ_VCiW0mkeaQ/edit#)


Non-human primates: download from [here][usc]:
- Baboon
- Bonobo
- Bushbaby
- Chimpanzee
- Crab eating macaque
- Denisova
- Gibbon
- Golden snub-nosed monkey
- Gorilla
- Green Monkey
- Marmoset
- Orangutan
- Proboscis monkey
- Rhesus macaque
- Squirrel monkey
- Tarsier

DATA - GWAS Analysis
-  We have the epigenetic data across all cell types
-  We have downloaded (?) the GWAS catalog, which has other SNPs that are found to be significant for other phenotypes

TOOLS - GWAS Analysis
-  `bedtools` - command line tool for searching through epigenomic data
    - useful commands may include: `intersect`, `count-overlap` (documentation has good explanations)
-  R package `genomicRegions` - can be used to do the same thing and can make nice graphics (need to install using bioconductor)

STEPS - GWAS Analysis
1.  DONE (Kari): Compile SNP gene locations from the Malaria GWAS study that we have
1.  DOING (Kari): Compile other SNPS from the GWAS catalog, which we use as our background
1.  DONE (Claire): Use `bedtools` or R to find the overlap of our Malaria SNPs with Roadmap epigenetic data
1.  TODO (Claire): Use `bedtools` or R to find overlap of background SNPs with Roadmap epigenetic data
1.  Calculate whether or not the difference is significant - p value, all that fun stats stuff
1.  Repeat step 3-5 but with the RBC gene expression data (first need to find out if the data we have is the raw expression data, or if it shows genes that are expressed only in RBCs. This part we use with gene names, rather than locations, because of the nature of the data, so may be best done in R)
1.  Down the line, we can also use the DAVID database to look at pathways that are enriched <-- think about this later

[usc]: http://hgdownload.cse.ucsc.edu/goldenPath/panPan2/bigZips/


