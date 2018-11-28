# 6.047 Final Project

[Proposal][proposal]

[Roadmap Epigenomics Project metadata][roadmap]

Non-human primates: download from [here][usc]:

-   Baboon
-   Bonobo
-   Bushbaby
-   Chimpanzee
-   Crab eating macaque
-   Denisova
-   Gibbon
-   Golden snub-nosed monkey
-   Gorilla
-   Green Monkey
-   Marmoset
-   Orangutan
-   Proboscis monkey
-   Rhesus macaque
-   Squirrel monkey
-   Tarsier

DATA - GWAS Analysis

-   We have the epigenetic data across all cell types
-   We have downloaded (?) the GWAS catalog, which has other SNPs that are found to be significant for other phenotypes

TOOLS - GWAS Analysis

-   `bedtools` - command line tool for searching through epigenomic data
    -   useful commands may include: `intersect`, `count-overlap` (documentation has good explanations)
-   R package `genomicRegions` - can be used to do the same thing and can make nice graphics (need to install using bioconductor)

STEPS - GWAS Analysis

1.  DONE (Kari): Compile SNP gene locations from the Malaria GWAS study that we have
1.  DONE (Kari): Compile other SNPS from the GWAS catalog, which we use as our background
1.  DONE (Claire): Use `bedtools` to find the overlap of our Malaria SNPs with Roadmap epigenetic data
1.  DOING (Claire): Use `bedtools` to find overlap of background SNPs with Roadmap epigenetic data
1.  DOING (Kobbie): Use python of R to develop figures of displaying the overlap of our Malaria SNPs with Roadmap epigenetic data, and the overlap of background SNPs with Roadmap epigenetic data
1.  Calculate whether or not the difference is significant - p value, all that fun stats stuff
1.  DONE (Claire): Find overlap of malaria SNPs with RBC gene expression data.
1.  TODO (Claire): Find overlap of background SNPs with RBC gene expression data.
1.  TODO (Kobbie): Make figures displaying overlap
1.  Down the line, we can also use the DAVID database to look at pathways that are enriched <-- think about this later

[proposal]: https://docs.google.com/document/d/1F0Ke9Pjggio1-GSsk4dtYaaRajI1zjJQ_VCiW0mkeaQ/edit#
[usc]: http://hgdownload.cse.ucsc.edu/goldenPath/panPan2/bigZips/
[roadmap]: https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15
