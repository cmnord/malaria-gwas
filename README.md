# 6.047 Final Project

[Proposal][proposal]

[User guide][muscle] to MUSCLE with instructions on how to download. I used this to make the multiple alignments and phylogenetic trees.

[This website][muscleviz] is useful for visualizing the alignments (muscle just spits out another fasta).

[Tools][drawtree] to make the phylogenetic trees: I used the Drawtree Java app.

[Roadmap Epigenomics Project metadata][roadmap]

INFECTED PRIMATES:

-   Homo Sapiens (human)
-   Pan Troglodytes (chimp)
-   Gorilla ... duh
-   Aotus
-   Saimiri
-   Pongo abelii (orangutan) <-- I think one of them has another Pongo, another orangutan species

Non-human primates: download from [here][usc]:
**Bolded** primates are infected by malaria.

| Scientific name           | Common name              |
| ------------------------- | ------------------------ |
| **Aotus nancymaae**       | Nancy Ma's night monkey  |
| **Aotus trivirgatus**     | Night monkey             |
| Callithrix jacchus        | Common marmoset          |
| Cebus capuchinus imitator | Capuchin monkey          |
| Cercocebus atys           | Sooty mangabey           |
| Chlorocebus aethiops      | Grivet                   |
| Galveopterus variegatus   | Malayan flying lemur     |
| **Gorilla**               | Gorilla gorilla          |
| Macaca fascularis         | Crab eating macaque      |
| Macaca mulatta            | Rhesus macaque           |
| Mandrillus leucophaeus    | Drill                    |
| Nomascus leucogenys       | Gibbon                   |
| Pan paniscus              | Bonobo                   |
| **Pan troglodytes**       | Chimpanzee               |
| Papio anubis              | Baboon (olive)           |
| Pilicolorbus tephrosceles | Ugandan red colobus      |
| **Pongo alebii**          | Orangutan                |
| Rhinopithecus roxellana   | Golden snub-nosed monkey |
| **Saimiri boliviensis**   | Squirrel monkey          |
| ?                         | Bushbaby                 |
| ?                         | Denisova                 |
| ?                         | Green Monkey             |
| ?                         | Marmoset                 |
| ?                         | Proboscis monkey         |
| ?                         | Tarsier                  |

DATA - GWAS Analysis

-   We have the epigenetic data across all cell types
-   We have downloaded (?) the GWAS catalog, which has other SNPs that are found to be significant for other phenotypes

TOOLS - GWAS Analysis

-   `bedtools` - command line tool for searching through epigenomic data
    -   useful commands may include: `intersect`, `count-overlap` (documentation has good explanations)
-   R package `genomicRegions` - can be used to do the same thing and can make nice graphics (need to install using bioconductor)
-   Python

STEPS - GWAS Analysis

1.  DONE (Kari): Compile SNP gene locations from the Malaria GWAS study that we have
1.  DONE (Kari): Compile other SNPS from the GWAS catalog, which we use as our background
1.  DONE (Claire): Use `bedtools` to find the overlap of our Malaria SNPs with Roadmap epigenetic data
1.  DONE (Claire): Use `bedtools` to find overlap of background SNPs with Roadmap epigenetic data
1.  DONE (Kobbie): Use Python to develop figures displaying the overlap between our malaria SNPs and the Roadmap epigenetic data
1.  Calculate whether or not the difference is significant - p value, all that fun stats stuff
1.  DONE (Claire): Find overlap of malaria SNPs with RBC gene expression data.
1.  DONE (Claire): Find overlap of background SNPs with RBC gene expression data.
1.  DONE (Kobbie): Use Python to update/normailze relevant figures, displaying the overlap between our malaria SNPs and the Roadmap epigenetic data normailized by the overlap between all SNPs (background and malaria) and the Roadmap epigenetic data.
1.  Down the line, we can also use the DAVID database to look at pathways that are enriched <-- think about this later

[proposal]: https://docs.google.com/document/d/1F0Ke9Pjggio1-GSsk4dtYaaRajI1zjJQ_VCiW0mkeaQ/edit#
[usc]: http://hgdownload.cse.ucsc.edu/goldenPath/panPan2/bigZips/
[roadmap]: https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15
[muscle]: http://www.drive5.com/muscle/muscle.html
[muscleviz]: https://www.ebi.ac.uk/Tools/msa/mview/
[drawtree]: http://evolution.genetics.washington.edu/phylip.html
