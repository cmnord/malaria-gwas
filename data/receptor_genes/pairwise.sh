######## for running mulltiple alignment on multiple files
########-------------------------------------------------
for file in duffy/outer/*
do muscle -in "$file" -out ${file:0:(-3)}.afa
done

######## get the pairwise fastas, to put in a pairwise folder
########-----------------------------------------------------
mkdir duffy/outer/pairwise1
python2 make_pairwise_fastas.py "duffy/outer/duffy_all_primates_outer1.fa"

######## for doing pairwise alighments
########----------------------------------------------------------------------------

for file in duffy/outer/pairwise1/*
do muscle -in "$file" >> duffy/outer/pairwise1/pairwise1.afa
done
python2 pairwise_similarities.py "duffy/outer/pairwise1/pairwise1.afa"
