######## for running multiple alignment on multiple files
######## runs python script to compare each sequence with 
########   human sequence
########-------------------------------------------------
for file in tfr1/outer/*;
do 
n=${file:0:(-3)}.$afa
echo ${n}
muscle -in "$file" -out ${n}.$afa
#python2 pairwise_similarities.py ${n}
done


# ######## get the pairwise fastas, to put in a pairwise folder
# ########-----------------------------------------------------
# mkdir duffy/outer/pairwise1
# python2 make_pairwise_fastas.py "duffy/outer/duffy_all_primates_outer1.fa"

# ######## for doing pairwise alighments 
# ########----------------------------------------------------------------------------

# for file in duffy/outer/pairwise1/*
# do muscle -in "$file" >> duffy/outer/pairwise1/pairwise1.afa
# done
# python2 pairwise_similarities.py "duffy/outer/pairwise1/pairwise1.afa"
