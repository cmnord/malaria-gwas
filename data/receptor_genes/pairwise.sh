######## for doing pairwise alighments and putting all the results in the same file
########----------------------------------------------------------------------------
for file in duffy/outer/pairwise1/*
do muscle -in "$file" >> duffy/outer/pairwise1/pairwise1.afa
done

######## for running mulltiple alignment on multiple files
########-------------------------------------------------
#for file in duffy/outer/*
#do muscle -in "$file" -out ${file:0:(-3)}.afa
#done