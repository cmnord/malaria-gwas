# for doing pairwise alighments and putting all the results in the same file
#for file in tfr1/pairwise/*
#do muscle -in "$file" >> tfr1/pairwise/pairwise.afa
#done

# for running mulltiple alignment on multiple files
for file in duffy/outer/*
do muscle -in "$file" -out ${file:0:(-3)}.afa
done