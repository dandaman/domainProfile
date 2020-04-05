## Examples
```bash
# using kmeans with k=6  
./domainProfile_kmeans.pl -n Example -o examples/kmeans/ -l -p -k 6 -r 0.8 -q 0.05 -a 0.9 -s examples/Example.fasta -f fasta examples/Example.meme_fimo.jalview  examples/Example.PfamScan.jalview
./domainProfile_hclust.pl -n ward -o examples/hclust/ -l -p -m ward -r 0.8 -q 0.05 -a 0.9 -s examples/Example.fasta -f fasta examples/Example.meme_fimo.jalview examples/Example.PfamScan.jalview
./domainProfile_hclust.pl -n ward -o examples/hclust/ -d 2 -r 0.8 -q 0.05 -a 0.9 -l -p -m ward -s examples/Example.fasta -f fasta examples/Example.meme_fimo.jalview examples/Example.PfamScan.jalview



```


## Parameters
```
    -h # for help
    -o outpath 
    -n profile_name 
    -k kmeans k
    -l # optional use domain length for penalty 
    -p # optional report pairs 
    -m hclust_method [ward|single|complete|average|mcquitty|median|centroid] #optional otherwise ward
    -r ratio_considered_as_representative 
    -q qvalue_cutoff_for_outliers # optional #0.05 
    -a qvalue_cutoff_for_absent # optional #0.9
    -s seq_file
    -f format(fasta) #optional 
    -i ignore_motif_regex #optional 
    -d static_cut_at_tree_height #optional 
    JALVIEW_ANNOT_FILE(s)
```
