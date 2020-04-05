# DomainProfile
Author: Daniel Lang

The two perl scripts perform the following tasks:

1) Needleman-Wunsch Alignment of sequence domain profiles 
2) Use the resulting pairwise scores for clustering (kmeans/hierarchical)
3) plot the resulting clusters

## Install
Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html).
```bash
conda env create -f domainProfile.yml
```


## Examples
```bash
# using kmeans with k=6  
./domainProfile_kmeans.pl -n Example -o examples/kmeans/ -l -p -k 6 -r 0.8 -q 0.05 -a 0.9 -s examples/Example.fasta -f fasta examples/Example.meme_fimo.jalview  examples/Example.PfamScan.jalview
./domainProfile_hclust.pl -n ward -o examples/hclust/ -l -p -m ward -r 0.8 -q 0.05 -a 0.9 -s examples/Example.fasta -f fasta examples/Example.meme_fimo.jalview examples/Example.PfamScan.jalview
./domainProfile_hclust.pl -n ward -o examples/hclust/ -d 2 -r 0.8 -q 0.05 -a 0.9 -l -p -m ward -s examples/Example.fasta -f fasta examples/Example.meme_fimo.jalview examples/Example.PfamScan.jalview
```

## Parameters
```
######help
    -h # for help
######parameters for execution
    -o outpath 
    -s seq_file # sequence file
    -n profile_name 
    -k kmeans_k
    -l # optional use domain length for penalty 
    -p # optional report pairs 
    -m hclust_method [ward|single|complete|average|mcquitty|median|centroid] #optional otherwise ward
    -r ratio_considered_as_representative # optional default: 0.9
    -q qvalue_cutoff_for_outliers # optional default:0.05 
    -a qvalue_cutoff_for_absent # optional default:0.9
    -f format #optional default: fasta [Bioperl Bio::SeqIO formats] 
    -d static_cut_at_tree_height # optional 
    -i ignore_motif_regex # optional 
```
