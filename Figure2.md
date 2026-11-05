# eLRR tree, domain architecture, and models for Figure 2

## eLRR Tree - Figure 2B

### All sequences and put into alphfold
All sequences were gathered and inputted into the 'Antagonism_review' excel sheet. Sequences were loaded into Geneious and interproscan was used for domain prediction. Mature proteins were then loaded into Alphafold 3 for modeling. In the case of AtBAK1 and AtFLS2, the 5LRRs adjacent to the membrane were used for modeling.

To select the best alpha fold models the script 'select_best_models.py' was written with help from Claude. 
To run the script:
```
python select_best_models.py --input_dir folds --output_dir sorted_folds 
```
The best models were moved to one file and renamed for simplicity. 

```
for f in fold_*_model_*.cif; do
    new="${f#fold_}"
    new="${new/_model_*/}.cif"   
     mv -n -- "$f" "$new"
done
```
Resulting models in folder "sorted_folds."

### Inputted into Fold Mason MSA

required package: foldmason from bioconda 

```
foldmason easy-msa *.cif result tmp --report-mode 1 --gap-open aa:25

#filter-msa 1 --gap-open aa:25 --gap-extend aa:2,nucl:2 --report-paths 0 --report-mode 2
```
 --gap-open aa:25 - gap open penalty is 25, makes gaps fewer and should help with closely related sequences


### iQtree 

Decided to use combined 3di and AA from FoldMason based on the FoldMason paper (https://doi.org/10.1126/science.ads6733). 

Required package: bioconda::iqtree

```
mkdir alignments && cp result_3di.fa alignments/3di.fa && cp result_aa.fa alignments/aa.fa
iqtree3 -p alignments/ -B 1000 -m MFP+MERGE
```
-p : specifies parts
-B 1000 : 1000 ultrafast boostraps
-m MFP+MERGE : model finder MFP+MERGE model finder with extended model selection and reconstructs tree with best partitioning scheme found https://iqtree.github.io/doc/Command-Reference
guide to partitioned analysis with mixed data by iQTREE - https://iqtree.github.io/doc/Advanced-Tutorial

### Moved to R to visualize tree

R script eLRRtree used for visualization. 

## Protein domain annotation

### Signal peptides 

Full length proteins were put into interproscan and the gff file downloaded. 

### LRR motifs

initially manually entered into LRRsearch - https://lrrsearch.com
then used phytolrr - https://github.com/phytolrr/predict-phytolrr

```
predict-phytolrr -f eLRRs.fasta
```
To convert the tsv file to gff
```
awk -F'\t' 'BEGIN{OFS="\t"; print "##gff-version 3"}
NR==1 {next}  # skip header
{
    seq_id = $1
    offset = $2
    sequence = $3
    score = $4
    motif = $5
    start = offset + 1
    end = start + length(sequence) - 1
    print seq_id, "predict-phytolrr", motif, start, end, score, ".", ".", "sequence=" sequence
}
' lrrpred.tsv > lrrpred.gff
```
### Qx motif

Searched full eLRR file for QxxS/T motif 
```
seqkit locate -p "Q..[TS]" -r -I eLRRs.fasta >> qxxst.tsv
```
Change to gff file. 
```
awk -F'\t' 'BEGIN{OFS="\t"; print "##gff-version 3"}
NR==1 {next}  # skip header
{
    seq_id = $1
    patternName = $2
    pattern = $3
    strand = $4
    start = $5
    end = $6
    matched = $7
    print seq_id, "manual", pattern, start, end, ".", strand, ".", "sequence=" matched
}
' qxxst.tsv > qxxst.gff
```
Note: I manually removed any motifs not in the last LRR after visualization. 

### GFF file compilation and R visualization

Went to R and compiled all gff files together.

## Graphed tree and domain together

R script eLRRtree used for visualization. 

## Figure 2C models 

Individual proteins modeled in alphafold3. Python code was used to select best model and FoldMason was used to align models. 
The cif files were downloaded from FoldMason and visualized in pymol using the following code. 

Disordered regions were trimmed to the first LRR in HaSERK1, PlsLRR and NtLRP1.

All models in "models" folder and along with fasta file containing sequences. 

#### HaSERK1 and PlsLRR

HaSERK1 selected as *Plenodomus lindquistii* infects *Helianthus annus*. 

```
set ray_trace_mode, 3
bg_color white
set antialias, 2
set ray_trace_color, black
set ray_shadow, 1
set cartoon_sampling, 10 
as cartoon
set specular, 0
set power, 0.2
set ambient, 0.8
set direct, 0.6
color brown, STRUCTURE_1 #PlsLRR
color gray70, STRUCTURE_0 #HaSERK1
```


#### OsFLS2 and OsLRR1
```
set ray_trace_mode, 3
bg_color white
set antialias, 2
set ray_trace_color, black
set ray_shadow, 1
set cartoon_sampling, 10 
as cartoon
set specular, 0
set power, 0.2
set ambient, 0.8
set direct, 0.6
color forest, STRUCTURE_1 #OsLRR1
color gray70, STRUCTURE_0 #OsFLS2
```

#### NtSERK NtLRP
```
set ray_trace_mode, 3
bg_color white
set antialias, 2
set ray_trace_color, black
set ray_shadow, 1
set cartoon_sampling, 10 
as cartoon
set specular, 0
set power, 0.2
set ambient, 0.8
set direct, 0.6
color forest, STRUCTURE_0 #ntlrp
color gray70, STRUCTURE_1 #ntserk

```
