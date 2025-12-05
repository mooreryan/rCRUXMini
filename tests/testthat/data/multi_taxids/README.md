Make the DB

```
makeblastdb -in ecoli_genes.fasta -parse_seqids -taxid_map taxid_map.tsv -dbtype nucl -title 'E. coli Genes' -out ecoli_genes
```
