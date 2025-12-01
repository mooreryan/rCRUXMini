Make the DB

```
makeblastdb -in ecoli_rnr.fasta -parse_seqids -taxid_map ecoli_rnr.accession2taxid -dbtype nucl -title 'E. coli RNR' -out ecoli_rnr
```

Test it

```
blastn -query primers.fasta -db ecoli_rnr -outfmt 6 -evalue 3e7 -num_alignments 10000000 -perc_identity 50 -qcov_hsp_perc 90 -reward 2 -task 'blastn-short' -word_size 7
```

No hits! This DB is for ensuring that the pipeline handles some DBs with no hits.
