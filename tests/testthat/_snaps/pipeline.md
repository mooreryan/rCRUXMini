# pipeline works

    Code
      system2("ls", output_directory_path, stdout = TRUE)
    Output
      [1] "amplicon_data.tsv"
      [2] "pipeline_results.rds"
      [3] "plausible_amplicons_coordinates_distinct_taxonomic_ranks.tsv"
      [4] "plausible_amplicons_coordinates_with_taxonomy.tsv"
      [5] "primer_blast.tsv"
      [6] "primer_blast_results.tsv"
      [7] "primers.fasta"

---

    Code
      output_data_files
    Output
      $amplicon_data.tsv
      $amplicon_data.tsv[[1]]
      blast_db_path	index	blast_ordinal_id	accession	sequence_hash_value	seq...

      $amplicon_data.tsv[[2]]
      1	1	11	sequence_11	0X7C670790	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCC...

      $amplicon_data.tsv[[3]]
      1	2	15	sequence_15	0X1027983E	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCC...

      $amplicon_data.tsv[[4]]
      1	3	17	sequence_17	0X58F529BF	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCC...

      $amplicon_data.tsv[[5]]
      1	4	18	sequence_18	0XAB0C9099	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCC...

      $amplicon_data.tsv[[6]]
      1	5	2	sequence_2	0XE97290F0	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCC...

      $amplicon_data.tsv[[7]]
      1	6	21	sequence_21	0XC55CED3B	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCC...

      $amplicon_data.tsv[[8]]
      1	7	22	sequence_22	0X71ABE309	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCC...

      $amplicon_data.tsv[[9]]
      1	8	27	sequence_27	0XDC7C7B5	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCC...

      $amplicon_data.tsv[[10]]
      1	9	28	sequence_28	0X5E814903	AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCC...


      $plausible_amplicons_coordinates_distinct_taxonomic_ranks.tsv
      $plausible_amplicons_coordinates_distinct_taxonomic_ranks.tsv[[1]]
      superkingdom	phylum	class	order	family	genus	species...

      $plausible_amplicons_coordinates_distinct_taxonomic_ranks.tsv[[2]]
      1	3	3	4	8	15	34...


      $plausible_amplicons_coordinates_with_taxonomy.tsv
      $plausible_amplicons_coordinates_with_taxonomy.tsv[[1]]
      accession	gi	staxids	forward_start	forward_stop	forward_mismatch	rever...

      $plausible_amplicons_coordinates_with_taxonomy.tsv[[2]]
      sequence_11	0	1861716	5724	5770	0	6043	5997	1	319	1861716	Diplosphaera...

      $plausible_amplicons_coordinates_with_taxonomy.tsv[[3]]
      sequence_15	0	1400020	6665	6711	1	7304	7258	0	639	1400020	Deuterostich...

      $plausible_amplicons_coordinates_with_taxonomy.tsv[[4]]
      sequence_17	0	1077733	156	202	0	678	632	0	522	1077733	Diplosphaera sp....

      $plausible_amplicons_coordinates_with_taxonomy.tsv[[5]]
      sequence_18	0	1909397	14	60	0	520	474	0	506	1909397	Prasiola cf. delic...

      $plausible_amplicons_coordinates_with_taxonomy.tsv[[6]]
      sequence_2	0	2735166	48	94	0	366	320	0	318	2735166	NA	NA	Viridiplantae...

      $plausible_amplicons_coordinates_with_taxonomy.tsv[[7]]
      sequence_21	0	2980489	264	310	0	701	655	1	437	2980489	Trimyema sp. str...

      $plausible_amplicons_coordinates_with_taxonomy.tsv[[8]]
      sequence_22	0	1461498	3579	3625	0	4149	4103	0	570	1461498	Stichococcus...

      $plausible_amplicons_coordinates_with_taxonomy.tsv[[9]]
      sequence_27	0	40805	7496	7542	0	7980	7934	1	484	40805	NA	NA	NA	Cilioph...

      $plausible_amplicons_coordinates_with_taxonomy.tsv[[10]]
      sequence_28	0	2749887	5319	5365	0	5890	5844	0	571	2749887	Stichococcus...


      $primer_blast.tsv
      $primer_blast.tsv[[1]]
      qseqid	sgi	saccver	mismatch	sstart	send	staxids...

      $primer_blast.tsv[[2]]
      forward_1	0	sequence_99	0	2380	2426	5987...

      $primer_blast.tsv[[3]]
      forward_1	0	sequence_97	0	172	218	2053061...

      $primer_blast.tsv[[4]]
      forward_1	0	sequence_95	0	50	96	1861716...

      $primer_blast.tsv[[5]]
      forward_1	0	sequence_95	1	49	95	1861716...

      $primer_blast.tsv[[6]]
      forward_1	0	sequence_95	1	51	97	1861716...

      $primer_blast.tsv[[7]]
      forward_1	0	sequence_94	0	184	230	1144366...

      $primer_blast.tsv[[8]]
      forward_1	0	sequence_91	0	413	459	3422551...

      $primer_blast.tsv[[9]]
      forward_1	0	sequence_90	0	88	134	1861739...

      $primer_blast.tsv[[10]]
      forward_1	0	sequence_88	0	32	78	2815379...


      $primer_blast_results.tsv
      $primer_blast_results.tsv[[1]]
      qseqid	sgi	saccver	mismatch	sstart	send	staxids...

      $primer_blast_results.tsv[[2]]
      forward_1	0	sequence_0	0	335	381	1861737...

      $primer_blast_results.tsv[[3]]
      forward_1	0	sequence_0	1	336	382	1861737...

      $primer_blast_results.tsv[[4]]
      reverse_1	0	sequence_0	0	1009	963	1861737...

      $primer_blast_results.tsv[[5]]
      reverse_1	0	sequence_0	1	1010	964	1861737...

      $primer_blast_results.tsv[[6]]
      forward_1	0	sequence_11	0	5724	5770	1861716...

      $primer_blast_results.tsv[[7]]
      reverse_1	0	sequence_11	1	6043	5997	1861716...

      $primer_blast_results.tsv[[8]]
      reverse_1	0	sequence_11	0	6044	5998	1861716...

      $primer_blast_results.tsv[[9]]
      forward_1	0	sequence_15	1	6665	6711	1400020...

      $primer_blast_results.tsv[[10]]
      forward_1	0	sequence_15	0	6666	6712	1400020...


      $primers.fasta
      $primers.fasta[[1]]
      >forward_1...

      $primers.fasta[[2]]
      AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCC...

      $primers.fasta[[3]]
      >reverse_1...

      $primers.fasta[[4]]
      GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTT...

---

    Code
      result
    Output
      $primer_blast_results
      # A tibble: 131 x 7
         qseqid    sgi   saccver     mismatch sstart  send staxids
         <chr>     <chr> <chr>          <int>  <int> <int> <chr>
       1 forward_1 0     sequence_0         0    335   381 1861737
       2 forward_1 0     sequence_0         1    336   382 1861737
       3 reverse_1 0     sequence_0         0   1009   963 1861737
       4 reverse_1 0     sequence_0         1   1010   964 1861737
       5 forward_1 0     sequence_11        0   5724  5770 1861716
       6 reverse_1 0     sequence_11        1   6043  5997 1861716
       7 reverse_1 0     sequence_11        0   6044  5998 1861716
       8 forward_1 0     sequence_15        1   6665  6711 1400020
       9 forward_1 0     sequence_15        0   6666  6712 1400020
      10 reverse_1 0     sequence_15        0   7304  7258 1400020
      # i 121 more rows

      $plausible_amplicons_coordinates_with_taxonomy
      # A tibble: 38 x 33
         accession   gi    staxids forward_start forward_stop forward_mismatch
         <chr>       <chr> <chr>           <int>        <int>            <int>
       1 sequence_11 0     1861716          5724         5770                0
       2 sequence_15 0     1400020          6665         6711                1
       3 sequence_17 0     1077733           156          202                0
       4 sequence_18 0     1909397            14           60                0
       5 sequence_2  0     2735166            48           94                0
       6 sequence_21 0     2980489           264          310                0
       7 sequence_22 0     1461498          3579         3625                0
       8 sequence_27 0     40805            7496         7542                0
       9 sequence_28 0     2749887          5319         5365                0
      10 sequence_31 0     732243             33           79                0
      # i 28 more rows
      # i 27 more variables: reverse_start <int>, reverse_stop <int>,
      #   reverse_mismatch <int>, product_length <int>, taxonomy_id <int>,
      #   species <chr>, superkingdom <chr>, kingdom <chr>, phylum <chr>,
      #   subphylum <chr>, superclass <chr>, class <chr>, subclass <chr>,
      #   order <chr>, family <chr>, subfamily <chr>, genus <chr>, infraorder <chr>,
      #   subcohort <chr>, superorder <chr>, superfamily <chr>, tribe <chr>, ...

      $plausible_amplicons_coordinates_distinct_taxonomic_ranks
      # A tibble: 1 x 7
        superkingdom phylum class order family genus species
               <int>  <int> <int> <int>  <int> <int>   <int>
      1            1      3     3     4      8    15      34

      $amplicon_blast_result
      # A tibble: 61 x 10
         query_accession subject_accession_version percent_identical_matches
         <chr>           <chr>                                         <dbl>
       1 sequence_11     sequence_11                                   100
       2 sequence_15     sequence_15                                   100
       3 sequence_15     sequence_16                                    98.7
       4 sequence_17     sequence_17                                   100
       5 sequence_18     sequence_18                                   100
       6 sequence_18     sequence_19                                    97.8
       7 sequence_2      sequence_2                                    100
       8 sequence_21     sequence_21                                   100
       9 sequence_22     sequence_22                                   100
      10 sequence_22     sequence_23                                    96.7
      # i 51 more rows
      # i 7 more variables: alignment_length <int>, expect_value <dbl>,
      #   subject_sequence_length <int>, subject_alignment_start <int>,
      #   subject_alignment_end <int>, subject_aligned_sequence <chr>,
      #   unique_subject_taxonomy_ids <chr>

      $parsed_amplicon_blast_result
      # A tibble: 61 x 39
         query_accession subject_accession_version percent_identical_matches
         <chr>           <chr>                                         <dbl>
       1 sequence_9      sequence_10                                    98.1
       2 sequence_99     sequence_100                                   98.3
       3 sequence_11     sequence_11                                   100
       4 sequence_15     sequence_15                                   100
       5 sequence_15     sequence_16                                    98.7
       6 sequence_17     sequence_17                                   100
       7 sequence_18     sequence_18                                   100
       8 sequence_18     sequence_19                                    97.8
       9 sequence_2      sequence_2                                    100
      10 sequence_21     sequence_21                                   100
      # i 51 more rows
      # i 36 more variables: alignment_length <int>, expect_value <dbl>,
      #   subject_sequence_length <int>, subject_alignment_start <int>,
      #   subject_alignment_end <int>, subject_aligned_sequence <chr>,
      #   unique_subject_taxonomy_ids <chr>, degapped_subject_aligned_sequence <chr>,
      #   degapped_alignment_length <int>, taxonomy_id <int>, species <chr>,
      #   superkingdom <chr>, kingdom <chr>, phylum <chr>, subphylum <chr>, ...
