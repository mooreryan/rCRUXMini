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
      

