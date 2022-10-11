#' @param targetfasta reference genome. the region we want to align to [character]
#' @param queryfasta the data that we want to align to something [character]
#' @param outpaf  [character/link] Path to the output paffile to be written.

run_minimap2 <-
  function(targetfasta,
           queryfasta,
           outpaf,
           nthreads = 4) {
    #system(paste0(minimap2loc," -x asm20 -c -z400,50 -s 0 -M 0.2 -N 100 -P --hard-mask-level ", fastatarget, " ", fastaquery, " > ", outpaf))
    
  
    minimap2loc = query_config("minimap2")
    
    
    # If we have a pre-computed coarse alignment, then we can use this to find out 
    # which region we are talking about. 
    if (use_paf_library) {
      # Pad-sequence
      start_end_pad = enlarge_interval_by_factor(start_x,
                                                 end_x,
                                                 params$xpad,
                                                 seqname_f = seqname_x,
                                                 conversionpaf_f = conversionpaf_link)
      start_x_pad = start_end_pad[1]
      end_x_pad = start_end_pad[2]
     
      
    #get coordinates in y
      #it will translate chromosome start-end to contig start-end (if there is translocation, maybe to more than one contigs? it cannot be multiple)
    coords_liftover = liftover_coarse(seqname_x, #sequence in the reference
                                      start_x_pad,
                                      end_x_pad,
                                      conversionpaf_link,
                                      lenfactor = aln_pad_factor,
                                      whole_chr = params$whole_chr)
    ??liftover_coarse
    
    # Get subseq-fastas in x and y
    extract_subseq_bedtools(genome_x_fa,
                            seqname_x,
                            start_x_pad,
                            end_x_pad,
                            outlinks$genome_x_fa_subseq)
    extract_subseq_bedtools(
      genome_y_fa,
      coords_liftover$lift_contig,
      coords_liftover$lift_start,
      coords_liftover$lift_end,
      outlinks$genome_y_fa_subseq
    )
  } else {
    system(paste0('cp ', genome_x_fa, ' ', outlinks$genome_x_fa_subseq))
    system(paste0('cp ', genome_y_fa, ' ', outlinks$genome_y_fa_subseq))
    
    start_x_pad = 0
    end_x_pad = 1
    start_x = 0
    end_x = 1
    system(paste0('cp ', genome_x_fa, ' ', outlinks$genome_x_fa_subseq))
    system(paste0('cp ', genome_y_fa, ' ', outlinks$genome_y_fa_subseq))
    
  }
    
    # Some self-defined parameters
    system(
      paste0(
        minimap2loc,
        " -x asm20 -P -c -s 0 -M 0.2 -t ", ##see what it's about, p.x. 0.2 can be a cutoff and add additional parameters (px only allow uniquely mapped reads)
        nthreads,
        " ",
        targetfasta,
        " ",
        queryfasta,
        " > ",
        outpaf
      )
    )
    
    # Check if that was successful. 
    stopifnot("Alignment error: Minimap2 has not reported any significant alignment. 
              Check if your input sequence is sufficiently long." = 
                file.size(outpaf) != 0)
    
    
   
  }
    
    
    
    
    
    
    
    
    
    
    #pbmm2: CCS/HIFI
    #system(paste0(minimap2loc," -k 19 -w 10 -u -o 5 -O 56 -e 4 -E 1 -A 2 -B 5 -z 400 -Z 50 -r 2000 -L 0.5 -g 5000", targetfasta, " ", queryfasta, " > ", outpaf))
    
    # W, 23rd Dec 2021. Since the dev of pbmm2, the minimap2 parameters have changed. I adapted everything to the new notation.
    # the -L paramter (formerly Long join flank ratio) is no longer existing, so I'm leaving it out.
    #system(paste0(minimap2loc," -k 19 -w 10 -O 5,56 -E 4,1 -A 2 -B 5 -z 400,50 -r 2000 -g 5000 ", targetfasta, " ", queryfasta, " > ", outpaf))
    
    
  