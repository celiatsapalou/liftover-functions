#' liftover_coarse
#' Wrapper function to find a liftover sequence from one assembly to the other.
#' Better description TBD
#' @export


#first I need to choose the assembly of interest
#then I need to choose the region of interest in the reference

chr= "chr16"
start = 14000000
end = 15500000


#infasta=assembly
#chunkleng=length of chop-up
#outfasta_chunk=Destination of output shreds.
# shred_function
shred_seq <- function(infasta,
                      outfasta_chunk,
                      chunklen) {
  scriptloc = query_config("shred")
  print(paste0(
    scriptloc,
    " in=",
    infasta,
    " out=",
    outfasta_chunk,
    " length=",
    chunklen
  ))
  system(
    paste0(
      scriptloc,
      " in=",
      infasta,
      " out=",
      outfasta_chunk,
      " length=",
      chunklen,
      " overwrite=true"
    )
  )
}





# infasta= an assembly
# outfasta chunk= directory for the output file
# chunk length: number of chopped up reads
source('shred_seq.R')
source('run_minimap2.R')


# new command
liftover_coarse<- function(chr,start,end,infasta_assembly,
                           chunklen, ref_fasta, outfasta_ref) {
  output1 <- shred_seq(infasta_assembly, chunklen, outfasta_chunk)
  #read specific region from fasta file
  output2<-system(paste0(" samtools faidx -r ", ref_fasta, "< (echo ", chr:start-end,  " )", outfasta_ref))
  output3<- run_minimap2 (outfasta_ref, output1 ) #THE OUTPUT NEEDS TO BE the chr, start, end inside a paf file

}

#samtools faidx -r <(echo "seq1:3-5") seq.fasta
#system(paste0(" samtools faidx -r ", ref_fasta, "< (echo ", chr:start-end,  " )", outfasta_ref))#

output2<-system(paste0(" samtools faidx -r < (echo ", chr:start-end,  " )", ref_fasta, ">", outfasta_ref))

#apply it
shred_seq(queryfasta, queryfasta_chunk, chunklen)

#targetfasta= part of the reference genome
#queryfasta= the chopped up assembly
#output= paf file


#' Submit a system command to run minimap2
#'
#' @description This is a helperfunction to run minimap2. Also check out:
#' https://github.com/PacificBiosciences/pbmm2/ . Minimap2 parameters:
#'  -k   k-mer size (no larger than 28). [-1]
#' -w   Minimizer window size. [-1]
#' -u   Disable homopolymer-compressed k-mer (compression is active for SUBREAD & UNROLLED presets).
#' -A   Matching score. [-1]
#' -B   Mismatch penalty. [-1]
#' -z   Z-drop score. [-1]
#' -Z   Z-drop inversion score. [-1]
#' -r   Bandwidth used in chaining and DP-based alignment. [-1]
#' -g   Stop chain enlongation if there are no minimizers in N bp. [-1]
#' #'
#' @param targetfasta [character/link] link to the 'target' single-sequence fasta (sometimes reference, e.g. chm13.)
#' @param queryfasta [character/link] link to the 'query' fasta. Can be single or multi-fasta
#' @param outpaf [character/link] Path to the output paffile to be written.
#' @param minimap2loc [character/link] link to minimap2 binary.

#' @return nothing. Only output files written.
#'
#' @author Wolfram Höps
#' @export
run_minimap2 <-
  function(targetfasta,
           queryfasta,
           outpaf,
           nthreads = 4) {
    #system(paste0(minimap2loc," -x asm20 -c -z400,50 -s 0 -M 0.2 -N 100 -P --hard-mask-level ", fastatarget, " ", fastaquery, " > ", outpaf))

    minimap2loc = query_config("minimap2")

    # Some self-defined parameters
    system(
      paste0(
        minimap2loc,
        " -x asm20 -P -c -s 0 -M 0.2 -t ",
        nthreads,
        " ",
        targetfasta, #referece
        " ",
        queryfasta, #assembly
        " > ",
        outpaf
      )
    )
    # Check if that was successful.
    stopifnot("Alignment error: Minimap2 has not reported any significant alignment.
              Check if your input sequence is sufficiently long." =
                file.size(outpaf) != 0)}
























# apply liftover_coarse
# Get coords in assembly
coords_liftover = liftover_coarse(seqname, start, end, conversionpaf_link, lenfactor = 1)


liftover_coarse <-
  function(seqname,
           start,
           end,
           lenfactor = 1.2,
           whole_chr = F) {



    # What is the majority vote for the target chromosome?
    winner_chr = names(sort(table(liftover_coords$seqname), decreasing = TRUE)[1])


    #
    # Make sure we don't exceed chromosome boundaries in the query.
    start_winners_cutoff = as.integer(max(0, start_winners))
    end_winners_cutoff =   as.integer(min(cpaf[cpaf$tname == winner_chr,][1, 'tlen'], end_winners))



    return(
      list(
        lift_contig = winner_chr,
        lift_start = start_winners_cutoff,
        lift_end = end_winners_cutoff
      )
    )

  }
