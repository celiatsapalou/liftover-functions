#' Extract Chromosome region
#'
#' @description extract only a specific genomic region from the reference
#'
#'
#' @param ref_fasta [character/link] reference fasta from which we want to extract specific genomic regions
#' @param chr [character] the desired chromosome to be extracted
#' @param start [numeric] the start position of the desired chromosome region in bp
#' @param end [numeric] the end position of the desired chromosome region in bp
#' @param outfasta_ref [character] the output fasta containing the extracted chr region
#' @return nothing. Only output files written.
#'
#' @export


check_shred<- function(ref_fasta,chr,start,end,outfasta_ref){
  locus <- paste0(chr,":",start,"-",end)
  filename<-paste0("bash -c 'samtools faidx -r <(echo ", locus, ") ",  ref_fasta, " > ", outfasta_ref, "'")
  print(filename)
  system(filename)
  }

#filename <- paste0("samtools faidx -r <(echo '",locus,"') ", ref_fasta, " > ", outfasta_ref)
#system("bash -c 'samtools faidx -r <(echo 'chr7:65500000-65531823') /home/celia/Downloads/assemblies/hg38.fa > /home/celia/Downloads/Inversions_2022/newfile3.fasta'")
#samtools faidx -r <(echo "seq1:3-5") seq.fasta

