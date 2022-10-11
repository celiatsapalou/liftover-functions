#' Chunkify query fasta
#'
#' @description This is a helper function calling an external script to
#' chop a query sequence into chunks.
#'
#' @param infasta [character/link] single-seq fasta to be chopped
#' @param outfasta_chunk [character/link] output chopped multi-seq fasta.
#' @param chunklen [numeric] length of sequence chunks in bp
#' @param scriptloc [character/link] link to shred.ss from bbmap.
#' @return nothing. Only output files written.
#'
#' @export


shred_seq <- function(infasta,
                      chunklen,
                      outfasta_chunk) {
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
