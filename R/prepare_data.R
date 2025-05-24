#' @title Extract Hi-C Bin Interactions with GC, Accessibility, and TE Counts
#'
#' @description
#' Reads a Hi-C contact matrix from a \code{.hic}, \code{.cool}, or \code{.h5} file,
#' subsets to a specified genomic region, rebins to a desired resolution (merging or
#' slicing as needed), and returns a complete grid of bin-pair interactions. For each
#' bin-pair it computes:
#' \itemize{
#'   \item GC content from a BSgenome
#'   \item DNase-I accessibility (if provided, lifted over)
#'   \item Transposable element (TE) overlap counts
#'   \item Interaction counts
#' }
#'
#' @usage
#' get_data(
#'   file_path, chr, start, end, resolution,
#'   genome_package, acc_wig = NULL,
#'   chain_file = NULL, te_granges = NULL
#' )
#'
#' @param file_path     Path to the Hi-C file (\code{.hic}, \code{.cool}, or \code{.h5}).
#' @param chr           Chromosome name (e.g. \code{"chr2L"}).
#' @param start         1-based start coordinate of the region.
#' @param end           End coordinate of the region.
#' @param resolution    Target bin size (bp); smaller native bins are merged, larger
#'                      ones are sliced.
#' @param genome_package
#'   Name of a BSgenome package (e.g.
#'   \code{"BSgenome.Dmelanogaster.UCSC.dm6"}) for GC calculations.
#' @param acc_wig       (Optional) Path to a DNase-I wig (bedGraph) file. Requires
#'                      \code{chain_file} to lift over.
#' @param chain_file    (Optional) LiftOver chain file for mapping \code{acc_wig}.
#' @param te_granges    (Optional) Path to a BED/GTF of TE annotations or a
#'                      \code{GRanges} object. Only TEs on \code{chr} are counted.
#'
#' @return
#' A \code{data.frame} with one row per bin-pair across \code{[start,end]} at
#' \code{resolution}, containing:
#' \describe{
#'   \item{\code{start}}{Bin1 start.}
#'   \item{\code{end.i.}}{Bin1 end.}
#'   \item{\code{start.j.}}{Bin2 start.}
#'   \item{\code{end}}{Bin2 end.}
#'   \item{\code{chrom}}{Chromosome name.}
#'   \item{\code{GC}}{Combined GC fraction of the two bins.}
#'   \item{\code{ACC}}{Mean DNase-I score per bin or \code{NA}.}
#'   \item{\code{TES}}{Sum of TE overlaps in both bins or \code{NA}.}
#'   \item{\code{interactions}}{Observed contact count (\code{0} if absent).}
#' }
#'
#' @examples
#' \dontrun{
#' bb <- get_data(
#'   file_path      = "C:/Users/name/Desktop/hicTransform.h5",
#'   chr            = "chr2L",
#'   start          = 1,
#'   end            = 40000,
#'   resolution     = 2000,
#'   genome_package = "BSgenome.Dmelanogaster.UCSC.dm6",
#'   acc_wig        = "other/BG3.r2c.dhs.density.wig",
#'   chain_file     = "dm3ToDm6.over.chain",
#'   te_granges     = "dm6_TEs.gtf"
#' )
#' }
#'
#' @seealso
#' \code{\link[strawr]{straw}}, \code{\link[rhdf5]{h5read}}, \code{\link[rtracklayer]{liftOver}}
#'
#' @importFrom Biostrings getSeq DNAString letterFrequency
#' @importFrom GenomicRanges GRanges seqnames countOverlaps findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer import.chain import liftOver
#' @importFrom strawr straw
#' @importFrom utils read.table
#' @importFrom stats aggregate
#' @importFrom rhdf5 h5ls h5read
#' @export
get_data <- function(file_path, chr, start, end, resolution,
                     genome_package = NULL,
                     acc_wig = NULL, chain_file = NULL,
                     te_granges = NULL) {
  ## 0) Core package check
  required_pkgs <- c("rhdf5", "strawr", "rtracklayer", "GenomicRanges")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs,
    requireNamespace,
    quietly = TRUE, FUN.VALUE = logical(1)
  )]
  if (length(missing_pkgs)) {
    stop(
      "Missing required packages:\n  - ",
      paste(missing_pkgs, collapse = "\n  - "),
      "\nPlease install them before running get_data()."
    )
  }

  ## 0b) Optional BSgenome + Biostrings load for real GC
  has_genome <- FALSE
  genome <- NULL
  if (!is.null(genome_package)) {
    if (requireNamespace("BSgenome", quietly = TRUE) &&
      requireNamespace("Biostrings", quietly = TRUE) &&
      requireNamespace(genome_package, quietly = TRUE)) {
      genome <- get(genome_package)
      has_genome <- TRUE
    } else {
      warning(
        "Requested real GC, but BSgenome/Biostrings or ",
        genome_package,
        " not installed. Falling back to GC = NA."
      )
    }
  }

  ## 1) Import DNase-I (ACC) and liftOver if provided
  acc_gr <- NULL
  if (!is.null(acc_wig) && !is.null(chain_file)) {
    dnase_df <- utils::read.table(acc_wig,
      header = FALSE, skip = 1,
      col.names = c("chr", "start", "end", "score"),
      stringsAsFactors = FALSE
    )
    if (!grepl("^chr", dnase_df$chr[1])) {
      dnase_df$chr <- paste0("chr", dnase_df$chr)
    }
    acc_gr <- GenomicRanges::GRanges(dnase_df$chr,
      IRanges::IRanges(dnase_df$start + 1, dnase_df$end),
      score = as.numeric(dnase_df$score)
    )
    chain <- rtracklayer::import.chain(chain_file)
    acc_gr <- unlist(rtracklayer::liftOver(acc_gr, chain))
    acc_gr <- acc_gr[GenomicRanges::seqnames(acc_gr) == chr]
  }

  ## 2) Import TE ranges if provided
  te_in <- NULL
  if (!is.null(te_granges)) {
    if (is.character(te_granges)) {
      te_in <- rtracklayer::import(te_granges)
    } else if (inherits(te_granges, "GRanges")) {
      te_in <- te_granges
    } else {
      stop("`te_granges` must be a file path or a GRanges object")
    }
    te_in <- te_in[GenomicRanges::seqnames(te_in) == chr]
  }

  ## 3) Detect .h5 schema and read raw interactions
  hc <- rhdf5::h5ls(file_path)
  has_bins <- any(hc$group == "/" & hc$name == "bins")
  has_intervals <- any(hc$group == "/" & hc$name == "intervals")

  read_cool <- function(fp, ch) {
    codes_raw <- rhdf5::h5read(fp, "bins/chrom")
    names_map <- rhdf5::h5read(fp, "chroms/name")
    if (is.factor(codes_raw)) {
      cl <- as.character(codes_raw)
      bin_chrom <- if (all(cl %in% names_map)) paste0("chr", cl) else cl
    } else {
      bin_chrom <- paste0("chr", names_map[codes_raw + 1])
    }
    bin_start <- rhdf5::h5read(fp, "bins/start") + 1
    bin_end <- rhdf5::h5read(fp, "bins/end")
    bins_tbl <- data.frame(
      bin_id = seq_along(bin_chrom) - 1,
      chrom = bin_chrom,
      start = bin_start,
      end = bin_end,
      stringsAsFactors = FALSE
    )
    px <- data.frame(
      bin1_id = rhdf5::h5read(fp, "pixels/bin1_id"),
      bin2_id = rhdf5::h5read(fp, "pixels/bin2_id"),
      count = rhdf5::h5read(fp, "pixels/count"),
      stringsAsFactors = FALSE
    )
    df1 <- merge(px, bins_tbl, by.x = "bin1_id", by.y = "bin_id", sort = FALSE)
    names(df1)[4:6] <- c("bin1_chrom", "bin1_start", "bin1_end")
    df2 <- merge(df1, bins_tbl, by.x = "bin2_id", by.y = "bin_id", sort = FALSE)
    names(df2)[7:9] <- c("bin2_chrom", "bin2_start", "bin2_end")
    subset(df2, bin1_chrom == ch & bin2_chrom == ch)
  }

  read_interval <- function(fp, ch) {
    ints <- rhdf5::h5read(fp, "intervals")
    ints$start_list <- ints$start_list + 1
    bins_tbl <- data.frame(
      bin_id = seq_along(ints$start_list) - 1,
      chrom = ints$chr_list,
      start = ints$start_list,
      end = ints$end_list,
      stringsAsFactors = FALSE
    )
    dat <- rhdf5::h5read(fp, "matrix")
    rows <- rep(seq_len(dat$shape[1]) - 1, diff(dat$indptr))
    px <- data.frame(
      bin1_id = rows,
      bin2_id = dat$indices,
      count = dat$data,
      stringsAsFactors = FALSE
    )
    df1 <- merge(px, bins_tbl, by.x = "bin1_id", by.y = "bin_id", sort = FALSE)
    names(df1)[4:6] <- c("bin1_chrom", "bin1_start", "bin1_end")
    df2 <- merge(df1, bins_tbl, by.x = "bin2_id", by.y = "bin_id", sort = FALSE)
    names(df2)[7:9] <- c("bin2_chrom", "bin2_start", "bin2_end")
    subset(df2, bin1_chrom == ch & bin2_chrom == ch)
  }

  raw_df <- if (has_bins) {
    read_cool(file_path, chr)
  } else if (has_intervals) {
    read_interval(file_path, chr)
  } else {
    stop("Unrecognized .h5 schema: neither /bins nor /intervals found")
  }

  ## 4) Subset to the requested region
  raw_df <- subset(
    raw_df,
    bin1_start <= end & bin1_end >= start &
      bin2_start <= end & bin2_end >= start
  )

  ## 5) Determine native resolution & re-bin if needed
  if (nrow(raw_df) > 0) {
    native_res <- unique(raw_df$bin1_end - raw_df$bin1_start + 1)
    native_res <- if (length(native_res) > 1) min(native_res) else native_res
    obs_df <- raw_df
  } else {
    native_res <- resolution
    obs_df <- raw_df[FALSE, ]
  }

  if (native_res != resolution && nrow(obs_df) > 0) {
    if (native_res < resolution) {
      obs_df$new_s1 <- floor((obs_df$bin1_start - 1) / resolution) * resolution + 1
      obs_df$new_e1 <- obs_df$new_s1 + resolution - 1
      obs_df$new_s2 <- floor((obs_df$bin2_start - 1) / resolution) * resolution + 1
      obs_df$new_e2 <- obs_df$new_s2 + resolution - 1
      agg <- stats::aggregate(count ~ new_s1 + new_e1 + new_s2 + new_e2,
        data = obs_df, FUN = sum
      )
      obs_df <- data.frame(
        bin1_start = agg$new_s1, bin1_end = agg$new_e1,
        bin2_start = agg$new_s2, bin2_end = agg$new_e2,
        count = agg$count,
        stringsAsFactors = FALSE
      )
    } else {
      obs_df <- do.call(rbind, apply(obs_df, 1, function(r) {
        s1 <- as.integer(r["bin1_start"])
        e1 <- as.integer(r["bin1_end"])
        s2 <- as.integer(r["bin2_start"])
        e2 <- as.integer(r["bin2_end"])
        cnt <- as.numeric(r["count"])
        s1v <- seq(s1, e1, by = resolution)
        e1v <- pmin(s1v + resolution - 1, e1)
        s2v <- seq(s2, e2, by = resolution)
        e2v <- pmin(s2v + resolution - 1, e2)
        comb <- expand.grid(i = seq_along(s1v), j = seq_along(s2v))
        out <- data.frame(
          bin1_start = s1v[comb$i], bin1_end = e1v[comb$i],
          bin2_start = s2v[comb$j], bin2_end = e2v[comb$j],
          stringsAsFactors = FALSE
        )
        w1 <- (out$bin1_end - out$bin1_start + 1) / (e1 - s1 + 1)
        w2 <- (out$bin2_end - out$bin2_start + 1) / (e2 - s2 + 1)
        out$count <- cnt * w1 * w2
        out
      }))
    }
  }

  ## 6) Build the full grid and seed zero counts
  bs <- seq(start, end, by = resolution)
  be <- pmin(bs + resolution - 1, end)
  full <- expand.grid(
    bin1_start = bs, bin2_start = bs,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  full$bin1_end <- be[match(full$bin1_start, bs)]
  full$bin2_end <- be[match(full$bin2_start, bs)]
  full$interactions <- 0
  if (nrow(obs_df) > 0) {
    keys_full <- paste(full$bin1_start, full$bin1_end,
      full$bin2_start, full$bin2_end,
      sep = "_"
    )
    keys_obs <- paste(obs_df$bin1_start, obs_df$bin1_end,
      obs_df$bin2_start, obs_df$bin2_end,
      sep = "_"
    )
    full$interactions[match(keys_obs, keys_full)] <- obs_df$count
  }

  ## 7) Compute GC content
  if (has_genome) {
    full$GC <- mapply(
      function(s1, e1, s2, e2) {
        seq1 <- Biostrings::getSeq(genome, chr, s1, e1)
        seq2 <- Biostrings::getSeq(genome, chr, s2, e2)
        dna <- Biostrings::DNAString(paste0(as.character(seq1), as.character(seq2)))
        mean(Biostrings::letterFrequency(dna, c("G", "C"), as.prob = TRUE))
      },
      full$bin1_start, full$bin1_end,
      full$bin2_start, full$bin2_end
    )
  } else {
    full$GC <- NA_real_
  }

  ## 8) Compute ACC
  full$ACC <- if (is.null(acc_gr)) {
    NA_real_
  } else {
    mapply(
      function(s1, e1, s2, e2) {
        g1 <- GenomicRanges::GRanges(chr, IRanges::IRanges(s1, e1))
        g2 <- GenomicRanges::GRanges(chr, IRanges::IRanges(s2, e2))
        hits1 <- GenomicRanges::findOverlaps(acc_gr, g1)
        hits2 <- GenomicRanges::findOverlaps(acc_gr, g2)
        m1 <- if (length(hits1) > 0) mean(acc_gr$score[S4Vectors::queryHits(hits1)], na.rm = TRUE) else NA_real_
        m2 <- if (length(hits2) > 0) mean(acc_gr$score[S4Vectors::queryHits(hits2)], na.rm = TRUE) else NA_real_
        mean(c(m1, m2), na.rm = TRUE)
      },
      full$bin1_start, full$bin1_end,
      full$bin2_start, full$bin2_end
    )
  }

  ## 9) Compute TES
  full$TES <- if (is.null(te_in)) {
    NA_integer_
  } else {
    mapply(
      function(s1, e1, s2, e2) {
        g1 <- GenomicRanges::GRanges(chr, IRanges::IRanges(s1, e1))
        g2 <- GenomicRanges::GRanges(chr, IRanges::IRanges(s2, e2))
        sum(GenomicRanges::countOverlaps(g1, te_in)) + sum(GenomicRanges::countOverlaps(g2, te_in))
      },
      full$bin1_start, full$bin1_end,
      full$bin2_start, full$bin2_end
    )
  }

  ## 10) Return final data.frame
  data.frame(
    start = full$bin1_start,
    `end.i.` = full$bin1_end,
    `start.j.` = full$bin2_start,
    end = full$bin2_end,
    chrom = chr,
    GC = full$GC,
    ACC = full$ACC,
    TES = full$TES,
    interactions = round(as.double(full$interactions)),
    stringsAsFactors = FALSE
  )
}
