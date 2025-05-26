test_that("returns a data.frame with the correct columns and types", {
  wig <- system.file("extdata", "DNaseI_BG3_gr_chr4.bedGraph", package = "HMRFHiC")
  chain <- system.file("extdata", "dm3ToDm6_chr4_only.chain", package = "HMRFHiC")
  te <- system.file("extdata", "dm6_TEs_chr4.gtf", package = "HMRFHiC")
  hic <- system.file("extdata", "BG3_WT_merged_hic_matrix_chr4_100Kb.cool", package = "HMRFHiC")

  expect_true(file.exists(wig))
  expect_true(file.exists(chain))
  expect_true(file.exists(te))
  expect_true(file.exists(hic))

  bb <- get_data(
    file_path      = hic,
    chr            = "chr4",
    start          = 1,
    end            = 100000,
    resolution     = 10000,
    genome_package = "BSgenome.Dmelanogaster.UCSC.dm6",
    acc_wig        = wig,
    chain_file     = chain,
    te_granges     = te
  )

  expect_s3_class(bb, "data.frame")
  expect_named(bb, c("start", "end.i.", "start.j.", "end", "chrom", "GC", "ACC", "TES", "interactions"))
  expect_type(bb$start, "double")
  expect_type(bb$GC, "double")
  expect_type(bb$interactions, "double")
})

test_that("builds full 20kb × 20kb grid over 1–40kb", {
  hic <- system.file("extdata", "BG3_WT_merged_hic_matrix_chr4_100Kb.cool", package = "HMRFHiC")
  bb <- get_data(
    file_path      = hic,
    chr            = "chr4",
    start          = 1,
    end            = 100000,
    resolution     = 20000,
    genome_package = "BSgenome.Dmelanogaster.UCSC.dm6"
  )
  expect_equal(nrow(bb), 5 * 5)
  expect_equal(bb$start[1], 1L)
  expect_equal(bb$end.i.[1], 20000L)
})

test_that("GC, ACC, and TES are in valid ranges or NA", {
  wig <- system.file("extdata", "DNaseI_BG3_gr_chr4.bedGraph", package = "HMRFHiC")
  chain <- system.file("extdata", "dm3ToDm6_chr4_only.chain", package = "HMRFHiC")
  te <- system.file("extdata", "dm6_TEs_chr4.gtf", package = "HMRFHiC")
  hic <- system.file("extdata", "BG3_WT_merged_hic_matrix_chr4_100Kb.cool", package = "HMRFHiC")

  bb <- get_data(
    file_path      = hic,
    chr            = "chr4",
    start          = 1,
    end            = 100000,
    resolution     = 10000,
    genome_package = "BSgenome.Dmelanogaster.UCSC.dm6",
    acc_wig        = wig,
    chain_file     = chain,
    te_granges     = te
  )
  expect_true(all(bb$GC >= 0 & bb$GC <= 1, na.rm = TRUE))
  expect_true(all(is.na(bb$ACC) | bb$ACC >= 0))
  expect_true(all(is.na(bb$TES) | bb$TES >= 0))
})
