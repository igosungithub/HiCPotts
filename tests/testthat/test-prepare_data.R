skip_if_no_extdata <- function() {
  hic <- system.file("extdata", "BG3_WT_merged_hic_matrix_chr4_100Kb.cool",
                     package = "HiCPotts")
  if (!nzchar(hic) || !file.exists(hic))
    testthat::skip("extdata .cool file not available in installed package")
  hic
}

test_that(".hic files are explicitly rejected (regression test)", {
  tmp <- tempfile(fileext = ".hic")
  file.create(tmp)
  on.exit(unlink(tmp))
  expect_error(
    get_data(tmp, chr = "chr1", start = 1, end = 1000, resolution = 100),
    "\\.hic"
  )
})

test_that("Returns a data.frame with the expected columns", {
  hic <- skip_if_no_extdata()
  bb <- get_data(
    file_path  = hic,
    chr            = "chr4",
    start          = 1,
    end            = 400000,
    resolution     = 200000,
  )
  expect_s3_class(bb, "data.frame")
  expect_true(all(c("start","end.i.","start.j.","end","chrom",
                    "GC","ACC","TES","interactions") %in% names(bb)))
  expect_equal(nrow(bb), 2 * 2)
  expect_equal(bb$start[1], 1L)
  expect_equal(bb$end.i.[1], 200000L)
})

test_that("Interactions are non-negative real values (no forced rounding)", {
  hic <- skip_if_no_extdata()
  bb <- get_data(
    file_path  = hic,
    chr            = "chr4",
    start          = 1,
    end            = 400000,
    resolution     = 200000,
  )
  expect_type(bb$interactions, "double")       # FIX: double, not round-to-int
  expect_true(all(is.finite(bb$interactions)))
  expect_true(all(bb$interactions >= 0))
})
