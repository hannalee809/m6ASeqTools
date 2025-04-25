test_that("Basic 3-group overlap comparison works", {
  test_list <- list(
    A = c("gene1", "gene2", "gene3"),
    B = c("gene2", "gene3", "gene4"),
    C = c("gene3", "gene4", "gene5")
  )

  result <- compare_m6a_distribution(test_list)
  summary <- attr(result, "summary")

  expect_setequal(result$unique_A, "gene1")
  expect_setequal(result$unique_B, character(0))
  expect_setequal(result$unique_C, "gene5")

  expect_setequal(result$A_B, c("gene2", "gene3"))
  expect_setequal(result$B_C, c("gene3", "gene4"))
  expect_setequal(result$A_C, "gene3")

  expect_setequal(result$common_all, "gene3")

  expect_setequal(result$A_B_only, "gene2")
  expect_setequal(result$B_C_only, "gene4")
  expect_setequal(result$A_C_only, character(0))
})

test_that("Disjunct cases work", {
  test_list <- list(
    A = c("tx1", "tx2"),
    B = c("tx3", "tx4"),
    C = c("tx5", "tx6")
  )
  result <- compare_m6a_distribution(test_list)
  summary <- attr(result, "summary")

  expect_setequal(result$common_all, character(0))
  expect_equal(summary$unique_A, 2)
  expect_equal(summary$unique_B, 2)
  expect_equal(summary$unique_C, 2)

  # All pairwise intersections and exclusives should be empty
  expect_true(all(sapply(result[names(result)[grepl("_only$", names(result))]], length) == 0))
})

test_that("Handles fully overlapping sets", {
  identical <- list(
    G1 = c("x", "y", "z"),
    G2 = c("x", "y", "z"),
    G3 = c("x", "y", "z")
  )

  result <- compare_m6a_distribution(identical)

  expect_setequal(result$common_all, c("x", "y", "z"))
  expect_equal(length(result$unique_G1), 0)
  expect_equal(length(result$G1_G2_only), 0)
  expect_equal(length(result$G2_G3_only), 0)
  expect_equal(length(result$G1_G3_only), 0)
})







