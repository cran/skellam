test_that("rskellam returns integers, like rpois  [N17]", {
  set.seed(1)
  z <- rskellam(10, 3, 2)
  expect_type(z, "integer")
  set.seed(1)
  expect_identical(z, rpois(10, 3) - rpois(10, 2))   # same draws as before
  expect_identical(rskellam(0, 3, 2), integer(0))
})

test_that("invalid rates give NA with a warning, like rpois", {
  expect_warning(z <- rskellam(3, c(1, -1, 2), 2), "NAs produced")
  expect_true(is.na(z[2]))
  expect_false(anyNA(z[-2]))
})
