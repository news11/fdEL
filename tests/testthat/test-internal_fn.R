test_that("mu_a_hat works", {
  expect_equal(mu_a_hat(Xt <- c(2,4,2,0), rep(0.25, 4), c(0, 2, 4)), c(0.75, 0.25, 0))
})
