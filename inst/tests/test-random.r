context("Random")

## Some tests:
test_that("Test RNG", {
    set.seed(1234)
    rnos <- sample(1:100, size = 5)
    expect_that(rnos, equals(c(12, 62, 60, 61, 83)))
})
