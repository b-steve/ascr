context("Testing plots")

## These check for lack of errors, but not of plot itself.

test_that("location plotting", {
    expect_that(locations(simple.hn.fit, 1), is_true())
    expect_that(locations(bearing.hn.fit, 2), is_true())
    expect_that(locations(bearing.hn.fit, 3, plot.arrows = FALSE), is_true())
    expect_that(locations(simple.hr.fit, 4, xlim = c(0, 5)), is_true())
    expect_that(locations(simple.hr.fit, 5, ylim = c(-5, 5)), is_true())
    expect_that(locations(simple.hn.fit, 6, levels = seq(0.1, 0.9, 0.1)), is_true())
    expect_that(locations(bearing.hn.fit, 7, nlevels = 3), is_true())
    expect_that(locations(simple.hr.fit, 8, density = TRUE), is_true())
    expect_that(locations(simple.hn.fit, 9, legend = FALSE), is_true())
    expect_that(locations(bearing.hn.fit, 10, add = TRUE), is_true())
})

test_that("detection function plotting", {
    ## To go here.
})
