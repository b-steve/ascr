context("Testing plots")

## These check for lack of errors, but not of plot itself.

test_that("location plotting", {
    expect_that(locations(example$fits$simple.hn, 1), is_true())
    expect_that(locations(example$fits$bearing.hn, 2), is_true())
    expect_that(locations(example$fits$bearing.hn, 3, plot.arrows = FALSE), is_true())
    expect_that(locations(example$fits$simple.hr, 4, xlim = c(0, 5)), is_true())
    expect_that(locations(example$fits$simple.hr, 5, ylim = c(-5, 5)), is_true())
    expect_that(locations(example$fits$simple.hn, 6, levels = seq(0.1, 0.9, 0.1)), is_true())
    expect_that(locations(example$fits$bearing.hn, 7, nlevels = 3), is_true())
    expect_that(locations(example$fits$simple.hr, 8, density = TRUE), is_true())
    expect_that(locations(example$fits$simple.hn, 9, show.legend = FALSE), is_true())
    expect_that(locations(example$fits$bearing.hn, 10, add = TRUE), is_true())
})

test_that("detection function plotting", {
    expect_null(show.detfn(example$fits$simple.hn))
    expect_null(show.detfn(example$fits$simple.hr, add = TRUE))
    expect_null(show.detfn(example$fits$bearing.hn, main = "A title"))
})

test_that("survey plotting", {
    expect_null(show.survey(example$fits$simple.hn), is_null())
    expect_null(show.survey(example$fits$bearing.hn), is_null())
})

test_that("detection surface plotting", {
    expect_null(show.detsurf(example$fits$simple.hn), is_null())
    expect_null(show.detsurf(example$fits$simple.hr), is_null())
    expect_null(show.detsurf(example$fits$bearing.hn), is_null())
})

test_that("inhomogeneous density surface plotting", {
    simple.capt <- example$capt[1]
    fit <- fit.ascr(capt = simple.capt, traps = example$traps,
                    mask = example$mask, fix = list(g0 = 1),
                    ihd.opts = list(model = ~ x + y))
    expect_null(show.Dsurf(fit), is_null())
})
