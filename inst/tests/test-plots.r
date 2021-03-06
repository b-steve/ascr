context("Testing plots")

## These check for lack of errors, but not of plot itself.

test_that("location plotting", {
    expect_true(locations(example.data$fits$simple.hn, 1))
    expect_true(locations(example.data$fits$bearing.hn, 2))
    expect_true(locations(example.data$fits$bearing.hn, 3, plot.arrows = FALSE))
    expect_true(locations(example.data$fits$simple.hr, 4, xlim = c(0, 5)))
    expect_true(locations(example.data$fits$simple.hr, 5, ylim = c(-5, 5)))
    expect_true(locations(example.data$fits$simple.hn, 6, levels = seq(0.1, 0.9, 0.1)))
    expect_true(locations(example.data$fits$bearing.hn, 7, nlevels = 3))
    expect_true(locations(example.data$fits$simple.hr, 8, density = TRUE))
    expect_true(locations(example.data$fits$simple.hn, 9, show.legend = FALSE))
    expect_true(locations(example.data$fits$bearing.hn, 10, add = TRUE))
})

test_that("detection function plotting", {
    expect_null(show.detfn(example.data$fits$simple.hn))
    expect_null(show.detfn(example.data$fits$simple.hr, add = TRUE))
    expect_null(show.detfn(example.data$fits$bearing.hn, main = "A title"))
})

test_that("survey plotting", {
    expect_null(show.survey(example.data$fits$simple.hn), is_null())
    expect_null(show.survey(example.data$fits$bearing.hn), is_null())
})

test_that("detection surface plotting", {
    expect_null(show.detsurf(example.data$fits$simple.hn), is_null())
    expect_null(show.detsurf(example.data$fits$simple.hr), is_null())
    expect_null(show.detsurf(example.data$fits$bearing.hn), is_null())
})

test_that("inhomogeneous density surface plotting", {
    simple.capt <- example.data$capt[1]
    fit <- fit.ascr(capt = simple.capt, traps = example.data$traps,
                    mask = example.data$mask, fix = list(g0 = 1),
                    ihd.opts = list(model = ~ x + y))
    expect_null(show.Dsurf(fit), is_null())
})

test_that("capture history plotting", {
    expect_null(show.capt(multi.example.data$capt, multi.example.data$traps, multi.example.data$mask, ask = FALSE))
})
