
test_that("test with dube data", {
  id_col.target <- 2
  t0 <- 2003

  # load datadf <- copy(dube)
  df <- copy(dube)
  disco <- DiSCo(df, id_col.target, t0, G = 1000, num.cores = 1,
                 permutation = TRUE, CI = TRUE, boots = 1000, graph = TRUE, simplex=TRUE, seed=1, q_max=0.9)

  discot <- DiSCoTEA(disco,  agg="quantileDiff", graph=TRUE)
  summary(discot)
  })

test_that("quantile works", {

  # disco
  Ts <- 2
  t0 <- 2
  df <- ex_gmm(Ts=Ts,  num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, CI=TRUE, boots=2)

  expect_no_error(discot <- DiSCoTEA(disco, agg="quantile"))
  expect_true(typeof(discot$plot) == "list")
  expect_no_error(summary(discot))
  expect_true(!is.null(discot$treats))
  expect_true(!is.null(discot$ses))
  expect_true(!is.null(discot$ci_lower))
  expect_true(!is.null(discot$ci_upper))


  # mixture
  Ts <- 2
  t0 <- 2
  df <- ex_gmm(Ts=Ts,  num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, CI=TRUE, boots=2, mixture=TRUE)


  expect_no_error(discot <- DiSCoTEA(disco, agg="quantile"))
  expect_true(typeof(discot$plot) == "list")
  expect_no_error(summary(discot))
  expect_true(!is.null(discot$treats))
  expect_true(!is.null(discot$ses))
  expect_true(!is.null(discot$ci_lower))
  expect_true(!is.null(discot$ci_upper))

})


test_that("cdf works", {

  # disco
  Ts <- 2
  t0 <- 2
  df <- ex_gmm(Ts=Ts,  num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1,  CI=TRUE, boots=2)

  expect_no_error(discot <- DiSCoTEA(disco, agg="cdf"))
  expect_true(typeof(discot$plot) == "list")
  expect_no_error(summary(discot))
  expect_true(!is.null(discot$treats))
  expect_true(!is.null(discot$ses))
  expect_true(!is.null(discot$ci_lower))
  expect_true(!is.null(discot$ci_upper))


  # mixture
  Ts <- 2
  t0 <- 2
  df <- ex_gmm(Ts=Ts,  num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1,  CI=TRUE, boots=2, mixture=TRUE)

  expect_no_error(discot <- DiSCoTEA(disco, agg="cdf"))
  expect_true(typeof(discot$plot) == "list")
  expect_no_error(summary(discot))
  expect_true(!is.null(discot$treats))
  expect_true(!is.null(discot$ses))
  expect_true(!is.null(discot$ci_lower))
  expect_true(!is.null(discot$ci_upper))

})


test_that("quantileDiff works", {

  # disco
  Ts <- 2
  t0 <- 2
  df <- ex_gmm(Ts=Ts,  num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1,  CI=TRUE, boots=2)

  expect_no_error(discot <- DiSCoTEA(disco, agg="quantileDiff"))
  expect_true(typeof(discot$plot) == "list")
  expect_no_error(summary(discot))
  expect_true(!is.null(discot$treats))
  expect_true(!is.null(discot$ses))
  expect_true(!is.null(discot$ci_lower))
  expect_true(!is.null(discot$ci_upper))


  # mixture
  Ts <- 2
  t0 <- 2
  df <- ex_gmm(Ts=Ts,  num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1,  CI=TRUE, boots=2, mixture=TRUE)

  expect_no_error(discot <- DiSCoTEA(disco, agg="quantileDiff"))
  expect_true(typeof(discot$plot) == "list")
  expect_no_error(summary(discot))
  expect_true(!is.null(discot$treats))
  expect_true(!is.null(discot$ses))
  expect_true(!is.null(discot$ci_lower))
  expect_true(!is.null(discot$ci_upper))

})


test_that("cdfDiff works", {

  # disco
  Ts <- 2
  t0 <- 2
  df <- ex_gmm(Ts=Ts,  num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, CI=TRUE, boots=2)

  expect_no_error(discot <- DiSCoTEA(disco, agg="quantileDiff"))
  expect_true(typeof(discot$plot) == "list")
  expect_no_error(summary(discot))
  expect_true(!is.null(discot$treats))
  expect_true(!is.null(discot$ses))
  expect_true(!is.null(discot$ci_lower))
  expect_true(!is.null(discot$ci_upper))


  # mixture
  Ts <- 2
  t0 <- 2
  df <- ex_gmm(Ts=Ts,  num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, CI=TRUE, boots=2, mixture=TRUE)

  expect_no_error(discot <- DiSCoTEA(disco, agg="quantileDiff"))
  expect_true(typeof(discot$plot) == "list")
  expect_no_error(summary(discot))
  expect_true(!is.null(discot$treats))
  expect_true(!is.null(discot$ses))
  expect_true(!is.null(discot$ci_lower))
  expect_true(!is.null(discot$ci_upper))

})

test_that("samples works", {
  Ts <- 2
  t0 <- 2
  df <- ex_gmm(Ts=Ts,  num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, CI=TRUE, boots=2)

  smpls <- c(0.1, 0.2, 0.5, 0.9)
  expect_no_error(discot <- DiSCoTEA(disco, samples=smpls))
  expect_equal(discot$agg_df$X_from, c(0, smpls))
  expect_equal(discot$agg_df$X_to, c(smpls, 1))

})

test_that("nonsensical samples throws error", {

  Ts <- 2
  t0 <- 2
  df <- ex_gmm(Ts=Ts,  num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, CI=TRUE, boots=2)

  smpls <- c(0.5, 100)
  expect_error(discot <- DiSCoTEA(disco, samples=smpls))

})

test_that("graph to FALSE works", {
  Ts <- 2
  t0 <- 2
  df <- ex_gmm(Ts=Ts,  num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, CI=FALSE, boots=2)

  expect_no_error(discot <- DiSCoTEA(disco, graph=FALSE))

})



test_that("t_plots works", {
  Ts <- 4
  t0 <- 2
  df <- ex_gmm(Ts=Ts,  num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, CI=FALSE, boots=2)

  expect_no_error(discot <- DiSCoTEA(disco, t_plot=c(4)))


})



