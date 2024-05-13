## -----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# test each estimation method with sample of
# Expected results: treatment effects = 1, p-value for pre-test
# uniformly distributed, ipw model is incorectly specified here
#-----------------------------------------------------------------------------


test_that("mixture weights sum up to 1", {
  Ts <- 2
  t0 <- 2
  df <- ex_gmm(Ts=Ts,  num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, CI=TRUE, boots=2, mixture=TRUE, num.cores=1)
  # period-specific weights
  for (t in 1:(t0-1)) {
    expect_equal(sum(disco$results.periods[[t]]$mixture$weights), 1, tolerance=1e-3)
  }
  expect_equal(sum(disco$weights), 1, tolerance=1e-3)

})




#-----

test_that("disco weights sum up to 1", {

  ## test for more than pre- and post-period
  Ts <- 5
  t0 <- 3
  df <- ex_gmm(Ts=Ts,  num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, num.cores=1)

  # period-specific weights
  for (t in 1:(t0-1)) {
    expect_equal(sum(disco$results.periods[[t]]$DiSCo$weights), 1)
  }

  # overall weights
  expect_equal(sum(disco$weights), 1)

  ## test for only pre- and post-period
  Ts <- 2
  t0 <- 2
  df <- ex_gmm(Ts=Ts, num.con=4)
  disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, num.cores=1)

  # period-specific weights
  for (t in 1:(t0-1)) {
    expect_equal(sum(disco$results.periods[[t]]$DiSCo$weights), 1)
  }

  # overall weights
  expect_equal(sum(disco$weights), 1)

})



test_that("wrong treatment period throws error", {

  Ts <- 2
  t0 <- 1
  df <- ex_gmm(Ts=Ts, num.con=4)
  expect_error(DiSCo(df=df, id_col.target=1, t0=t0, seed=1))

  Ts <- 2
  t0 <- 3
  df <- ex_gmm(Ts=Ts)
  expect_error(DiSCo(df=df, id_col.target=1, t0=t0, seed=1))

})

test_that("simplex results in weakly positive weights", {

  Ts <- 2
  num.con <- 4
  t0 <- 2
  df <- ex_gmm(Ts=Ts, num.con=num.con)

  # test simplex=TRUE
  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, num.cores=1, simplex=TRUE))
  # expect true up to some margin of error
  expect_true(all(disco$weights > -1e-10))
})


test_that("simplex results in weakly positive weights for mixture", {

  Ts <- 2
  num.con <- 4
  t0 <- 2
  df <- ex_gmm(Ts=Ts, num.con=num.con)

  # test simplex=TRUE
  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, num.cores=1, simplex=TRUE, mixture=TRUE))
  # expect true up to some margin of error
  expect_true(all(disco$weights > -1e-5))
})


test_that("returned quantile functions and cdfs are weakly increasing", { # TODO left off here

  ## test mixture
  Ts <- 2
  num.con <- 4
  t0 <- 2
  df <- ex_gmm(Ts=Ts, num.con=num.con)

  # test simplex=TRUE
  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, num.cores=1, simplex=TRUE, mixture=TRUE))

  t_cdf <- disco$results.periods[[1]]$target$cdf
  expect_true(all(t_cdf[-1] - t_cdf[-length(t_cdf)] > -1e-5))

  disco_cdf <-  disco$results.periods[[1]]$DiSCo$cdf
  expect_true(all(disco_cdf[-1] - disco_cdf[-length(disco_cdf)] > -1e-5))

  ## test disco
  Ts <- 2
  num.con <- 4
  t0 <- 2
  df <- ex_gmm(Ts=Ts, num.con=num.con)

  # test simplex=TRUE
  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, num.cores=1, simplex=TRUE, mixture=FALSE))

  t_cdf <- disco$results.periods[[1]]$target$quantiles
  expect_true(all(t_cdf[-1] - t_cdf[-length(t_cdf)] > -1e-5))

  disco_cdf <-  disco$results.periods[[1]]$DiSCo$quantile
  expect_true(all(disco_cdf[-1] - disco_cdf[-length(disco_cdf)] > -1e-5))


})


test_that("alternative qmethods work", {
  Ts <- 2
  t0 <- 2
  num.con <- 4
  df <- ex_gmm(Ts=Ts, num.con=num.con)

  # qmethod = "qkden"
  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, qmethod="qkden", seed=1, num.cores=1))

  # qmethod = "extreme"
  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, qmethod="extreme", seed=1, num.cores=1))
})


test_that("discrete support works for mixture", {
  Ts <- 2
  t0 <- 2
  num.con <- 4
  df <- ex_gmm(Ts=Ts, num.con=num.con)

  df[, y_col := round(y_col)]
  sprt <- unique(df$y_col)


  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0,  seed=1, num.cores=1, mixture=TRUE, grid.cat=sprt))

  # check that the length of the cdf is the same as the length of the support
  expect_true(length(disco$results.periods[[1]]$DiSCo$cdf) == length(sprt))

  # check that disco quantile function is weakly increasing
  t_cdf <- disco$results.periods[[1]]$DiSCo$quantile
  expect_true(all(t_cdf[-1] - t_cdf[-length(t_cdf)] > -1e-5))

})


# test that all arguments of DiSCo work
test_that("test that variations of other arguments work",  {
  Ts <- 2
  t0 <- 2
  num.con <- 4
  df <- ex_gmm(Ts=Ts, num.con=num.con)

  # M too small
  M <- num.con - 1
  expect_error(DiSCo(df=df, id_col.target=1, t0=t0, M=M, seed=1))

  # M just large enough
  M <- num.con
  expect_no_error(DiSCo(df=df, id_col.target=1, t0=t0, M=M, seed=1))

  # G too small
  G <- 1
  expect_error(DiSCo(df=df, id_col.target=1, t0=t0, G=G, seed=1))

  # G just large enough
  G <- 2
  expect_no_error(DiSCo(df=df, id_col.target=1, t0=t0, G=G, seed=1))

  # id_col.targt not in df
  id_col.target <- max(df$id_col) + 1
  expect_error(DiSCo(df=df, id_col.target=id_col.target, t0=t0, seed=1))

  # # num.cores too large
  # num.cores <- 100
  # expect_error(DiSCo(df=df, id_col.target=1, t0=t0, num.cores=num.cores, seed=1))

  # # parallel cores
  # num.cores <- 1
  # expect_no_error(DiSCo(df=df, id_col.target=1, t0=t0, num.cores=num.cores, seed=1))

  # permutation TRUE
  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, permutation=TRUE, seed=1))
  expect_true(!is.null(disco$perm))


  # permutation TRUE and mixture TRUE
  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, permutation=TRUE, seed=1, mixture=TRUE))
  expect_true(!is.null(disco$perm))


  # permutation TRUE and mixture TRUE and grid.cat
  temp <- df
  temp[, y_col := round(y_col)]
  sprt <- unique(temp$y_col)
  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, permutation=TRUE, seed=1, mixture=TRUE, grid.cat=sprt))
  expect_true(!is.null(disco$perm))

  # permutation FALSE
  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, permutation=FALSE, seed=1))
  expect_true(is.null(disco$perm))

  # q_min too small
  q_min <- -1
  expect_error(DiSCo(df=df, id_col.target=1, t0=t0, q_min=q_min, seed=1))

  # q_max too large
  q_max <- 2
  expect_error(DiSCo(df=df, id_col.target=1, t0=t0, q_max=q_max, seed=1))

  # q_min larger than q_max
  q_min <- 0.8
  q_max <- 0.2
  expect_error(DiSCo(df=df, id_col.target=1, t0=t0, q_min=q_min, q_max=q_max, seed=1))

  # q_min and q_max normal range
  q_min <- 0.2
  q_max <- 0.8
  expect_no_error(DiSCo(df=df, id_col.target=1, t0=t0, q_min=q_min, q_max=q_max, seed=1))

  # CI FALSE
  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, CI=FALSE, seed=1))
  expect_true(is.null(disco$CI$quantile))

  # CI without replacement
  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, CI=TRUE, replace=FALSE, boots=2, seed=1))
  expect_true(!is.null(disco$CI$quantile))

  # boots 0
  expect_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, boots=0, seed=1, CI=TRUE, num.cores=1))


  # cl = 0
  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, cl=0, seed=1, CI=TRUE, num.cores=1, boot=2))
  expect_true(all(disco$CI$quantile$lower == disco$CI$quantile$upper))

  # cl < 0
  expect_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, cl=-1, seed=1, CI=TRUE, num.cores=1, boot=2))

  # cl > 1
  expect_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, cl=2, seed=1, CI=TRUE, num.cores=1, boot=2))


  # test seed
  expect_no_error(disco <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, num.cores=1))
  expect_no_error(disco2 <- DiSCo(df=df, id_col.target=1, t0=t0, seed=1, num.cores=1))
  expect_equal(disco$weights, disco2$weights)


})


