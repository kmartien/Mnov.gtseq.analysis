# 18March2025 This is not working as expected. I can't use predict for a single
# sample, which is what I need in order to do LOOCV. Worse, when I use predict
# with multiple samples, the result for each sample varies depending on what
# other samples it is predicted with. That shouldn't be the case.

# num.ind <- getNumInd(g.stratified)
# cv.ids <- sample(
#   1:num.ind,
#   size = 0.1 * num.ind
# )
# cv.genind <- genind.strat[cv.ids]
# train.genind <- genind.strat[-cv.ids]
# n.da <- 8
# n.pca <- 100


fitTrainDAPC <- function(train.genind, n.da, n.pca){
  fit <- dapc(
    x = train.genind, 
    n.da = n.da, 
    n.pca = n.pca
    )
}

predictTestDAPC <- function(fit, cv.genind){
  pred <- predict.dapc.loov(
    fit, 
    newdata = cv.genind
    )
  
  tibble(
    swfsc.id = dimnames(tab(cv.genind))[[1]],
    obs.stratum = pop(cv.genind),
    pred.stratum = pred$assign
  )
}

predictAllIDsDAPC <- function(genind.all, n.da, n.pca){
  lapply(1:length(pop(genind.all)), function(cv.ind){
    fitTrainDAPC(train.genind = genind.all[-cv.ind], n.da = n.da, n.pca = n.pca) |> 
      predictTestDAPC(cv.genind = genind.all[cv.ind])
  }) |> 
    bind_rows()
}

predict.dapc.loov <- function (object, newdata, prior = object$prior, dimen, method = c("plug-in", 
                                                                                        "predictive", "debiased"), ...) 
{
  if (!inherits(object, "dapc")) 
    stop("x is not a dapc object")
  method <- match.arg(method)
  x <- as.lda(object)
  if (!missing(newdata)) {
    if (is.null(object$pca.loadings)) 
      stop("DAPC object does not contain loadings of original variables. \nPlease re-run DAPC using 'pca.loadings=TRUE'.")
    if (is.genind(newdata)) {
      newdata <- matrix(tab(newdata, freq = TRUE, NA.method = "mean"), nrow = 1) #this is the line I had to change
    }
    else if (inherits(newdata, "genlight")) {
      newdata <- as.matrix(newdata)/ploidy(newdata)
    }
    else {
      newdata <- as.matrix(newdata)
    }
    if (ncol(newdata) != nrow(object$pca.loadings)) 
      stop("Number of variables in newdata does not match original data.")
    for (i in 1:nrow(newdata)) {
      newdata[i, ] <- (newdata[i, ] - object$pca.cent)/object$pca.norm
    }
    newdata[is.na(newdata)] <- 0
    XU <- newdata %*% as.matrix(object$pca.loadings)
  }
  else {
    XU <- object$tab
  }
  colnames(XU) <- colnames(object$tab)
  if (!missing(dimen)) {
    if (dimen > object$n.da) 
      stop(paste("Too many dimensions requested. \nOnly", 
                 object$n.da, "discriminant functions were saved in DAPC."))
  }
  else {
    dimen <- object$n.da
  }
  temp <- predict(x, XU, prior, dimen, method, ...)
  res <- list()
  res$assign <- temp$class
  res$posterior <- temp$posterior
  res$ind.scores <- temp$x
  return(res)
}
