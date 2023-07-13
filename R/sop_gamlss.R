#require(SOP)
################################################################################
################################################################################
################################################################################
################################################################################
so <- SOP <- function (formula, control = sop.control(...), ...) 
{
#------------------------------------------
# function starts here
#------------------------------------------   
  scall <- deparse(sys.call(), width.cutoff = 500L)
  if (!is(formula, "formula")) 
    stop("formula argument in SOP() needs a formula starting with ~")
# get where "gamlss" is in system call, it can be in gamlss() or predict.gamlss()
  rexpr <- grepl("gamlss", sys.calls())
  for (i in length(rexpr):1) {
    position <- i
    if (rexpr[i] == TRUE) 
      break
  }
  gamlss.env <- sys.frame(position) #gamlss or predict.gamlss
## get the data
## this has been modified on the 12-12-14 to make sure that 
##  if model.frame.gamlss() is used as for example in term.plot() the
## function does not fail (It need to implemented to all smoother using formulae?)  
  if (sys.call(position)[1] == "predict.gamlss()") {# if predict is used 
    Data <- get("data", envir = gamlss.env)
  }
  else if (sys.call(position)[1] == "gamlss()") {# if gamlss() is used
    if (is.null(get("gamlsscall", envir = gamlss.env)$data)) { # if no data argument but the formula can be interpreted
      Data <- model.frame(formula)
    }
    else {# data argument in gamlss
      Data <- get("gamlsscall", envir = gamlss.env)$data
    }
  }
  else {
    Data <- get("data", envir = gamlss.env)
  }
  Data <- data.frame(eval(substitute(Data), envir=gamlss.env)) #<--- changed here
#----------------------------------------------
     fOrmula <- as.formula(paste0("Y.var", deparse(formula, width.cutoff = 500L)))
  Data$Y.var <- rep(0, nrow(Data))
           G <- sop(fOrmula, data=Data, fit = FALSE)
#----------------------------------------------
# get a random name to use it in the gamlss() environment
#--------
             sl <- sample(letters, 4)
    fourLetters <- paste(paste(paste(sl[1], sl[2], sep=""), sl[3], sep=""),sl[4], sep="")
startLambdaName <- paste("start.Lambda",fourLetters, sep=".")
           ## put the starting values in the gamlss()environment
           #--------
assign(startLambdaName, rep_len(1,length(G$G)+1), envir=gamlss.env)
           #--------           
                       xvar <- rep(0, dim(Data)[1])
      attr(xvar, "formula") <- formula
  #   attr(xvar, "control") <- control
   attr(xvar, "gamlss.env") <- gamlss.env
attr(xvar, "NameForLambda") <- startLambdaName
   attr(xvar, "gamlss.env") <- gamlss.env
         attr(xvar, "data") <- as.data.frame(Data)
         attr(xvar, "call") <- substitute(gamlss.so(data[[scall]], z, w, ...))
        attr(xvar, "class") <- "smooth"
            attr(xvar, "G") <- G
  xvar
}
################################################################################
################################################################################
################################################################################
################################################################################
gamlss.so  <- function (x, y, w, xeval = NULL, ...) 
{
  if (is.null(xeval)) {
            Y.var <- y
            W.var <- w
                G <- attr(x, "G")
         G$data$y <- y
         G$data$w <- w
  startLambdaName <- as.character(attr(x, "NameForLambda"))
       gamlss.env <- as.environment(attr(x, "gamlss.env"))
          lambdaS <- get(startLambdaName, envir=gamlss.env)
              fit <- sop.fit(X = G$X, Z = G$Z, G = G$G, y = G$data$y, 
                                  weights = G$data$w, vcstart=lambdaS,
                             control=sop.control(trace=FALSE))
  assign(startLambdaName, fit$out$vc, envir=gamlss.env)   
               fv <- fitted(fit)
        residuals <- y - fv
          fit$lin <- G$lin
       fit$random <- G$random
            fit$f <- G$f    
            fit$G <- G$G
      fit$formula <- G$formula
         fit$vars <- G$vars
        fit$terms <- G$terms
    fit$na.action <- G$na.action 
  fit$names.terms <- G$names.terms
  fit$model.terms <- G$model.terms
         fit$data <- G$data
           fit$gg <- G$gg
           fit$np <- G$np
       fit$nterms <- G$nterms
              lvc <- length(fit$out$vc)
           lambda <- fit$out$vc[1]/fit$out$vc[2:lvc]
              edf <- sum(fit$out$edf)+length(fit$b.fixed)-1
list(fitted.values = fv, residuals = residuals, nl.df = edf, 
            lambda = lambda, coefSmo = fit, 
               var = NA)
  }
  else {
    gamlss.env <- as.environment(attr(x, "gamlss.env"))
            obj <- get("object", envir = gamlss.env)
             TT <- get("TT", envir = gamlss.env)
             SL <- get("smooth.labels", envir = gamlss.env)
            fit <- eval(parse(text = paste("obj$", get("what", envir = gamlss.env), ".coefSmo[[", as.character(match(TT, SL)), "]]", sep = "")))
         OData <- attr(x, "data")
            ll <- dim(OData)[1]
          pred <- predict(fit, newdata = OData[seq(length(y) + 1, ll), ])
  }
}
################################################################################
################################################################################
################################################################################
################################################################################
