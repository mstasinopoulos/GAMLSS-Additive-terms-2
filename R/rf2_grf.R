################################################################################
################################################################################
################################################################################
################################################################################
# require(grf)
################################################################################
################################################################################
################################################################################
################################################################################
rf2 <- function(formula, ...){
	control <- as.list(match.call())
#---------------------------------
	scall <- deparse(sys.call(), width.cutoff = 500L) # check the formula
	if (!is(formula, "formula")) stop("formula argument in grf() needs a formula starting with ~")

	# get where "gamlss" is in system call
	# it can be in gamlss() or predict.gamlss()
	rexpr <- grepl("gamlss", sys.calls()) ##
	for (i in length(rexpr):1) {
		position <- i # get the position, we are geting the fist from the last
		if (rexpr[i] == TRUE) break
	}

	gamlss.env <- sys.frame(position) # gamlss or predict.gamlss
	if (sys.call(position)[1] == "predict.gamlss()") { # if predict is used
		Data <- get("data", envir = gamlss.env)
	} else if (sys.call(position)[1] == "gamlss()") { # if gamlss() is used
		if (is.null(get("gamlsscall", envir = gamlss.env)$data)) { # if no data argument but the formula can be interpreted
			Data <- model.frame(formula)
		} else { # data argument in gamlss
			Data <- get("gamlsscall", envir = gamlss.env)$data
		}
	} else {
		Data <- get("data", envir = gamlss.env)
	}

	Data <- data.frame(eval(substitute(Data)))

	# =====
	len <- dim(Data)[1] # get the length of the data
	xvar <- rep(0, len)
	attr(xvar, "formula") <- formula
	attr(xvar, "control") <- control
	attr(xvar, "gamlss.env") <- gamlss.env
	attr(xvar, "data") <- as.data.frame(Data)
	attr(xvar, "call") <- substitute(gamlss.rf2(data[[scall]], z, w, ...))
	attr(xvar, "class") <- "smooth"
	xvar
}
################################################################################
################################################################################
################################################################################
################################################################################
# the definition of the backfitting additive function
gamlss.rf2 <- function(x, y, w, xeval = NULL, ...) {
	formula <- attr(x, "formula")
	formula <- as.formula(paste("y", deparse(formula, width.cutoff = 500L), sep = ""))
	control <- as.list(attr(x, "control"))
	gamlss.env <- as.environment(attr(x, "gamlss.env"))

	OData <- attr(x, "data")
	Data <- if (is.null(xeval)) {
		OData # this is for prediction
	} else {
		OData[seq(1, length(y)), ]
	}

	Data <- data.frame(eval(substitute(Data)), y, w)
	rexpr <- regexpr("gamlss", sys.calls())

	nt <- if(is.null(control$num.trees)){ # The number of trees (num.trees)
		500 # Uses 500 as a default value, like ranger does
	} else {
		control$num.trees
	}

	if(length(formula) == 3){
		ks <- all.vars(formula[[3]])
	} else {
		ks <- all.vars(formula[[2]])
	}

	fit <- grf::regression_forest(X = Data[, ks], Y = Data[,'y'],
		sample.weights = w, num.trees = nt,
		seed = 1)

	if (is.null(xeval)) { # This is for fit/not-prediction
		df <- NA

		fv <- predict(fit, Data[,ks])$predictions

		residuals <- Data[,'y'] - fv # Raw residuals

		list(fitted.values = fv, residuals = residuals,
#			nl.df = df - 1, lambda = NA, # I'm not sure how to calculate df for this case
			nl.df = 1, lambda = NA,
			coefSmo = fit, var = NA)
	} else { # this is for prediction/not-fit
		ndata <- subset(OData, source == "newdata")
		pred <- predict(fit, ndata[,ks])$predictions
		pred
	}
}
################################################################################
################################################################################
################################################################################
################################################################################

