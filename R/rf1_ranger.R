################################################################################
################################################################################
################################################################################
################################################################################
#require(ranger)
################################################################################
################################################################################
################################################################################
################################################################################
rf1 <- function(formula, ...){
	control <- as.list(match.call())
	#---------------------------------
	scall <- deparse(sys.call(), width.cutoff = 500L) # check the formula
	if (!is(formula, "formula")) stop("formula argument in ranger() needs a formula starting with ~")
	# get where "gamlss" is in system call
	# it can be in gamlss() or predict.gamlss()
	rexpr <- grepl("gamlss", sys.calls()) ##
for (i in length(rexpr):1) 
  {
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
      attr(xvar, "call") <- substitute(gamlss.rf1(data[[scall]], z, w, ...))
     attr(xvar, "class") <- "smooth"
	xvar
}
################################################################################
################################################################################
################################################################################
################################################################################
# the definition of the backfitting additive function
gamlss.rf1 <- function(x, y, w, xeval = NULL, ...) {
 	  formula <- attr(x, "formula")
  	formula <- as.formula(paste("y", deparse(formula, width.cutoff = 500L), 
  	                            sep = ""))
	  control <- as.list(attr(x, "control"))
 gamlss.env <- as.environment(attr(x, "gamlss.env"))
     	OData <- attr(x, "data")
       Data <- if (is.null(xeval)) 
         {
		 OData # this is for prediction
	       } else {
		OData[seq(1, length(y)), ]
	       }
    	Data <- data.frame(eval(substitute(Data)), y, w)
	   rexpr <- regexpr("gamlss", sys.calls())
	      nt <- if(is.null(control$num.trees))
	        {
		        500 # Uses 500 as a default value, like ranger does
	        } else 
	        {
		control$num.trees
	        }
importance <- if(is.null(control$importance))
          {
	          NULL # Uses 500 as a default value, like ranger does
        	} else 
        	{
	  match.arg(control$importance, choices=c('none', 'impurity', 'impurity_corrected', 'permutation'))
	        }
# needs to control options 
	fit <- ranger::ranger(formula, data = Data, case.weights = w, num.trees = nt,
		write.forest = TRUE, oob.error = FALSE, seed = 1, 
		regularization.factor=1, importance=importance)
	# mtry how many x's to try
	# min.node.size 5 for regression
	#  write.forest 
	#  probability
	#  min.node.size 5 for regression
	#  max.depth
	#  # min.bucket Maximal tree depth. A value of NULL or 0 (the default) 
	#            corresponds to unlimited depth, 1 to tree stumps 
	#            (1 split per tree).
	# replace = TRUE
	# sample.fraction = ifelse(replace, 1, 0.632)
	# case.weights = NULL
	# class.weights = NULL
	# splitrule = NULL,
	# num.random.splits = 1,
	# alpha = 0.5,
	# minprop = 0.1,
	# split.select.weights = NULL,
	# always.split.variables = NULL,
	# respect.unordered.factors = NULL
	# Unordered factor covariates can be handled in 3 different ways by using                 respect.unordered.factors: For 'ignore' all factors are regarded 
	# ordered, for 'partition' all possible 2-partitions are considered for 
	#  splitting. For 'order' and 2-class classification the factor levels are
	#   ordered by their proportion falling in the second class, 
	#   for regression by their mean response, as described in 
	#   Hastie et al. (2009), chapter 9.2.4. For multiclass classification the
	#    factor levels are ordered by the first principal component of the weighted
	#     covariance matrix of the contingency table (Coppersmith et al. 1999), 
	#     for survival by the median survival (or the largest available quantile 
	#     if the median is not available). The use of 'order' is recommended, 
	#     as it computationally fast and can handle an unlimited number of factor
	#      levels. Note that the factors are only reordered once and not again in
	#       each split.
	# scale.permutation.importance = FALSE,
	# local.importance = FALSE,
	# Regularization works by penalizing new variables by multiplying the 
	# splitting criterion by a factor, see Deng & Runger (2012) for details. 
	# If regularization.usedepth=TRUE, is used, where f is the regularization 
	# factor and d the depth of the node. If regularization is used, 
	# multithreading is deactivated because all trees need access to the list 
	# of variables that are already included in the model.
	# regularization.factor = 1,
	# regularization.usedepth = FALSE,
	# keep.inbag = FALSE
	# inbag = NULL,
	# holdout = FALSE
	# quantreg = FALSE,
	# oob.error = TRUE,
	# num.threads = NULL,
	# save.memory = FALSE,
	# verbose = TRUE,
	# seed = NULL,
	# dependent.variable.name = NULL,
	# status.variable.name = NULL,
	# classification = NULL,
	# x = NULL,
	# y = NULL
if (is.null(xeval)) { # This is for fit/not-prediction
		df <-  sum(sapply(fit$forest$split.varIDs, function(x) sum(x)/nt))
	#	cat("df=", df, "\n")
		#  sum(sapply(fit$forest$split.varIDs, length))  # Number of break point - I wonder if this is the best way to calculate df
		fv <- predict(fit, data = Data, predict.all = FALSE, type = "response", seed = 1)$predictions
		residuals <- Data[,'y'] - fv # Raw residuals
		list(fitted.values = fv, residuals = residuals,
			nl.df = df, #, lambda = NA, # I'm not sure how to calculate df for this case
		  lambda = NA,
			coefSmo = fit, var = NA)
	} else { # this is for prediction/not-fit
		ndata <- subset(OData, source == "newdata")
		pred <- predict(fit, data = ndata, predict.all = FALSE, type = "response", seed = 1)$predictions
		pred
	}
}
################################################################################
################################################################################
################################################################################
################################################################################
