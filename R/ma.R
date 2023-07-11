# This is the interface to use the earth() function 
# within gamlss()
# require(earth)
# Author Vlasios Voudouris
# latest change 22 September 2015
ma <- function(formula,control=ma.control(...),...) 
{ 
  #------------------------------------------
  # function starts here
  #------------------------------------------
  scall=gsub("  ", "", Reduce(paste, deparse(sys.call(),width.cutoff = 500L)))
  if (!is(formula, "formula")) stop("formula argument in ma() needs a formula starting with ~")
  # get where "gamlss" is in system call
  # it can be in gamlss() or predict.gamlss()  
  rexpr <- grepl("gamlss",sys.calls()) ## 
  for (i in length(rexpr):1)
  { 
    position <- i # get the position
    if (rexpr[i]==TRUE) break
  }
  # 
  gamlss.env <- sys.frame(position) #gamlss or predict.gamlss
  ##---
  ## get the data
  if (sys.call(position)[1]=="predict.gamlss()")
  { # if predict is used 
    Data <- get("data", envir=gamlss.env)
  }
  else { # if gamlss() is used
    #stop("the option data in gamlss() is required for lo() to work")
    if (is.null(get("gamlsscall", envir=gamlss.env)$data)) 
    { # if no data argument but the formula can be interpreted
      Data <- model.frame(formula)	
    }
    else
    {# data argument in gamlss 
      Data <- get("gamlsscall", envir=gamlss.env)$data
    }
  }
     Data <- data.frame(eval(substitute(Data)))
  #===== 
      len <- dim(Data)[1] # get the lenth of the data
  ## out
     xvar <- rep(0,  len) # model.matrix(as.formula(paste(paste("~",ff, sep=""), "-1", sep="")), data=Data) #
  attr(xvar,"formula")     <- formula
   attr(xvar,"control")    <- control
  attr(xvar, "gamlss.env") <- gamlss.env
  attr(xvar, "data")       <- as.data.frame(Data)  
  attr(xvar, "call")       <- substitute(gamlss.ma(data[[scall]], z, w, ...)) 
  attr(xvar, "class")      <- "smooth"
  xvar
}
#*********************************************************************
#*********************************************************************
ma.control <-function ( wp = NULL, 
                   pmethod = "backward",
                    keepxy = FALSE, 
                     trace = 0,  
                    degree = 1,
                    nprune = NULL,
                    ncross = 1, 
                     nfold = 0, 
                  stratify = TRUE,
                   Scale.y = 1,... ) 
{
  
  list(wp=wp, pmethod=pmethod, keepxy=keepxy,trace=trace, degree=degree,
       nprune=nprune, ncross=ncross, nfold=nfold, stratify =stratify, 
       Scale.y =Scale.y ,...)
}
#****************************************************************
#****************************************************************
gamlss.ma <- function(x, y, w, xeval = NULL, ...)
{      	      
      form <- attr(x,"formula")
  formtemp <- gsub("  ", "", Reduce(paste, deparse(form,width.cutoff = 500L)))
      form <- as.formula(paste("Y.var",formtemp, sep=""))
   control <- as.list(attr(x, "control"))  
  #gamlss.env <- as.environment(attr(x, "gamlss.env"))
     OData <- attr(x,"data") 
      Data <- if (is.null(xeval)) OData #the trick is for prediction
             else  OData[seq(1,length(y)),]
    Y.var <- y
    W.var <- w
     Data <- data.frame(eval(substitute(Data)),Y.var,W.var)
  # 
      fit <- earth(form, data=Data, weights=W.var, 
                    wp = control$wp,
               pmethod = control$pmethod, 
                keepxy = control$keepxy,
                 trace = control$trace, 
                degree = control$degree,
                nprune = control$nprune,
                ncross = control$ncross,
                 nfold = control$nfold,
              stratify = control$stratify,
               Scale.y = control$Scale.y)
 #earth(f   wp = NULL, subset = NULL,
       df <- length(coef(fit))
       fv <- fitted(fit) 
residuals <- y-fv #resid(fit)
  if (is.null(xeval)) #proper fit
  {
    list(fitted.values=fv, residuals=residuals,
         nl.df = df, lambda=NA, #
         coefSmo = fit, var=hatvalues(fit))    # var=fv has to fixed
  }
  else #predict
  {
    ll<-dim(OData)[1]
    pred <- predict(fit,newdata = OData[seq(length(y)+1,ll),])
  }         
}
