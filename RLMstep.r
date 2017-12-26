RLMstep <- function(yinput, xinput,f.to.enter = 2,f.to.leave=1,bca.alpha = .95,is.bca = F){
## This function is based on the academic paper Robust Stepwise by Claudio Agostinelli
## The goal of this function is to preform hybrid selection for weighted least squares 
  if(is.bca){
    require(boot)
  }
  if(bca.alpha >= 1 || bca.alpha <= 0){
    geterrmessage("Pick bca.alpha between 0 and 1. This alpha is used to construct BCA conf intervals")
  }
  if(f.to.enter <= f.to.leave){
    geterrmessage("Please pick f.to.enter > f.to.leave")
  }
  if(xinput[,1] == 1){
    xinput<-xinput[,-1]
  }
  model.list <- list()
  max.num.pred <- ncol(xinput)    # maximum number of predictors in full model = # of predictors
  n <- nrow(xinput)
  error.final.models <- rep(NA,max.num.pred+1)     # error terms for each "Best" model with p predictors. +1 to include model with just intercept
  col.ind.in.data <- NULL    # tracks which columns ARE in model, INITIALLY NULL 
  col.ind.not.data <- 2:(max.num.pred+1)     # tracks which columns are NOT in model, INITIALLY EVERYTHING OTHER THEN ITERCEPT
  data.in.model <- NULL     # matrix that is BUILT UP with each iteration as variables are ADDED to the model. INITIALLY NULL
  
  # We initially calculate the intercept model
  rlm.fit <- rlm(yinput ~ 1) 
  error.final.models[1] <-rss(rlm.fit)
  model.list[[1]] <- rlm.fit # Include as first model
  col.ind.in.data <- append(col.ind.in.data,1) # Include in the data in model. It is already excluded from col.not.
  xinput <-cbind(1,xinput)
  data.in.model <- xinput[,1]
  #print(data.in.model)
  
  #We begin the selection process
  for(i in 1:max.num.pred){
    didAdd <- FALSE     # bool for adding variable
    error.add.jth.var <-  NULL      # Max size is everything NOT in data null since we are starting with forward
    error.remove.jth.var <- NULL  # Max size of elements we may remove from the model
    model.compare.rss <- error.final.models[i] # The model which we will use for F test comparison
    old.col.ind.data <- col.ind.in.data
    #FORWARD SELECTION
    weight1 <- NULL # Define a weight matrix
    for (j in col.ind.not.data) {
      model.build.data <- cbind(xinput[,col.ind.in.data], xinput[,j])     # loop through and create a dataset with last iteration PLUS new variable to test
      rlm.fit <- rlm(yinput ~ model.build.data-1) # Fit the model
      weight1 <- cbind(weight1,rlm.fit$w) # Gather weights
      error.add.jth.var<- append(error.add.jth.var,rss(rlm.fit)) # Append the errors to the vector.
    }
    
    ind.smallest <- which.min(error.add.jth.var)#give me the index of the variable which caused the biggest reduction in RSS
    min1 <- min(error.add.jth.var) # Find the min RSS
    test.1 <- model.forward.ratio(model.compare.rss,min1,length(col.ind.in.data),weight1[,ind.smallest]) # F test described by Agostinelli
    
    
    #print(c("Test Result Forward", test.1))
    if(test.1 > f.to.enter){ # If the folowing condition holds we add the model
      var.to.add <- col.ind.not.data[ind.smallest]  # Give me the actual variable 
      col.ind.in.data <-append(col.ind.in.data,var.to.add) #Append new variable
      col.ind.not.data <- col.ind.not.data[!col.ind.not.data==var.to.add] # Remove from not list 
      data.in.model <- xinput[,col.ind.in.data,drop=FALSE] # Input the variable
      model.compare.rss <- min1 # Record new comparison rss
      didAdd <- TRUE # We added
    }
    #Backward Selection
    weight2 <- NULL #Define a second weight matrix
    if(i>1 & didAdd){ #Entery Conditions
      for (j in old.col.ind.data){ 
        model.build.data <- xinput[,col.ind.in.data[col.ind.in.data!=j]] # Iterativly remove predictors
        rlm.fit <- rlm(yinput ~ model.build.data-1)  # Regress them
        weight2 <- cbind(weight2,rlm.fit$w) # Retrieve rlm weights
        error.remove.jth.var <-append(error.remove.jth.var,rss(rlm.fit)) # store rss
      }  
      min2 <- min(error.remove.jth.var) # Min RSS
      ind.min2 <- which.min(error.remove.jth.var) # Index min 
      test.2 <- model.backward.ratio(min2,min1,length(old.col.ind.data),weight2[,ind.min2]) # Agostelli test 
      #print(c("Test Result Backward",test.2))
      if( test.2 < f.to.leave ){ # If the following holds remove
        var.to.remove <- old.col.ind.data[ind.min2] # find var to remove
        col.ind.in.data <-col.ind.in.data[col.ind.in.data != var.to.remove] # remove it from col.in.data
        col.ind.not.data<- append(col.ind.not.data,var.to.remove) # add it to col.not
        data.in.model <- xinput[,col.ind.in.data,drop=FALSE] # update data.in.model
        #Exist condition described by agostinelli to ensure convergence
        if(min2<min1*(1+f.to.leave/(sum(weight2)-length(old.col.ind.data)))/(1+f.to.enter/(sum(weight2))-length(old.col.ind.data)-1)){break} 
      }
    }
    rlm.fit <- rlm(yinput ~ data.in.model-1)     #refit model, store it, and extract needed error term  
    model.list[[i+1]] <- rlm.fit ## Figure out where to store
    error.final.models[i+1] <- rss(rlm.fit)
    if(!didAdd){break}
    else{if(i>2)if((test.1 <= f.to.enter) & (test.2 >= f.to.leave)){break}}
  }
  
  
  
  rlm.final.model <- rlm(yinput ~ data.in.model-1) # Compute the final model
  regression.fucntion <- function(formula, data, index){ # Boot helper function for BCa
    rlm.new <- rlm(formula,data=data[index,])
    return(coef(rlm.new))
  }
  if(is.bca){
    # Compute the BCa confidence intervals
    data.frame.1 <- cbind(data.frame(yinput),data.frame(data.in.model))
    boot.rlm <- boot(data = data.frame.1,statistic =regression.fucntion ,R = 10000,formula =yinput ~ .-1)
    bca<-lapply(1:(ncol(data.in.model)),function(x)return(boot.ci(boot.out =boot.rlm,index = x,type ="bca",conf = bca.alpha)))
  }
  # Store the results
  df.temp.model.results <- list()
  df.temp.model.results$predictor <- data.in.model
  df.temp.model.results$predictor.names <- colnames(data.in.model)
  if(is.bca){df.temp.model.results$bca <- bca}else{df.temp.model.results$bca <- "Set is.bca = TRUE"}
  df.temp.model.results$coef <- rlm.final.model$coef
  return(df.temp.model.results)
  
  #df.model.results <<- rbind(df.model.results, df.temp.model.results)
}

# Helper functions
model.forward.ratio <- function(rss.a,rss.a1,da,w){
  return((rss.a-rss.a1)/(rss.a1/(sum(w)-da-1)))
}
model.backward.ratio <- function(rss.a,rss.a1,da,w){
  return((rss.a-rss.a1)/(rss.a1/(sum(w)-da)))
}
rss <- function(rl){return(sum(rl$w*rl$resid^2) )}

