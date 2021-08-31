
#GENERATES LAGGED ENVIRONMENTAL DATA BEFORE EVERY FISH CATCHES SAMPLE
## Functions from Blas Benito (https://github.com/BlasBenito/BaleFire) https://zenodo.org/record/2599859#.YCTzbGhKiUk 

backwardLags <- function(lags, reference.data, data.to.lag){
  
  #df to store the lagged data
  lag.data = data.frame(biomass=double(), environment=double(), lag=integer())
  
  #iterates through erica lines
  for (catch.case in 1:nrow(catches)){
    
    #take a line of the erica dataframe and replicate it as many times as lags are
    catch.value = rep(catches[catch.case, "biomass"], max(lags))

    #get the age of the replicated line
    catch.case.year=catches[catch.case, "year"]
    catch.case.year.plus.lags=catch.case.year + (1 * max(lags))
    
    #if beyond maximum age
    if (catch.case.year.plus.lags > max(env$year)){break}
    
    #get from char.interpolated the lines with age > age.erica && age <= age.erica + lags
    env.temp=env[which(env$year > catch.case.year & env$year <= catch.case.year.plus.lags), "environment"]
    
    #put the data together
    env.temp=data.frame(biomass=catch.value, environment=env.temp, lag=lags)
    
    #put them in the final table
    lag.data=rbind(lag.data, env.temp)
    
  }#end of iterations
  
  #remove stuff we don't need
  rm(env.temp, catch.case, catch.case.year, catch.case.year.plus.lags, catch.value)
  
  #order by lag
  lag.data=lag.data[order(lag.data$lag),]
  
  #standardize data
  #get lags column
  lags.column=lag.data$lag
  
  #standardize
  lag.data=scale(lag.data[, c("biomass", "environment")])
  
  #add lag
  lag.data=data.frame(lag.data, lag=lags.column)
  
  #lag as factor
  lag.data$lag=as.factor(lag.data$lag)
  
  return(lag.data)
}


#COMPUTES GLS MODELS BETWEEN A RESPONSE AND A VARIABLE IN A LAGGED DATASET. PROVIDES A DATAFRAME
modelLagData <- function(model.formula, lagged.data){
  
  lags=as.numeric(sort(unique(lagged.data$lag)))
  model.formula=as.formula(model.formula)
  response = all.vars(model.formula)[1]
  
  #list to store results
  results.list = list()
  
  #fitting a model per lag
  for (lag in lags){
    results.list[[lag]] = gls(model.formula, data=lagged.data[lagged.data$lag==lag,])
  }
  
  #list to store pseudo R2
  results.list.R2=list()
  
  #computing pseudo R2 per lag
  for (lag in lags){
    results.list.R2[[lag]] = cor(lagged.data[lagged.data$lag==lag, response], predict(results.list[[lag]]))^2
  }
  
  #gathering coefficients
  results.list.coef = lapply(results.list, function(x){summary(x)$tTable[2,1]})
  
  #gathering the standard error of the coefficients
  results.list.coef.se = lapply(results.list, function(x){summary(x)$tTable[2,2]})
  
  #gathering p-value of coefficients
  results.list.pvalue = lapply(results.list, function(x){summary(x)$tTable[2,4]})
  
  #to data frame
  output.df = as.data.frame(do.call(rbind, results.list.coef))
  output.df = data.frame(lag=as.numeric(rownames(output.df)), Coefficient=output.df$V1)
  output.df[, "p-value"] = as.data.frame(do.call(rbind, results.list.pvalue))$V1
  output.df$R2 = as.data.frame(do.call(rbind, results.list.R2))$V1
  output.df$lag=output.df$lag

  #se of coefficients
  output.df.se = as.data.frame(do.call(rbind, results.list.coef.se))
  output.df.se$lower = output.df$Coefficient - output.df.se$V1
  output.df.se$upper = output.df$Coefficient + output.df.se$V1
  output.df.se$V1 = NULL
  
  #to long format for plotting
  output.df.long = gather(output.df, variable, value, 2:ncol(output.df))
  
  #adding the errors
  output.df.long$lower = c(output.df.se$lower, output.df.long[output.df.long$variable=="p-value", "value"], output.df.long[output.df.long$variable=="R2", "value"])
  output.df.long$upper = c(output.df.se$upper, output.df.long[output.df.long$variable=="p-value", "value"], output.df.long[output.df.long$variable=="R2", "value"])
  
  return(output.df.long)
}


#COMPUTING NULL
#randomizing lag column of both datasets 999 times and computing model with randomized data
modelRandomLagData <- function(lagged.data, model.formula, iterations){
  
  #getting response column
  response = all.vars(as.formula(model.formula))[1]
  
  #list to store results
  results = list()
  
  #computing one model per iteration
  for(i in 1:iterations){
    
    #randomization of the response column
    lagged.data[, response] = lagged.data[sample(1:nrow(lagged.data)), response]
    
    #computing model and storing results in list
    results[[i]] = modelLagData(model.formula=model.formula, lagged.data=lagged.data)
    
  }#end of iterations
  
  #preparing the data
  #getting lags and variable names
  null.model.df = data.frame(lag=results[[1]]$lag, variable=results[[1]]$variable, stringsAsFactors = FALSE)
  
  #getting the value column of each dataframe in each list
  null.model.values = sapply(results, `[[`, 3)
  
  null.model.df$value = apply(null.model.values,1, quantile, na.rm = TRUE, probs=0.5)
  null.model.df$upper = apply(null.model.values,1, quantile, na.rm = TRUE, probs=0.95)
  null.model.df$lower = apply(null.model.values,1, quantile, na.rm = TRUE, probs=0.05)
  
  return(null.model.df)
  
}
