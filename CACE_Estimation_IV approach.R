
#Practical Guide to Estimate Complier Average Causal Effect in the Context of Low Treatment Compliance: Utility of Instrumental Variable Approach

#Walk through of implementation using R code- January 2024

#Step 1
library(data.table)
library(simstudy)
library(ggplot2)
library(ivreg)# package to quickly estimate CACE in R using IV approach

#Step 2
# formulaa: Below we have created an object that saves our formula for distribution of always takers, never takers, compliers. For example, if we set proportion of compliers as 20%, then remaining 80% will be divided between never takers (40%) and always takers (40%). The distribution is mainly defined by proportion of compliers and this formula needs to be changed every time when we calculate CACE, ITT, and true causal effect at varying proportion of compliers. 
# Check: This object verifies at what proportion of compliers we need to calculate the parameters of interest (CACE, ITT,and true causal effect). Since, our proportion of compliers varies between 2 to 100%, this will range between 0.02 to 1.
#Cols: This object will define the names of columns for the generated data set. 
#Rows: Since we planned to run 20 simulations, the first 20 rows in the data set will save the data from 20 simulations followed by saving the data on descriptive statistics such as mean, minimum, maximum, median, lower and upper 95% CI values.
#Samples: define the number of simulations for the given distribution of always takers, never takers, compliers.
#Sim: This object stores the number of simulations to be used in "if" conditions later in the code (see below)

formulaa <- "0.005; 0.005; 0.99";
check <- 0.99;
cols <- c("Truth","ITT","CACE","Prop_C","IV");
rows <- c(1:20, "Mean","Min","Max","SD","Median", "LL", "UL");
Samples <- 1:20;
sim <- 20;

#Step 3
# assign: This function is used to automatically assign a name to a data set based on the variables created above. This will generate 5 cols and 27 rows as defined above; however, without names of the cols and rows. 
#data:This object is a dummy object, which will temporarily assign the same data 
#rlabels: A dummy object
#clabels:A dummy object
#rownames: This will assign names to all 27 rows 
#colnames: Will assign names to 5 columns
#assign function is used to assign names of columns and rows from data to generated dataset.
#rm: since we do not need rlables, clabels, and data any more, we remove these dummy objects
assign(paste0("data_", check, sep=""), matrix(0, nrow = sim+7, ncol = 5))

data <- get(paste0("data_", check, sep=""))
rlabels <- rows
clabels <- cols
rownames(data)<-rlabels
colnames(data)<-clabels
assign(paste0("data_", check, sep=""),data)
rm(rlabels,clabels,data)

#Step 4: Using “for loop” to create simulated data set for an RCT and perform 20 simulations on a data set of 10,000 observations to compute true causal effect, ITT, CACE, and IV estimand followed by generating their descriptive statistics at varying proportion of compliance starting at 2% and ending at 99%.

for (i in Samples){
  
  simulated <- defDataAdd(varname = "Subgroup", 
                          formula = formulaa, dist = "categorical")
  
  simulated <- defDataAdd(simulated, varname = "X0", 
                          formula = "(Subgroup == 1) * 1", dist = "nonrandom")
  simulated <- defDataAdd(simulated, varname = "X1", 
                          formula = "(Subgroup!= 2) * 1", dist = "nonrandom")
  
  simulated <- defDataAdd(simulated, varname = "x", 
                          formula = "(z==0) * X0 + (z==1) * X1", dist = "nonrandom")
  
  sumi <- genData(100000)
  sumi <- trtAssign(sumi, n=2, grpName = "z")
  sumi <- addColumns(simulated, sumi)
  sumi[, AStatus := factor(Subgroup, 
                           labels = c("Always-taker","Never-taker", "Complier"))]
  
  sumi[Subgroup == 1, Y0 := rnorm(.N, 1.0, sqrt(0.25))]
  sumi[Subgroup == 2, Y0 := rnorm(.N, 0.0, sqrt(0.36))]
  sumi[Subgroup == 3, Y0 := rnorm(.N, 0.1, sqrt(0.16))]
  
  sumi[Subgroup == 1, Y1 := rnorm(.N, 1.0, sqrt(0.25))]
  sumi[Subgroup == 2, Y1 := rnorm(.N, 0.0, sqrt(0.36))]
  sumi[Subgroup == 3, Y1 := rnorm(.N, 0.9, sqrt(0.49))]
  
  sumi[, y := (x == 0) * Y0 + (x == 1) * Y1]
  
# The above code chunk is used to create the variables of X1 and X0, which represent the value of received treatment after making random assignment of the treatment. We generate a data set, Sumi, with 10,000 observations and 7 variables (i.e., Id, z, subgroup, x1, x0, x, and AStatus). Id represents a unique identification number for a given row, Z represents the random assignment to the treatment or intervention in the context of individual RCT, X1 represents the potential outcome (status of treatment received) if treatment was assigned Z=1, and X0 represents the potential outcome (status of treatment received) if treatment was not assigned, Z=0. X represents the actual observed outcome under assigned treatment, Z=0 or Z=1. AStatus is the variable representing names of groups, compliers, always takers, and never takers. Finally, we added three additional variables i.e. Y1, Y0, and Y, which represent potential outcomes (Y1, Y0) and observed outcome (Y) after actually receiving the treatment (X=1 and X=0). This part of the code is modified using the original code from website: Goldfeld K. Complier average causal effect? Exploring what we learn from an RCT with participants who do not do what they are told. | R-bloggers. Accessed April 23, 2021. https://www.r-bloggers.com/2017/09/complier-average-causal-effect-exploring-what-we-learn-from-an-rct-with-participants-who-dont-do-what-they-are-told/. 
  
  
# “If condition” in the following code chunk helps to run 20 simulations in the data set (sumi) at proportion of compliers of 2%. 
  
  if(check == 0.02){
    
  data_0.02[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
  data_0.02[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
  data_0.02[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
  data_0.02[i,3] <- data_0.02[i,2]/data_0.02[i,4];
  IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
  data_0.02[i,5] <- IV_coef$coefficients[2];
  
# After running 20 simulations, the following code chunck helps to create mean, minimum, maximum, standard deviation, median, and lower and upper limits of 95% CIs of five parameters; Truth, ITT, CACE, Proportion of compliers, and IV estimand (CACE created using IVreg package in R). These parameters are generated at following distribution: Compliers (2%), Never takers (49%), and Always Takers (49%).
  
  if(i == sim){
    for(j in 1:5){
    data_0.02[21,j] <- mean(data_0.02[1:20,j]) 
    data_0.02[22,j] <- min(data_0.02[1:20,j])
    data_0.02[23,j] <- max(data_0.02[1:20,j])
    data_0.02[24,j] <- sd(data_0.02[1:20,j])
    data_0.02[25,j] <- median(data_0.02[1:20,j])
    CI_LUL <- t.test(data_0.02[1:20,j])
    data_0.02[26:27,j] <- CI_LUL$conf.int
    
    }
   rm(i,j,CI_LUL,IV_coef)
  }
  
  }
  
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 5%. 
  if(check == 0.05){
    
    data_0.05[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.05[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.05[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.05[i,3] <- data_0.05[i,2]/data_0.05[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.05[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.05[21,j] <- mean(data_0.05[1:20,j]) 
        data_0.05[22,j] <- min(data_0.05[1:20,j])
        data_0.05[23,j] <- max(data_0.05[1:20,j])
        data_0.05[24,j] <- sd(data_0.05[1:20,j])
        data_0.05[25,j] <- median(data_0.05[1:20,j])
        CI_LUL <- t.test(data_0.05[1:20,j])
        data_0.05[26:27,j] <- CI_LUL$conf.int
        
      }
    rm(i,j,CI_LUL,IV_coef)
   }
   
  }
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 10%. 
  if(check == 0.1){
    
    data_0.1[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.1[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.1[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.1[i,3] <- data_0.1[i,2]/data_0.1[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.1[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.1[21,j] <- mean(data_0.1[1:20,j]) 
        data_0.1[22,j] <- min(data_0.1[1:20,j])
        data_0.1[23,j] <- max(data_0.1[1:20,j])
        data_0.1[24,j] <- sd(data_0.1[1:20,j])
        data_0.1[25,j] <- median(data_0.1[1:20,j])
        CI_LUL <- t.test(data_0.1[1:20,j])
        data_0.1[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 15%. 
  if(check == 0.15){
    
    data_0.15[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.15[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.15[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.15[i,3] <- data_0.15[i,2]/data_0.15[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.15[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.15[21,j] <- mean(data_0.15[1:20,j]) 
        data_0.15[22,j] <- min(data_0.15[1:20,j])
        data_0.15[23,j] <- max(data_0.15[1:20,j])
        data_0.15[24,j] <- sd(data_0.15[1:20,j])
        data_0.15[25,j] <- median(data_0.15[1:20,j])
        CI_LUL <- t.test(data_0.15[1:20,j])
        data_0.15[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
  
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 20%. 
  if(check == 0.2){
    
    data_0.2[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.2[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.2[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.2[i,3] <- data_0.2[i,2]/data_0.2[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.2[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.2[21,j] <- mean(data_0.2[1:20,j]) 
        data_0.2[22,j] <- min(data_0.2[1:20,j])
        data_0.2[23,j] <- max(data_0.2[1:20,j])
        data_0.2[24,j] <- sd(data_0.2[1:20,j])
        data_0.2[25,j] <- median(data_0.2[1:20,j])
        CI_LUL <- t.test(data_0.2[1:20,j])
        data_0.2[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
  
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 25%.  
  if(check == 0.25){
    
    data_0.25[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.25[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.25[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.25[i,3] <- data_0.25[i,2]/data_0.25[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.25[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.25[21,j] <- mean(data_0.25[1:20,j]) 
        data_0.25[22,j] <- min(data_0.25[1:20,j])
        data_0.25[23,j] <- max(data_0.25[1:20,j])
        data_0.25[24,j] <- sd(data_0.25[1:20,j])
        data_0.25[25,j] <- median(data_0.25[1:20,j])
        CI_LUL <- t.test(data_0.25[1:20,j])
        data_0.25[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
  
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 30%.  
  if(check == 0.3){
    
    data_0.3[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.3[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.3[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.3[i,3] <- data_0.3[i,2]/data_0.3[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.3[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.3[21,j] <- mean(data_0.3[1:20,j]) 
        data_0.3[22,j] <- min(data_0.3[1:20,j])
        data_0.3[23,j] <- max(data_0.3[1:20,j])
        data_0.3[24,j] <- sd(data_0.3[1:20,j])
        data_0.3[25,j] <- median(data_0.3[1:20,j])
        CI_LUL <- t.test(data_0.3[1:20,j])
        data_0.3[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 35%. 
  if(check == 0.35){
    
    data_0.35[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.35[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.35[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.35[i,3] <- data_0.35[i,2]/data_0.35[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.35[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.35[21,j] <- mean(data_0.35[1:20,j]) 
        data_0.35[22,j] <- min(data_0.35[1:20,j])
        data_0.35[23,j] <- max(data_0.35[1:20,j])
        data_0.35[24,j] <- sd(data_0.35[1:20,j])
        data_0.35[25,j] <- median(data_0.35[1:20,j])
        CI_LUL <- t.test(data_0.35[1:20,j])
        data_0.35[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
  
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 40%.  
  if(check == 0.4){
    
    data_0.4[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.4[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.4[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.4[i,3] <- data_0.4[i,2]/data_0.4[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.4[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.4[21,j] <- mean(data_0.4[1:20,j]) 
        data_0.4[22,j] <- min(data_0.4[1:20,j])
        data_0.4[23,j] <- max(data_0.4[1:20,j])
        data_0.4[24,j] <- sd(data_0.4[1:20,j])
        data_0.4[25,j] <- median(data_0.4[1:20,j])
        CI_LUL <- t.test(data_0.4[1:20,j])
        data_0.4[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 45%.   
  if(check == 0.45){
    
    data_0.45[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.45[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.45[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.45[i,3] <- data_0.45[i,2]/data_0.45[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.45[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.45[21,j] <- mean(data_0.45[1:20,j]) 
        data_0.45[22,j] <- min(data_0.45[1:20,j])
        data_0.45[23,j] <- max(data_0.45[1:20,j])
        data_0.45[24,j] <- sd(data_0.45[1:20,j])
        data_0.45[25,j] <- median(data_0.45[1:20,j])
        CI_LUL <- t.test(data_0.45[1:20,j])
        data_0.45[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
  
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 50%.   
  if(check == 0.5){
    
    data_0.5[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.5[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.5[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.5[i,3] <- data_0.5[i,2]/data_0.5[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.5[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.5[21,j] <- mean(data_0.5[1:20,j]) 
        data_0.5[22,j] <- min(data_0.5[1:20,j])
        data_0.5[23,j] <- max(data_0.5[1:20,j])
        data_0.5[24,j] <- sd(data_0.5[1:20,j])
        data_0.5[25,j] <- median(data_0.5[1:20,j])
        CI_LUL <- t.test(data_0.5[1:20,j])
        data_0.5[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
  
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 55%.  
  if(check == 0.55){
    
    data_0.55[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.55[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.55[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.55[i,3] <- data_0.55[i,2]/data_0.55[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.55[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.55[21,j] <- mean(data_0.55[1:20,j]) 
        data_0.55[22,j] <- min(data_0.55[1:20,j])
        data_0.55[23,j] <- max(data_0.55[1:20,j])
        data_0.55[24,j] <- sd(data_0.55[1:20,j])
        data_0.55[25,j] <- median(data_0.55[1:20,j])
        CI_LUL <- t.test(data_0.55[1:20,j])
        data_0.55[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
  
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 60%.    
  if(check == 0.6){
    
    data_0.6[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.6[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.6[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.6[i,3] <- data_0.6[i,2]/data_0.6[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.6[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.6[21,j] <- mean(data_0.6[1:20,j]) 
        data_0.6[22,j] <- min(data_0.6[1:20,j])
        data_0.6[23,j] <- max(data_0.6[1:20,j])
        data_0.6[24,j] <- sd(data_0.6[1:20,j])
        data_0.6[25,j] <- median(data_0.6[1:20,j])
        CI_LUL <- t.test(data_0.6[1:20,j])
        data_0.6[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
  
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 65%.  
  if(check == 0.65){
    
    data_0.65[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.65[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.65[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.65[i,3] <- data_0.65[i,2]/data_0.65[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.65[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.65[21,j] <- mean(data_0.65[1:20,j]) 
        data_0.65[22,j] <- min(data_0.65[1:20,j])
        data_0.65[23,j] <- max(data_0.65[1:20,j])
        data_0.65[24,j] <- sd(data_0.65[1:20,j])
        data_0.65[25,j] <- median(data_0.65[1:20,j])
        CI_LUL <- t.test(data_0.65[1:20,j])
        data_0.65[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 70%.   
  if(check == 0.7){
    
    data_0.7[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.7[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.7[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.7[i,3] <- data_0.7[i,2]/data_0.7[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.7[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.7[21,j] <- mean(data_0.7[1:20,j]) 
        data_0.7[22,j] <- min(data_0.7[1:20,j])
        data_0.7[23,j] <- max(data_0.7[1:20,j])
        data_0.7[24,j] <- sd(data_0.7[1:20,j])
        data_0.7[25,j] <- median(data_0.7[1:20,j])
        CI_LUL <- t.test(data_0.7[1:20,j])
        data_0.7[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
  
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 75%.   
  if(check == 0.75){
    
    data_0.75[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.75[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.75[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.75[i,3] <- data_0.75[i,2]/data_0.75[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.75[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.75[21,j] <- mean(data_0.75[1:20,j]) 
        data_0.75[22,j] <- min(data_0.75[1:20,j])
        data_0.75[23,j] <- max(data_0.75[1:20,j])
        data_0.75[24,j] <- sd(data_0.75[1:20,j])
        data_0.75[25,j] <- median(data_0.75[1:20,j])
        CI_LUL <- t.test(data_0.75[1:20,j])
        data_0.75[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 80%.   
  if(check == 0.8){
    
    data_0.8[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.8[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.8[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.8[i,3] <- data_0.8[i,2]/data_0.8[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.8[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.8[21,j] <- mean(data_0.8[1:20,j]) 
        data_0.8[22,j] <- min(data_0.8[1:20,j])
        data_0.8[23,j] <- max(data_0.8[1:20,j])
        data_0.8[24,j] <- sd(data_0.8[1:20,j])
        data_0.8[25,j] <- median(data_0.8[1:20,j])
        CI_LUL <- t.test(data_0.8[1:20,j])
        data_0.8[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 85%.   
  if(check == 0.85){
    
    data_0.85[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.85[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.85[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.85[i,3] <- data_0.85[i,2]/data_0.85[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.85[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.85[21,j] <- mean(data_0.85[1:20,j]) 
        data_0.85[22,j] <- min(data_0.85[1:20,j])
        data_0.85[23,j] <- max(data_0.85[1:20,j])
        data_0.85[24,j] <- sd(data_0.85[1:20,j])
        data_0.85[25,j] <- median(data_0.85[1:20,j])
        CI_LUL <- t.test(data_0.85[1:20,j])
        data_0.85[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 90%.  
  if(check == 0.9){
    
    data_0.9[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.9[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.9[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.9[i,3] <- data_0.9[i,2]/data_0.9[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.9[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.9[21,j] <- mean(data_0.9[1:20,j]) 
        data_0.9[22,j] <- min(data_0.9[1:20,j])
        data_0.9[23,j] <- max(data_0.9[1:20,j])
        data_0.9[24,j] <- sd(data_0.9[1:20,j])
        data_0.9[25,j] <- median(data_0.9[1:20,j])
        CI_LUL <- t.test(data_0.9[1:20,j])
        data_0.9[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 95%.   
  if(check == 0.95){
    
    data_0.95[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.95[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.95[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.95[i,3] <- data_0.95[i,2]/data_0.95[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.95[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.95[21,j] <- mean(data_0.95[1:20,j]) 
        data_0.95[22,j] <- min(data_0.95[1:20,j])
        data_0.95[23,j] <- max(data_0.95[1:20,j])
        data_0.95[24,j] <- sd(data_0.95[1:20,j])
        data_0.95[25,j] <- median(data_0.95[1:20,j])
        CI_LUL <- t.test(data_0.95[1:20,j])
        data_0.95[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
# If condition in the following code chunk helps to run 20 simulations in the data set sumi at proportion of compliers of 99%.     
  if(check == 0.99){
    
    data_0.99[i,1] <- sumi [AStatus == "Complier", mean (Y1 - Y0)];
    data_0.99[i,2] <- sumi[z==1, mean(y)] - sumi[z==0, mean(y)];
    data_0.99[i,4] <- sumi [z==1, mean(x)] - sumi [z==0, mean(x)];
    data_0.99[i,3] <- data_0.99[i,2]/data_0.99[i,4];
    IV_coef <- ivreg(formula = y ~ x | z, data = sumi, x = TRUE)
    data_0.99[i,5] <- IV_coef$coefficients[2];
    
    if(i == sim){
      for(j in 1:5){
        data_0.99[21,j] <- mean(data_0.99[1:20,j]) 
        data_0.99[22,j] <- min(data_0.99[1:20,j])
        data_0.99[23,j] <- max(data_0.99[1:20,j])
        data_0.99[24,j] <- sd(data_0.99[1:20,j])
        data_0.99[25,j] <- median(data_0.99[1:20,j])
        CI_LUL <- t.test(data_0.99[1:20,j])
        data_0.99[26:27,j] <- CI_LUL$conf.int
        
      }
      rm(i,j,CI_LUL,IV_coef)
    }
    
  }
  
};


#Step 5: Using mean values of true causal effect, ITT, CACE, and proportion of compliers from 21 simulated data sets (compliance between 2 to 99%) to generate a data set for creating a line graph.


df_line <- matrix(0, nrow = 21, ncol = 4);
rownames(df_line) <- c(1:21);
colnames(df_line) <- c("Prop-C", "Truth", "ITT", "CACE");
df_line[1:21,1] <- as.numeric(c(0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99));
df_data <- data.frame(df_line)


  df_data[1,2] <- data_0.02[21,1];
  df_data[1,3] <- data_0.02[21,2];
  df_data[1,4] <- data_0.02[21,3];
  df_data[2,2] <- data_0.05[21,1];
  df_data[2,3] <- data_0.05[21,2];
  df_data[2,4] <- data_0.05[21,3];
  df_data[3,2] <- data_0.1[21,1];
  df_data[3,3] <- data_0.1[21,2];
  df_data[3,4] <- data_0.1[21,3];
  df_data[4,2] <- data_0.15[21,1];
  df_data[4,3] <- data_0.15[21,2];
  df_data[4,4] <- data_0.15[21,3];
  df_data[5,2] <- data_0.2[21,1];
  df_data[5,3] <- data_0.2[21,2];
  df_data[5,4] <- data_0.2[21,3];
  df_data[6,2] <- data_0.25[21,1];
  df_data[6,3] <- data_0.25[21,2];
  df_data[6,4] <- data_0.25[21,3];
  df_data[7,2] <- data_0.3[21,1];
  df_data[7,3] <- data_0.3[21,2];
  df_data[7,4] <- data_0.3[21,3];
  df_data[8,2] <- data_0.35[21,1];
  df_data[8,3] <- data_0.35[21,2];
  df_data[8,4] <- data_0.35[21,3];
  df_data[9,2] <- data_0.4[21,1];
  df_data[9,3] <- data_0.4[21,2];
  df_data[9,4] <- data_0.4[21,3];
  df_data[10,2] <- data_0.45[21,1];
  df_data[10,3] <- data_0.45[21,2];
  df_data[10,4] <- data_0.45[21,3];
  df_data[11,2] <- data_0.5[21,1];
  df_data[11,3] <- data_0.5[21,2];
  df_data[11,4] <- data_0.5[21,3];
  df_data[12,2] <- data_0.55[21,1];
  df_data[12,3] <- data_0.55[21,2];
  df_data[12,4] <- data_0.55[21,3];
  df_data[13,2] <- data_0.6[21,1];
  df_data[13,3] <- data_0.6[21,2];
  df_data[13,4] <- data_0.6[21,3];
  df_data[14,2] <- data_0.65[21,1];
  df_data[14,3] <- data_0.65[21,2];
  df_data[14,4] <- data_0.65[21,3];
  df_data[15,2] <- data_0.7[21,1];
  df_data[15,3] <- data_0.7[21,2];
  df_data[15,4] <- data_0.7[21,3];
  df_data[16,2] <- data_0.75[21,1];
  df_data[16,3] <- data_0.75[21,2];
  df_data[16,4] <- data_0.75[21,3];
  df_data[17,2] <- data_0.8[21,1];
  df_data[17,3] <- data_0.8[21,2];
  df_data[17,4] <- data_0.8[21,3];
  df_data[18,2] <- data_0.85[21,1];
  df_data[18,3] <- data_0.85[21,2];
  df_data[18,4] <- data_0.85[21,3];
  df_data[19,2] <- data_0.9[21,1];
  df_data[19,3] <- data_0.9[21,2];
  df_data[19,4] <- data_0.9[21,3];
  df_data[20,2] <- data_0.95[21,1];
  df_data[20,3] <- data_0.95[21,2];
  df_data[20,4] <- data_0.95[21,3];
  df_data[21,2] <- data_0.99[21,1];
  df_data[21,3] <- data_0.99[21,2];
  df_data[21,4] <- data_0.99[21,3];
  
## Line chart can be generated in R using df_data set and ggplot.  


