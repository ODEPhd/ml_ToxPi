#Clear the global environment
rm(list=ls())

# RF regression
# Load packages

library(tidyverse)
library(caret)
library(Boruta)

# Load the data
ToxPi_variables <- read.csv(file.choose(), header = T)
row.names(ToxPi_variables) <- ToxPi_variables$X
ToxPi_variables <- ToxPi_variables[,-1]

# partition data into training and test sets
set.seed(111)
inTrain <- createDataPartition(ToxPi_variables$ToxPi_scores, p = 0.7, list = FALSE)
training <- ToxPi_variables[inTrain,]
test <- ToxPi_variables[-inTrain,]


##################################################                       
#Feature Selection by boruta algorithm
##################################################

#Load the required package
set.seed(111)
boruta.initial <- Boruta(ToxPi_scores ~ ., data = training, doTrace = 1, maxRuns = 500,
                         pValue=0.01, mcAdj = TRUE)
print(boruta.initial)
plOt1 <- plot(boruta.initial, las = 2, cex.axis = 0.8)
# plot2  <- plotImpHistory(boruta.initial)

### Save boruta plot1 as a high quality png file ###
png(filename = "Boruta importance plot_originally confirmed.png", width = 10, height = 8, units = "in", res = 600)
plot(boruta.initial, whichShadow = c(TRUE, TRUE, TRUE),
     las = 2, cex.axis = 0.5, boxwex = 0.6, xlab = "", main = "")
dev.off()



#############################################
# randomForest (package/method) in caret
##############################################

# fit regression model
set.seed(111)
system.time(rf_fit <- train(ToxPi_scores ~ ., data = training, method = 'rf', metric = "RMSE",
                           importance=TRUE, ntrees=500, preProcess = c("center", "scale"),
                            trainControl = trainControl(method = 'repeatedcv',
                                                        number = 5, repeats = 10)))
print(rf_fit)


# Tune the rf model hyperparameters
set.seed(111)
grid_rf <- expand.grid(.mtry = seq(1, 5, by=1))
system.time(rf_fit <- train(ToxPi_scores ~ ., data = training, method = 'rf', metric = "RMSE",
                importance = TRUE, ntrees=500,
                preProcess = c("center", "scale"), tuneGrid = grid_rf,
                trainControl = trainControl(method = 'repeatedcv',
                                            number = 5, repeats = 10)))

print(rf_fit) #print the training data rf model prediction results 
plot(rf_fit) #plot the training data rf model prediction results 

### Save rf_fit plot as a high quality png file ###
png(filename = "rf_fit plot.png", width = 10, height = 8, units = "in", res = 600)
plot(rf_fit)
dev.off()

varImp(rf_fit) #print the variable importance scores
plot(varImp(rf_fit)) #plot the variable importance scores

### Save varImp plot as a high quality png file ###
png(filename = "varImp plot.png", width = 10, height = 8, units = "in", res = 600)
plot(varImp(rf_fit))
dev.off()


# Furthermore, we can find the standard deviation around the
# Rsquared value by examining the R-squared from each fold.
sd(rf_fit$resample$Rsquared)


# model performance on the test data
testPred_rf_fit <- predict(rf_fit, newdata=test, mtry=3)
rmse_rf_fit<- RMSE(testPred_rf_fit,test$ToxPi_scores)
r2_rf_fit <- R2(testPred_rf_fit,test$ToxPi_scores)

print(rmse_rf_fit)
print(r2_rf_fit)



##################################################################
#Create partial dependence plots (pdp) for all variables_with rugs
##################################################################

library(ggplot2)
library(pdp)
library(gridExtra)

# Function to create and save partial dependence plots with rug
save_partial_dependence_plots <- function(rf_fit, filename) {
  
  # List of variables to plot
  variables <- c("DF", "Water_solubility", "Peak_areas", "Half_lives", "LogKow", 
                 "Exposure_category", "Koc", "Bioactivity_ratio", "IARC_category")
  
  # Titles for the x-axes with units
  x_labels <- c("Detection frequency (%)", 
                "Water solubility (µg/L)", 
                "Peak areas (arb.unit)", 
                "Half-lives (days)", 
                "LogKow (unitless)", 
                "Exposure category (unitless)", 
                "Koc (L/kg)", 
                "Bioactivity ratio (unitless)", 
                "IARC category (unitless)")
  
  # Response variable label with unit/context
  y_label <- "Predicted ToxPi Score (unitless)"
  
    # List to store plots
  plots <- list()
  
  # Loop through each variable
  for (i in seq_along(variables)) {
    var <- variables[i]
    x_label <- x_labels[i]
    
    # Get partial dependence data
    pd_data <- pdp::partial(rf_fit, pred.var = var, plot = FALSE)
    
    # Create the plot
    p <- ggplot(pd_data, aes_string(x = var, y = "yhat")) +
      geom_smooth(method = "auto", se = FALSE, fullrange = FALSE, linewidth = 3) +
      geom_rug(sides = "b", color = "black") +  # Add rug plot
      theme_bw() +
      theme(legend.position = "none") +
      xlab(x_label) +
      theme(axis.title.x = element_text(family = "arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue")) +
      theme(axis.title.y = element_text(family = "arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue")) +
      theme(axis.text.x = element_text(size = 12, color = "Black", face = "bold", angle = 0, vjust = 0.5, hjust = 0.5)) +
      theme(axis.text.y = element_text(size = 12, color = "Black", face = "bold", angle = 0, vjust = 0.5, hjust = 1))
    
    # Conditionally add y-axis label only to the fourth plot (index 4)
    if (i == 1 || i == 4 || i == 7) {
      p <- p + ylab(y_label)
    } else {
      p <- p + theme(axis.title.y = element_blank())
    }
    
    # Add the plot to the list
    plots[[i]] <- p
  }
  
  # Arrange all plots in a grid
  grid_plots <- do.call(grid.arrange, c(plots, ncol = 3))
  
  # Save the grid of plots to a file
  ggsave(filename = filename, plot = grid_plots, width = 12, height = 11)
}

# Assuming rf_fit is your pre-trained random forest model
# Call the function to save all partial dependence plots as a grid
save_partial_dependence_plots(rf_fit, "partial_dependence_plots_with_rug.png")



# ######################################################################
# #Create partial dependence plots (pdp) for all variables_without rugs
#Uncomment the relevant codes if you need to use them 
# #####################################################################
# 
# library(ggplot2)
# library(pdp)
# library(gridExtra)
# 
# # Function to create and save partial dependence plots
# save_partial_dependence_plots <- function(rf_fit, filename) {
#   
#   # List of variables to plot
#   variables <- c("DF", "Water_solubility", "Peak_areas", "Half_lives", "LogKow", 
#                  "Exposure_category", "Koc", "Bioactivity_ratio", "IARC_category")
#   
#   # Titles for the x-axes with units
#   x_labels <- c("Detection frequency (%)", 
#                 "Water solubility (µg/L)", 
#                 "Peak areas (arb.unit)", 
#                 "Half-lives (days)", 
#                 "LogKow (unitless)", 
#                 "Exposure category (unitless)", 
#                 "Koc (L/kg)", 
#                 "Bioactivity ratio (unitless)", 
#                 "IARC category (unitless)")
#   
#   # Response variable label with unit/context
#   y_label <- "Predicted ToxPi Score (unitless)"
#   
#   # List to store plots
#   plots <- list()
#   
#   # Loop through each variable
#   for (i in seq_along(variables)) {
#     var <- variables[i]
#     x_label <- x_labels[i]
#     
#     # Get partial dependence data
#     pd_data <- pdp::partial(rf_fit, pred.var = var, plot = FALSE)
#     
#     # Create the plot
#     p <- ggplot(pd_data, aes_string(x = var, y = "yhat")) +
#       geom_smooth(method = "auto", se = FALSE, fullrange = FALSE, linewidth = 3) +
#       theme_bw() +
#       theme(legend.position = "none") +
#       xlab(x_label) +
#       theme(axis.title.x = element_text(family = "arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue")) +
#       theme(axis.title.y = element_text(family = "arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 14, color = "darkblue")) +
#       theme(axis.text.x = element_text(size = 12, color = "Black", face = "bold", angle = 0, vjust = 0.5, hjust = 0.5)) +
#       theme(axis.text.y = element_text(size = 12, color = "Black", face = "bold", angle = 0, vjust = 0.5, hjust = 1))
#     
#     
#     # Conditionally add y-axis label only to the fourth plot (index 4)
#     if (i == 1 || i == 4 || i == 7) {
#       p <- p + ylab(y_label)
#     } else {
#       p <- p + theme(axis.title.y = element_blank())
#     }
#     
#     
#     # Add the plot to the list
#     plots[[i]] <- p
#   }
#   
#   # Arrange all plots in a grid
#   grid_plots <- do.call(grid.arrange, c(plots, ncol = 3))
#   
#   # Save the grid of plots to a file
#   ggsave(filename = filename, plot = grid_plots, width = 12, height = 11, dpi = 600)
# }
# 
# # Considering that rf_fit is our pre-trained random forest model
# # We call the function to save all partial dependence plots as a grid
# save_partial_dependence_plots(rf_fit, "partial_dependence_plots.png")
# 


####################################################
# support vector machines (package/method) in caret
###################################################

# fit regression model
set.seed(111)
system.time(svm_fit <- train(ToxPi_scores ~ ., data = training, method = 'svmRadial', metric = "RMSE",
                            preProcess = c("center", "scale"),
                            trainControl = trainControl(method = 'repeatedcv',
                                                        number = 5, repeats = 10)))
print(svm_fit)


# tune the svm model hyperparameters
set.seed(111)
grid_svm <- expand.grid(.sigma = c(0.1, 0.15, 0.2, 0.25, 0.3), .C=seq(1, 10, by=1))
system.time(svm_fit <- train(ToxPi_scores ~ ., data = training, method = 'svmRadial', metric = "RMSE",
                            importance = TRUE, preProcess = c("center", "scale"), tuneGrid = grid_svm,
                            trainControl = trainControl(method = 'repeatedcv',
                                                        number = 5, repeats = 10)))

print(svm_fit) #print the training data svm model prediction results 
plot(svm_fit) #plot the training data svm model prediction results 

### Save svm_fit plot as a high quality png file ###
png(filename = "svm_fit plot.png", width = 10, height = 8, units = "in", res = 600)
plot(svm_fit)
dev.off()

# model performance 
testPred_svm_fit <- predict(svm_fit, newdata=test, Sigma=0.1, C=3)
rmse_svm_fit<- RMSE(testPred_svm_fit,test$ToxPi_scores)
r2_svm_fit <- R2(testPred_svm_fit,test$ToxPi_scores)

print(rmse_svm_fit)
print(r2_svm_fit)

plot(testPred_svm_fit, test$ToxPi_scores,
     main = "Actual vs Predicted",
     xlab = "Test Set Value", ylab = "Predicted Value")

abline(lm(test$ToxPi_scores ~ testPred_svm_fit), col="black")


#############################################
# Neural networks (package/method) in caret
##############################################

# fit regression model
set.seed(111)
system.time(nnet_fit <- train(ToxPi_scores ~ ., data = training, method = 'nnet',
                            metric="RMSE", preProcess = c("center", "scale"),
                             ## Reduce the amount of printed output
                             trace=FALSE,
                             ## Expand the number of iterations to find
                             ## parameter estimates
                             maxit=500, linout=TRUE,
                             trainControl = trainControl(method = 'repeatedcv',
                                                         number = 5, repeats = 10)))
print(nnet_fit)


# tune the nnet model hyperparameters
set.seed(111)
grid_nnet <- expand.grid(decay = c(0.01, 0.015, 0.1, 0.15, 0.2, 0.25, 0.3),
                         size=seq(1, 10, by=1))

system.time(nnet_fit <- train(ToxPi_scores ~ ., data = training, method = 'nnet',
                              metric="RMSE", preProcess = c("center", "scale"),
                              ## Reduce the amount of printed output
                              trace=FALSE,
                              ## Expand the number of iterations to find
                              ## parameter estimates
                              maxit=500, linout=TRUE, tuneGrid = grid_nnet,
                              trainControl = trainControl(method = 'repeatedcv',
                                                          number = 5, repeats = 10)))

print(nnet_fit) #print the training data svm model prediction results 
plot(nnet_fit) #plot the training data svm model prediction results 

### Save nnet_fit plot as a high quality png file ###
png(filename = "nnet_fit plot.png", width = 10, height = 8, units = "in", res = 600)
plot(nnet_fit)
dev.off()

# model performance 
testPred_nnet_fit <- predict(nnet_fit, newdata=test, size=2, decay=0.01)
rmse_nnet_fit<- RMSE(testPred_nnet_fit,test$ToxPi_scores)
r2_nnet_fit <- R2(testPred_nnet_fit,test$ToxPi_scores)

print(rmse_nnet_fit)
print(r2_nnet_fit)

plot(testPred_nnet_fit, test$ToxPi_scores,
     main = "Actual vs Predicted",
     xlab = "Test Set Value", ylab = "Predicted Value")

abline(lm(test$ToxPi_scores ~ testPred_nnet_fit), col="black")

