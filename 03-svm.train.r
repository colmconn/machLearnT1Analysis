if ( 1 == 1 ) { 
model.formula=as.formula(paste("Grp", "~", paste(colnames(feature.df)[grep(roi.label.regexp, colnames(feature.df))], collapse=" + ")))

## load caret
library(caret)

## enable parallel processing
library(doMC)
registerDoMC(cores = 40)

set.seed(12345)
cat("*** Creating training and testing partitions\n")
in.training <- createDataPartition(feature.df$Grp, p = .80, list = FALSE)
training <- feature.df[ in.training,]
testing  <- feature.df[-in.training,]

cat("*** Creating ML fit control list\n")
svm.control =
    trainControl(
        method="LOOCV",
        verboseIter=FALSE,
        classProbs=TRUE)

svm.grid = expand.grid(
    C=seq.int(100, 1000, 10))

train.svm.model <- function (in.model.formula, in.training, in.train.control, in.grid) {
    set.seed(12345)
    cat("*** Starting to train the model at:", date(), "\n")
    start=Sys.time()
    
    fit <- train(in.model.formula, data=in.training,
                                method="svmRadialCost",
                                trControl = in.train.control,
                                tuneGrid=in.grid)
    end=Sys.time()
    cat("*** Finished training the model at:", date(), "\n")
    suppressMessages(library(chron))
    cat("*** Computation took", format(as.chron(end) - as.chron(start)), "\n")

    return(fit)
}

## no sampling to deal with class imbalance
cat("*** No sampling to deal with class imbalance\n")
no.sampling.svm.fit = train.svm.model(model.formula, training, svm.control, svm.grid)

## down sampling to deal with class imbalance
cat("*** Down sampling to deal with class imbalance\n")
svm.control$sampling="down"
down.sampling.svm.fit = train.svm.model(model.formula, training, svm.control, svm.grid)

## up sampling to deal with class imbalance
cat("*** Up sampling to deal with class imbalance\n")
svm.control$sampling="up"
up.sampling.svm.fit = train.svm.model(model.formula, training, svm.control, svm.grid)


## ROSE sampling to deal with class imbalance
cat("*** ROSE sampling to deal with class imbalance\n")
svm.control$sampling="rose"
rose.sampling.svm.fit = train.svm.model(model.formula, training, svm.control, svm.grid)

## smote sampling to deal with class imbalance
cat("*** SMOTE sampling to deal with class imbalance\n")
svm.control$sampling="smote"
smote.sampling.svm.fit = train.svm.model(model.formula, training, svm.control, svm.grid)




## cat("*** Testing predictions\n")
## test.pred = predict(svm.fit, testing)
## cm=confusionMatrix  (test.pred, testing[complete.cases(testing), "Grp"])
## print(cm)


}

## fit.control <-
##     trainControl(
##         ## 10-fold CV
##         method = "repeatedcv",
##         number = 10,
##         ## repeated ten times
##         repeats = 10)

## gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9),
##                         n.trees = (1:100)*50,
##                         shrinkage = 0.1,
##                         n.minobsinnode = 20)
## gbmFit1 <- train(model.formula, data = training,
##                  method = "gbm",
##                  trControl = fit.control,
##                  ## This last option is actually one
##                  ## for gbm() that passes through
##                  verbose = FALSE,
##                  tuneGrid = gbmGrid)
