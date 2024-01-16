rm(list = ls())
library(randomForest)

##################FUNCTION DEFINITION##################
# Function to perform one-hot encoding: converts from categorial to numeric data 
one_hot_encode <- function(categories) {
  unique_categories <- unique(categories)
  num_categories <- length(unique_categories)
  
  encoded_matrix <- matrix(0, nrow = length(categories), ncol = num_categories,
                           dimnames = list(NULL, unique_categories))
  
  for (i in 1:length(categories)) {
    category_index <- match(categories[i], unique_categories)
    encoded_matrix[i, category_index] <- 1
  }
  
  return(encoded_matrix)
}

iucn_classification<-function(input.file){
    sep=":"
    thr           <-0.5
    cat("Reading the input data\n")
    #read the input data
    input.data <- read.csv(input.file)
    columnsToExclude<-c("SpecCode","Class","Family","Genus","Species","IUCN_Code","AquariumShorta","ImportanceShort","AnaCatShort","IUCN_Short")
    Selected.data.completecolumns<-input.data
    Selected.data.completecolumns<-Selected.data.completecolumns[,!( names(Selected.data.completecolumns) %in% columnsToExclude ) ]
    
    completeCases<-which(complete.cases(Selected.data.completecolumns))
    Original.data.completecolumns<-input.data[completeCases,]
    Selected.data.completecolumns <- Selected.data.completecolumns[completeCases, ]
    
    nonemptycases<-apply(Selected.data.completecolumns, 1, function(row) !any(row == ""))
    Selected.data.completecolumns <- Selected.data.completecolumns[nonemptycases, ]
    
    Original.data.completecolumns<-Original.data.completecolumns[nonemptycases,]
    test.data<-data.frame(matrix(data = NA, nrow = dim(Selected.data.completecolumns) [1], ncol = 0))
    test.data.colnames<-c()
  
    #collect the different data types
    categorial.columns <- c()
    numeric.columns <- c()
    boolean.columns <-c()
    cat("loading required scaling file data.scaling.bin\n")
    load("data.scaling.bin")
    
    for (c in names(Selected.data.completecolumns)){
      if (is.numeric(Selected.data.completecolumns[,c]))
        numeric.columns<-c(numeric.columns,c)
      else if (is.logical(Selected.data.completecolumns[,c]))
        boolean.columns<-c(boolean.columns,c)
      else
        categorial.columns<-c(categorial.columns,c)
    }
    
    for (c in numeric.columns){
      #standardise the numeric column and add it to the data frame
      column<-Selected.data.completecolumns[,c]
      ref_col<-data_for_scaling[,c]
      column<-(column-mean(ref_col))/sd(ref_col)
      
      test.data<-cbind(test.data, column)
      test.data.colnames<-c(test.data.colnames,c)
    }
    
    for (c in categorial.columns){
      #add one-hot encoded categorial features
      one_hot_encoded <- one_hot_encode(Selected.data.completecolumns[,c])
      test.data<-cbind(test.data, one_hot_encoded[,,drop=F])
      test.data.colnames<-c(test.data.colnames, paste0(c,sep,colnames(one_hot_encoded)) )
    }
    
    for (c in boolean.columns){
      #convert T/F to 1/0
      test.data<-cbind(test.data, as.numeric(Selected.data.completecolumns[,c]))
      test.data.colnames<-c(test.data.colnames,c)
    }
    
    names(test.data)<-test.data.colnames
    cat("loading required model file iucn_rf_model.bin\n")
    load(file = "iucn_rf_model.bin")
    test.1<-test.data[,test.data.colnames]
    testcolumnset<-paste0("i",seq(1,length(test.data.colnames)))
    names(test.1)<-testcolumnset
    options(warn = -1)
    set.seed(20)
    pr.nn1 <- predict(rf_model, newdata = test.1)
    pr.nn.t<-pr.nn1
    pr.nn.t[which(pr.nn1>=thr)]<-1
    pr.nn.t[which(pr.nn1<thr)]<-0
    pr.nn.t<-as.numeric(pr.nn.t)
    Original.data.completecolumns$threatened<-pr.nn.t
    
    outputfile <- gsub(pattern=".csv", x=input.file, replacement = "_iucn_ml_assessed.csv")
    
    cat("Saving the result to",outputfile,"including the new column \"threatened\" containing a binary IUCN status (0=Least Concern; 1=Threatened)\n")
    write.csv(file = outputfile, x = Original.data.completecolumns, row.names = F)

}
##################END OF FUNCTION DEFINITION##################

input.file<-"AAA_RF_IUCN_240115_test.csv"
iucn_classification(input.file = input.file)