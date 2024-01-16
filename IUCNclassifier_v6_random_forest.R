rm(list = ls())
library(randomForest)
library(ClusterR)

#####FEATURE SELECTION parameters
feature_selection_pca<-F #activate (T) / deactivate (F) feature selection
#select the Principal Components (PCs) explaining up to 95% of the variance
eigenvector_threshold<-0.95
#select the features contributing up to 90% to the PCs
feature_threshold<-0.9
#enable characterisation analysis on the most important features
importance.cumulative.threshold<-0.8
do.characterisation.analysis<-F

#####NUMERIC FEATURE TRANFORMATION
transform.numeric.columns.to.categorial<-F

#####ERROR ANALYSIS: extracts confusing species with different classifications but overlapping feature values
do.error.analysis<-F

#####RF parameters
crossvalidate <- T  # Do cross-validation (T/F)
k             <- 20 # k-fold cross validation
thr           <-0.5 # dichotomic decision threshold on the RF output
n.trees       <-500 #number of decision trees from which an ensemble assessment should be extracted
CV.proportion <- (1-(1/k)) # Set to 0.995 for LOOCV

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

cat("Reading the input data\n")
#read the input data
Reference.data <- read.csv("AAA_RF_IUCN_240115.csv")

cat("Filtering the input data\n")
#NOTE: adjust this part to include/exclude more features
#exclude the records with at least an NA and an empty string
Selected.data.completecolumns<-Reference.data
#test for excluding one feature at a time

#columnsToExclude<-c("SpecCode","Class","Order","Family","Genus","Species","IUCN_Code","AquariumShorta","Importance")
columnsToExclude<-c("SpecCode","Class","Family","Genus","Species","IUCN_Code","AquariumShorta","ImportanceShort","AnaCatShort")
#all-fulfilled features:
#columnsToExclude<-c("SpecCode","Class","Family","Genus","Order", "Species", "IUCN_Code","AquariumShorta","ImportanceShort","PD50","AnaCatShort","MaxLength","Troph","AquariumShort","GameFish","Importance","RepGuild1","CountOfC_Code","CountOfAreaCode")

if(!do.error.analysis){ #exclude the taxonomic ranks only if we are not conducting an error analysis of the overlaps between differently classified data
  Selected.data.completecolumns<-Selected.data.completecolumns[,!( names(Selected.data.completecolumns) %in% columnsToExclude ) ]
}
Selected.data.completecolumns <- Selected.data.completecolumns[complete.cases(Selected.data.completecolumns), ]
Selected.data.completecolumns <- Selected.data.completecolumns[apply(Selected.data.completecolumns, 1, function(row) !any(row == "")), ]
Selected.data.completecolumns.withSpp<-Selected.data.completecolumns
#the following line is for subsetting only - comment it for regime usage
#columnsToExclude<-c("GameFish","Brack","AquariumShort")
#columnsToExclude<-c(columnsToExclude,"Brack")
#columnsToExclude<-c("SpecCode","Class","Family","Genus","Order", "Species", "IUCN_Code","AquariumShorta","ImportanceShort","PD50","AnaCatShort","MaxLength","Troph","AquariumShort","GameFish","Importance","RepGuild1","CountOfC_Code","CountOfAreaCode")
Selected.data.completecolumns<-Selected.data.completecolumns[,!( names(Selected.data.completecolumns) %in% columnsToExclude ) ]

data_for_scaling<-Selected.data.completecolumns
save(file="data.scaling.bin",data_for_scaling)
cat("NOTE - TRAINING DATA FOR SCALING SAVED TO =","data.scaling.bin","\n")

#collect the different data types
categorial.columns <- c()
numeric.columns <- c()
boolean.columns <-c()

for (c in names(Selected.data.completecolumns)){
  if (is.numeric(Selected.data.completecolumns[,c]))
    numeric.columns<-c(numeric.columns,c)
  else if (is.logical(Selected.data.completecolumns[,c]))
    boolean.columns<-c(boolean.columns,c)
  else
    categorial.columns<-c(categorial.columns,c)
}

cat("Transforming input data into numeric data\n")
#build a numeric training set through standardised numeric columns and one-hot encoded features
training.data<-data.frame(matrix(data = NA, nrow = dim(Selected.data.completecolumns) [1], ncol = 0))
training.data.colnames<-c()

if(transform.numeric.columns.to.categorial){
  for (c in numeric.columns){
    #standardise the numeric column and add it to the data frame
    nc<-Selected.data.completecolumns[,c]
    if (!grepl(c,"Troph")){
      nc<-log(nc)
      cat("column",c,"log scaled before categorisation\n")
    }
    q<-quantile(nc)
    nc.categ<-c()
    for (n in nc){
      categ<-""
      if (n>q[4])
        categ<-"High"
      else if (n>q[2])
        categ<-"Medium"
      else 
        categ<-"Low"
      nc.categ<-c(nc.categ,categ)
    }
    newcolumnname<-paste0(c,"_categ")
    Selected.data.completecolumns[,newcolumnname]<-nc.categ
    categorial.columns<-c(categorial.columns,newcolumnname)
}
}else{
  for (c in numeric.columns){
    #standardise the numeric column and add it to the data frame
    training.data<-cbind(training.data, scale(Selected.data.completecolumns[,c]))
    training.data.colnames<-c(training.data.colnames,c)
  }
}

sep=":"

for (c in categorial.columns){
  #add one-hot encoded categorial features
  one_hot_encoded <- one_hot_encode(Selected.data.completecolumns[,c])
  training.data<-cbind(training.data, one_hot_encoded[,,drop=F])
  training.data.colnames<-c(training.data.colnames, paste0(c,sep,colnames(one_hot_encoded)) )
}

for (c in boolean.columns){
  #convert T/F to 1/0
  training.data<-cbind(training.data, as.numeric(Selected.data.completecolumns[,c]))
  training.data.colnames<-c(training.data.colnames,c)
}
names(training.data)<-training.data.colnames

#rename the input and out
targetcolumnset<-c(paste0("IUCN_Short",sep,"Threatened"),paste0("IUCN_Short",sep,"Least concern"))
targetcolumnset.orig<-targetcolumnset
if (feature_selection_pca){
##########FEATURE SELECTION##########
cat("\n########PCA FEATURE SELECTION START########\n")
cat("Computing PCA\n")
df.scaled<-as.matrix(training.data[, which(!(names(training.data) %in% targetcolumnset))])
df<-df.scaled
pc <- prcomp(df.scaled,
             center = F,
             scale. = F)

#getting eigenvectors and eigenvalues
eigenvectors<-as.matrix(pc$rotation)
eigenvalues<-pc$sdev
importance<-as.vector(summary(pc)$importance[2,])
eigenperc<-importance

cat("Retrieving only eigenvectors with eigenvalues explaining",eigenvector_threshold*100,"% of the data variance in the PCA space\n")
cumulative.contrib<-0
eigenvalidx<-1
for (i in 1:length(eigenperc)){
  cumulative.contrib = cumulative.contrib+eigenperc[i]
  if (cumulative.contrib>=eigenvector_threshold){
    eigenvalidx=i
    cat("Stopping at eigenvalue",eigenvalidx,"containing cumulative contribution = ",(cumulative.contrib),">=",eigenvector_threshold,"\n")
    break
  }
}
#select the most important eigenvectors
reduced.eigenvectors<-eigenvectors[,1:eigenvalidx]
reduced.eigenperc<-eigenperc[1:eigenvalidx]
if (is.null(dim(reduced.eigenvectors))){
  reduced.eigenvectors<-data.frame(v1=reduced.eigenvectors)
}else{
  reduced.eigenvectors<-data.frame(reduced.eigenvectors)
}

cat("Evaluating information loss after excluding the less important eigenvectors\n")
# projected the original dataset in the reduced PC space
df.proj.reduced <- df.scaled %*% as.matrix(reduced.eigenvectors)
#reproject the dataset backwards
df.proj.reduced.reproj <- df.proj.reduced %*% t(reduced.eigenvectors)
#if information loss was low we should obtain the same data as the original data
error.committed<-mean(apply(abs(df.scaled-df.proj.reduced.reproj),2,mean))
cat("Reprojection error with the reduced space:",error.committed,"\n")

cat("Selecting the features with",feature_threshold*100,"% contribution to the most important eigenvectors\n")
#traspose the reduced eigenvectors and represent them as a data frame
reduced.eigenvectors.features<-t(reduced.eigenvectors)
reduced.eigenvectors.features<-as.data.frame(reduced.eigenvectors.features)
names(reduced.eigenvectors.features)<-names(df)

#project each length vector to the original space to discover its contributions from the original features (=its components in the original space)
reduced.eigenvectors.features.reproj<-sapply(1:dim(reduced.eigenvectors.features)[2], 
                                               function(i){
                                                 column<-reduced.eigenvectors.features[,i]
                                                 column<-column*reduced.eigenperc
                                                 return(column)
                                               },simplify = T)
#produce a data frame of the matrix of vector components in the original space
reduced.eigenvectors.features.weig<-as.data.frame(matrix(reduced.eigenvectors.features.reproj,nrow = dim(reduced.eigenvectors.features)[1],ncol=dim(reduced.eigenvectors.features)[2]))
names(reduced.eigenvectors.features.weig)<-names(df)
#calculate the sum of each components across the eigenvectors in the original space
reduced.eigenvectors.features.abs.sum<-apply(abs(reduced.eigenvectors.features.weig),2,"sum")
#normalise the component sums = contribution of each features to the eigenvectors
reduced.eigenvectors.features.abs.sum.norm<-reduced.eigenvectors.features.abs.sum/sum(reduced.eigenvectors.features.abs.sum)
#sort the contributions
reduced.eigenvectors.features.df <- data.frame(c = reduced.eigenvectors.features.abs.sum.norm)
reduced.eigenvectors.features.df.ordered<-reduced.eigenvectors.features.df[order(reduced.eigenvectors.features.df$c, decreasing=T),, drop=FALSE]

#select the features contributing to the reduced eigenvectors for more than the feature_threshold (e.g. 90%)
cumulative.contrib.f<-0
featurevalidx<-1
for (i in 1:length(reduced.eigenvectors.features.df.ordered$c)){
  
  cumulative.contrib.f = cumulative.contrib.f+reduced.eigenvectors.features.df.ordered$c[i]
  if (cumulative.contrib.f>feature_threshold){
    featurevalidx=i
    cat("Stopping at features",i,"over",length(reduced.eigenvectors.features.df.ordered$c)," containing cumulative contribution = ",(cumulative.contrib.f),">",feature_threshold,"\n")
    break
  }
}

#Report the highest-contribution features
f.names.selected<-as.numeric(row.names(reduced.eigenvectors.features.df.ordered)[1:featurevalidx])
#report the dataset after eliminating all lower-importance features
f.names.selected<-colnames(df[,f.names.selected])
cat("\n####Most important features:\n")
print(as.data.frame(f.names.selected))
f.names.selected<-c(targetcolumnset,f.names.selected)

discarded.features<-names(training.data)[!(names(training.data)%in%f.names.selected)]
cat("\n####Discarded features:\n")
print(as.data.frame(discarded.features))
cat("\n########FEATURE SELECTION FINISHED########\n")
}else{
  f.names.selected<-names(training.data)
}

cat("\n####RF training\n")
#prepare the training set
training.data.important<-training.data[,f.names.selected]
names(training.data.important)[which(names(training.data.important) %in% targetcolumnset.orig)]<-c("o1","o2")
targetcolumnset<-c("o1","o2")
traincolumnset.idx<-which(!(names(training.data.important) %in% targetcolumnset))
traincolumnset<-paste0("i",seq(1,length(traincolumnset.idx)))
names(training.data.important)[traincolumnset.idx]<-traincolumnset
training.data.important<-training.data.important[,c(targetcolumnset,traincolumnset)]
train.1 <- training.data.important
options(warn = -1)
set.seed(20)

#use only one output variable (0/1) for the assessment
f <- as.formula(paste("o1", "~", paste(traincolumnset, collapse = " + ") ))
rf_model <- randomForest(f, data = train.1,ntree=n.trees)
pr.nn1 <- predict(rf_model, newdata = train.1[,traincolumnset])

pr.nn.t<-pr.nn1
pr.nn.t[which(pr.nn1>=thr)]<-1
pr.nn.t[which(pr.nn1<thr)]<-0
pr.nn.t<-as.numeric(pr.nn.t)
original_values <- (train.1[, "o1"]) 
accuracy   <- length(which(pr.nn.t == original_values))/length(pr.nn.t)

cat("Accuracy selftest =", accuracy,"\n")
cat("NOTE - MODEL SAVED TO =","iucn_rf_model.bin","\n")
save(file="iucn_rf_model.bin",rf_model)

#####FEATURE IMPORTANCE ACCORDING TO RF
if (do.characterisation.analysis){
  
  cat("\n#CLASSES CHARACTERISATION")
  feature_importance <- importance(rf_model)
  original_f<-names(training.data[,f.names.selected])[-which(names(training.data[,f.names.selected]) %in% targetcolumnset.orig)]
  rownames(feature_importance)<-original_f
  sorted_features <- as.data.frame(feature_importance[order(-feature_importance[, 1]),])
  sorted_features.norm<-sorted_features/sum(sorted_features)
  
  cumulative<-0
  featureN<-1
  while (cumulative<importance.cumulative.threshold){
    cumulative<-cumulative+sorted_features.norm[featureN,1]
    featureN<-featureN+1  
  }
  
  most.important.features<-rownames(sorted_features.norm)[1:featureN]
  training.data.most.important.selected<-training.data[,c(most.important.features,targetcolumnset.orig)]
  threatened.spp<-training.data.most.important.selected[which(training.data.most.important.selected$`IUCN_Short:Threatened`==1),]
  lc.spp<-training.data.most.important.selected[which(training.data.most.important.selected$`IUCN_Short:Threatened`==0),]
  
  calcCentroid<-function(spp.tab, training.data.most.important.selected,numeric.columns){
    centroid<-c()
    features.centroid<-c()
    for (fe in names(spp.tab)){
      if (fe %in% numeric.columns){
        LOW<-0
        MED<-0
        HIGH<-0
        qfe<-quantile(training.data.most.important.selected[,fe])
        spp.col<-as.numeric(spp.tab[,fe])
        LOW<-length(spp.col[which(spp.col<qfe[3])])
        HIGH<-length(spp.col[which(spp.col>qfe[3])])
        MED<-length(spp.col)-LOW-HIGH
        cat(fe,": L M H: ",LOW,MED,HIGH,"\n")
        interpretation<-"Medium"
        if (LOW>MED && LOW>HIGH)
          interpretation<-"Medium-Low"
        else if (HIGH>MED && HIGH>LOW)
          interpretation<-"Medium-High"
        centroid<-c(centroid,interpretation)
        features.centroid<-c(features.centroid,fe)
      }else{
        spp.col<-as.numeric(spp.tab[,fe])
        ZERO<-length(spp.col[which(spp.col==0)])
        ONE<-length(spp.col[which(spp.col==1)])
        if (ZERO == 0 && ONE == 0){
          ZERO<-length(spp.col[which(spp.col<0.5)])
          ONE<-length(spp.col[which(spp.col>=0.5)]) 
        }
        cat(fe," 0 1: ",ZERO,ONE,"\n")
        
        if (ZERO>ONE)
          interpretation<-"0"
        else if (ZERO<ONE)
          interpretation<-"1"
        else                 
          interpretation<-"Balanced"
        features.centroid<-c(features.centroid,fe)
        if (interpretation == "1"){
          centroid<-c(centroid,fe)
        }else if (interpretation == "0"){
          if (grepl(pattern = ":", x = fe))
              centroid<-c(centroid,gsub(pattern = ":",x = fe,":Not-"))
          else
            centroid<-c(centroid,paste0("Not-",fe))
        }else{
          centroid<-c(centroid,paste0("Equal-",fe))
        }
        
      }
    }
    
    centroid<-as.data.frame(t(centroid))
    names(centroid)<-features.centroid
    return (centroid)
  }
  
  cat("\n#calculating centroid for threatened species:\n")
  centroid.threatened<-calcCentroid(spp.tab=threatened.spp,training.data.most.important.selected,numeric.columns)
  cat("\n#calculating centroid for lc species:\n")
  centroid.lc<-calcCentroid(spp.tab=lc.spp,training.data.most.important.selected,numeric.columns)
  
  characterise<-function(centroid.1,centroid.2){
    characterising<-c()
    #characterize the threatened species
    for (fe in names(centroid.1)){
      characteristic.1<-centroid.1[1,fe]
      if(fe %in% names(centroid.2)){
        characteristic.2<-centroid.2[1,fe]
        if (characteristic.1!=characteristic.2){
          if (grepl(pattern = sep, x = characteristic.1) ){
            characterising<-c(characterising,characteristic.1)
          }else{
            characterising<-c(characterising,paste0(fe,sep,characteristic.1))
          }
        }
      }else{
        if (grepl(pattern = sep, x = characteristic.1) ){
          characterising<-c(characterising,characteristic.1)
        }else{
          characterising<-c(characterising,paste0(fe,sep,characteristic.1))
        }
      }
    }
    return (characterising)
  }
  
  cat("\n#Characterisation of the two classes:\n")
  characterisation.threatened<-characterise(centroid.threatened,centroid.lc)
  
  cat("Characterisation of threatened species:",paste(characterisation.threatened,collapse = ", "),"\n\n")
  characterisation.lc<-characterise(centroid.lc,centroid.threatened)
  
  cat("Characterisation of least concern species:",paste(characterisation.lc, collapse = ", "),"\n\n")
  
  cat("#END OF CLASSES CHARACTERISATION\n\n")
  
  cat("#SEARCHING FOR POLARISATION PATTERNS IN THE DATA\n\n")
  training.data.clustering<-training.data.most.important.selected
  training.data.clustering$classification_score<-pr.nn1
  training.data.clustering.thr<-training.data.clustering[which(training.data.clustering$`IUCN_Short:Threatened`==1),c(most.important.features,"classification_score")]
  training.data.clustering.lc<-training.data.clustering[which(training.data.clustering$`IUCN_Short:Threatened`==0),c(most.important.features,"classification_score")]
  k_start <- 2
  k_max <- 5
  
  #library(ClusterR)
  silhouettes<-Optimal_Clusters_KMeans(training.data.clustering.thr, max_clusters=10, plot_clusters=F, criterion="silhouette")
  k_best<-which(silhouettes==max(silhouettes))
  km<-kmeans(training.data.clustering.thr, centers = k_best)
  km.centers<-as.data.frame(km$centers)
  km.centers<-km.centers[which(km.centers$classification_score>0.8),]
  km.centers<-km.centers[1,]
  centroidInterpreted.thr<-calcCentroid(spp.tab=km.centers,training.data.most.important.selected,numeric.columns)
  cat("Centroid Interpretation for polarised Threatened Species:\n")
  print(centroidInterpreted.thr)
  write.csv(file="Polarised Threatened Species Centroid.csv",x = centroidInterpreted.thr)
  
  silhouettes<-Optimal_Clusters_KMeans(training.data.clustering.lc, max_clusters=10, plot_clusters=F, criterion="silhouette")
  k_best<-which(silhouettes==max(silhouettes))
  km<-kmeans(training.data.clustering.lc, centers = k_best)
  km.centers<-as.data.frame(km$centers)
  km.centers<-km.centers[which(km.centers$classification_score<0.05),]
  km.centers<-km.centers[1,]
  centroidInterpreted.lc<-calcCentroid(spp.tab=km.centers,training.data.most.important.selected,numeric.columns)
  write.csv(file="Polarised Least Concern Species Centroid.csv",x = centroidInterpreted.lc)
  
  cat("Centroid Interpretation for polarised Least Concern Species:\n")
  print(centroidInterpreted.lc)
  
  cat("#END OF SEARCHING FOR POLARISATION PATTERNS IN THE DATA\n\n")
}

###ERROR ANALYSIS: extracts confusing species with different classifications but overlapping feature values
if (do.error.analysis){
  cat("Running error assessment\n")
  error.indices<-(which(pr.nn.t != original_values))
  correct.indices<-(which(pr.nn.t == original_values))
  nproblems<-length(error.indices)
  
  comparisonTable<- data.frame(matrix(ncol = dim(Selected.data.completecolumns.withSpp)[2], nrow = 0))
  for (i in error.indices){
    errvector<-train.1[i,traincolumnset]
    err.interpretation<-pr.nn.t[i]
    #data.with.same.interpr<-train.1[which(train.1$o1 == err.interpretation),traincolumnset]
    data.with.same.interpr<-train.1[,traincolumnset]
    #trace distances
    diff<-sweep(data.with.same.interpr, 2, as.numeric(errvector), "-")
    diff<-diff^2
    diff<-as.numeric(rowSums(diff))
    dd<-data.frame(diff = diff, o1=train.1$o1)
    dd$diff[which(dd$o1!=err.interpretation)]<-.Machine$double.xmax
    dd$diff[error.indices]<-.Machine$double.xmax
    min.j<-which(dd$diff==min(dd$diff))
    closest.vector<-data.with.same.interpr[min.j,]
    
    error.row.completeinfo<-Selected.data.completecolumns.withSpp[i,]
    error.row.completeinfo$classification<-err.interpretation
    closest.original<-Selected.data.completecolumns.withSpp[min.j,]
    closest.original$classification<-err.interpretation
    #added the error and closest correctly-recognised vector to the table
    comparisonTable<-rbind(comparisonTable,error.row.completeinfo)
    comparisonTable<-rbind(comparisonTable,closest.original)
    
  }
  cat("error assessment finished\n")
  #problematic.training.data<-problematic.training.data[order(problematic.training.data$absolute_error, decreasing=T),]
  #write.csv(problematic.training.data,"problematicSpp_rf.csv")
  write.csv(comparisonTable,"problematicComparisonSpp_rf.csv")
}

if (crossvalidate==T){
  # Results from cv
  cat("Cross-validating..\n")
  outs <- NULL
  # Train test split proportions
  #Crossvalidation
  
  for(i in 1:k) {
    cat(i," ")
    index    <- sample(1:nrow(train.1), round(CV.proportion*nrow(train.1)))
    train_cv <- train.1[index, ]
    test_cv  <- train.1[-index, ]
    train_cv <-train_cv[,c(traincolumnset,targetcolumnset)]
    test_cv  <-test_cv[,c(traincolumnset,targetcolumnset)]
    f <- as.formula(paste("o1", "~", paste(traincolumnset, collapse = " + ") ))
    rf_model.cv <- randomForest(f, data = train_cv, ntree=n.trees)
    
    pr.nn1.cv <- predict(rf_model, newdata = test_cv[,traincolumnset])
    pr.nn.t.cv<-pr.nn1.cv
    pr.nn.t.cv[which(pr.nn1.cv>=thr)]<-1
    pr.nn.t.cv[which(pr.nn1.cv<thr)]<-0
    pr.nn.t.cv<-as.numeric(pr.nn.t.cv)
    original_values.cv <- (test_cv[, "o1"])
    accuracy.cv   <- length(which(pr.nn.t.cv == original_values.cv))/length(pr.nn.t.cv)
    outs[i]  <- accuracy.cv
  }
  outs<-round(outs,digits=2)
  accuracy3 <- mean(outs)
  cat("\nCrossvalidation accuracy =",accuracy3,"(min",min(outs),", max",max(outs),")\n\n")
}


