
# IUCN Classifier

A script using machine learning to predict the **IUCN** **threatened** or **least concern** status based on FishBase data. 

## Application

 1. Download the IUCNclassifier_application.R script
 2. Download the **iucn_rf_model.bin** model file and **data.scaling.bin** data-scaling file and put them in the same folder as the script
 3. Load the script or execute it by changing the input file name. E.g., at the end of the IUCNclassifier_application.R script you can see:

   > library(randomForest) 
   > input.file<-"AAA_RF_IUCN_240115_test.csv"
   > iucn_classification(input.file = input.file)

 4. The input file name should follow the template of the [AAA_RF_IUCN_240115_test.csv](https://github.com/cybprojects65/IUCNClassifier/blob/main/AAA_RF_IUCN_240115_test.csv) file
 5. The mandatory fields of the test file are:
> "Order", "PD50", "Fresh", "Brack", "Saltwater", "BodyShapeShort",
> "DemPelShort", "MaxLength", "Troph", "AquariumShort", "GameFish",
> "Importance", "Protected", "ClimateZone", "Resilience", "RepGuild1",
> "CountOfC_Code", "CountOfAreaCode"

 7. The script will produce a new CSV file (with the "_iucn_ml_assessed.csv" suffix) containing the additional column "threatened" containing a binary IUCN status assessment (0=Least Concern; 1=Threatened).

## Training

 1. Download the IUCNclassifier_vX_random_forest.R script
 2. Change the input data in

    Reference.data <- read.csv("AAA_RF_IUCN_240115.csv")

 3. Check that mandatory fields are provided:
> "Order", "PD50", "Fresh", "Brack", "Saltwater", "BodyShapeShort",
> "DemPelShort", "MaxLength", "Troph", "AquariumShort", "GameFish",
> "Importance", "Protected", "ClimateZone", "Resilience", "RepGuild1",
> "CountOfC_Code", "CountOfAreaCode"
	Plus the "IUCN_Short" column containing the "Threatened" or "Least concern" value. See the [reference example](https://github.com/cybprojects65/IUCNClassifier/blob/main/AAA_RF_IUCN_240115.csv).

 4. Execute the script. It will produce a model file (**iucn_rf_model.bin**) and a data-scaling reference (**data.scaling.bin**) to be used for application.

**Parameter adjustements:**

 1. To explore the features contribution to the variance in the data (i.e., those carrying most of the information) you can activate Principal Component Analysis:

> feature_selection_pca<-F

 2. To explore the features characterising each IUCN class and detect small groups (clusters) of species with well defined feature combinations associable with the class, you can enable the characterisation analysis:

 > do.characterisation.analysis<-F
    
3. Numeric features can be transformed into categorial features to semantically align with the other categorial features. This usually results in performance loss due to information flattening. However, it can be activated through:

> transform.numeric.columns.to.categorial<-F
4. Error analysis selects the species that were wrongly classified during the model training and self-testing, and produces a CSV file (problematicComparisonSpp_rf.csv) containing, for each row, the misclassified species data and, in the next row, the most similar species correctly classified in the complementary class of the first species. This analysis can be acrivated through

> do.error.analysis<-F
5. The script internally uses Random Forest, a model that randomly extracts many subgroups of features and, for each, creates a decision tree, whose nodes represent questions about one feature at a time (it is 0 or 1? Is it high/low? etc.). The "leaves" are decisions about the species' "Threatened" or "Least Concern" status. The current Random Forest model generates 500 decision trees (the optimal number associated with the highest cross-validation accuracy) and decides the species' status based on the largest consensus among the decision trees. The Random Forest learns to optimise the decision trees based on the feature values associated with the examples (training set). It can later apply these trees to previously unseen features (test set) to generate new decisions.
The Random Forest parameters can be adjusted through the following lines:
> n.trees       <-500 #number of decision trees from which an ensemble assessment should be extracted

and the performance assessment can be regulated through the following parameters:

   > crossvalidate <- T  #enable/disable cross validation
   > k             <- 20 # k-fold cross validation
   > thr           <-0.5 # dichotomic decision threshold on the RF output to distinguish between threatened (>thr) or least concern (<thr) species


	
