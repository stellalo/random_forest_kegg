<h1>Random Forest Kegg Analysis Scripts</h1>

<h2>👩🏻‍💻 Description</h2>
This repository contains scripts that analyze the KEGG gene expression of 1154 yeast species using ML algorithm Random Forest in R. The objective of this project is to identify the KEGGs that have the highest possiblity of predicting pathogenic yeast using different characteristics such as survival at temperature 37 Celcius and resistance to 9 anti-fungal treatments.
<br />

<h2>🪐 Language and Packages</h2>

- <b>R Markdown</b>

The following packages are required to run the scripts: 
```ruby
library(reshape2)
library(dplyr)
library(tidyverse)
library(randomForest)
library(caret)
library(ranger)
```

<h2>🦠 RF_temp.Rmd Walk-Through</h2>
This script runs the RF (random forest) algorithm on yeast species' KEGG gene expression using survival at temperature 37+ as environmental factor. We first perform RF without any tuning, then perform undersampling to improve OOB err rate. 
<br />
<br />

- First, read files that contains all tAI values (codom optimization values) and the kegg annotation stored in "eTAI_kegg_wide.txt" and the file with temperature information stored in "thirtyseven_data.txt"
<br />

```ruby
kegg <- read.delim("etAI_kegg_wide.txt",sep = "\t",header = T)
env <- read.delim("thirtyseven_data.txt",header = T)
```
<br/>

- We want to remove columns with more than 80% of N/A data from the kegg data frame and change the remaining N/As to 0. Once done, combine the 2 data frames by species id. Clean the combined data frame by removing species with missing temperature data, then remove rows that all are N/As but have temperature value:
<br/>

```ruby
kegg <- kegg[, na_proportion < 0.8]
kegg <- kegg %>% replace(is.na(.), 0)
comb_mat <- full_join(kegg, env,by = c("spec" = "assembly_fullID_updated"))
comb_mat <- comb_mat[complete.cases(comb_mat$X37_EVIDENCE), ]
comb_mat <- comb_mat[complete.cases(comb_mat[,-ncol(comb_mat)]), ]
```
<br/>

- Make temperature data a "factor":

```ruby
comb_mat$X37_EVIDENCE <- as.character(comb_mat$X37_EVIDENCE)
comb_mat$X37_EVIDENCE <- as.factor(comb_mat$X37_EVIDENCE)
```

<br/>

- Run random forest without tuning. Set training data to 70% and testing data to 30% of all data:

```ruby
set.seed(120)
ind <- sample(2, nrow(comb_mat), replace = TRUE, prob = c(0.7, 0.3))
train <- comb_mat[ind==1,]
test <- comb_mat[ind==2,]

rf <- randomForest(X37_EVIDENCE~., data=train, proximity=TRUE, importance=TRUE)
```

<br/>

- Now plot the data from the RF algorithm:
  
```ruby
plot(rf)
p1 <- predict(rf, train)
p1.cm<-confusionMatrix(p1, train$X37_EVIDENCE)
View(p1.cm$byClass)

p2 <- predict(rf, test)
p2.cm<-confusionMatrix(p2, test$ X37_EVIDENCE)
View(p2.cm$byClass)
write.table(p2.cm$table, file="p2.cm.temp.txt")


hist(treesize(rf),
     main = "No. of Nodes for the Trees",
     col = "green")
#Variable Importance
varImpPlot(rf,
           sort = T,
           n.var = 10,
           main = "Top 10 - Variable Importance")
```

- Resulting graphs & tables: <br/>

<br/>

<b>OOB Error</b>

<img width="161" alt="oob_err" src="https://github.com/stellalo/random_forest_kegg/assets/89308696/3cad9cd9-8a1e-4772-87dc-e020c3d37b76">


![oob](https://github.com/stellalo/random_forest_kegg/assets/89308696/230b9193-4c3a-46ec-945f-aa7a594ed1a0)
![kegg](https://github.com/stellalo/random_forest_kegg/assets/89308696/364bdffe-fdd2-4163-9617-c1a4ed948d08)
![tree](https://github.com/stellalo/random_forest_kegg/assets/89308696/03ba59f0-b1c9-4591-bae2-bec9ddc6b14f)

<br/>

- Since the oob error is too high (~0.75 for testing set), we are going to perform undersampling to balance the data set. Once undersampling is completed, perform RF again on the balanced data set:

<br/>

```ruby
data_balanced <- downSample(x = comb_mat[, -1], y = comb_mat$X37_EVIDENCE)

set.seed(123) 
train_index <- createDataPartition(data_balanced$X37_EVIDENCE, p = 0.7, list = FALSE)
train_balanced <- data_balanced[train_index, ]
test_balanced <- data_balanced[-train_index, ]

# Train the random forest classifier on the balanced training set
rf_balanced <- randomForest(X37_EVIDENCE~., data=train_balanced, proximity=TRUE, importance=TRUE)
```

- Resulting graphs & tables after undersampling: <br/>

<b>OOB Error</b>

<img width="175" alt="oob_err1" src="https://github.com/stellalo/random_forest_kegg/assets/89308696/c42451c1-f94d-45bc-8e22-45d8b36a557f">
  
![oob1](https://github.com/stellalo/random_forest_kegg/assets/89308696/1dbf9db6-8d10-435c-9885-c0a4c2ee6995)
![kegg1](https://github.com/stellalo/random_forest_kegg/assets/89308696/e19c3c8e-d176-48bd-b14e-4cd37b0b1244)
![tree1](https://github.com/stellalo/random_forest_kegg/assets/89308696/25a61bc5-5e01-4eb6-a85a-34ea921e7143)

<br/>

<h2>🦠 RF_loop.Rmd Walk-Through</h2>
This script is very similar to RF_temp, but it loops through 100 seeds and then runs the RF model 100 times.

<br/>

- We first generate 100 non-repeated seed values and store it in a vector:

<br/>

```ruby
used_seeds <- c()
for(j in 1:100){
  #Generate a random seed value
  repeat {
  seed <- sample(1:1000, 1)
  if (!(seed %in% used_seeds)) {
    used_seeds <- c(used_seeds, seed)
    break
   }
  }
```

<br/>

- We then pass the 100 non-repeated seed values to the function set.seed() and run the RF algorithm after downsampling the dataset. Note that here, we are using the tuneRF function which searches for the optimal mtry value and setting doBest = TRUE to run the RF with the optimal mtry found.

<br/>

```ruby
  # Set the seed
  set.seed(seed)
  
  rf <-tuneRF(train_balanced[, -ncol(train_balanced)], train_balanced$X37_EVIDENCE, ntreeTry = 500, stepFactor = 1.5, improve = 0.01, trace = FALSE, plot = FALSE,doBest = TRUE,proximity=TRUE, importance=TRUE)
```

<br/>

<h2>🦠 RF_resistance.Rmd Walk-Through</h2>
This script is very similar to both RF_loop and RF_temp. It performs undersampling, loops through RF 100 times using the optimal mtry found by tuneRF, but the environmental (factor) variable changed from survival at temperature 37 to resistance-related characteristic. 

<br/>

```ruby
env <- read.delim("resistance.txt",header = T)
```
<br/>

<h2>🦠 KEGG_leave_out Walk-Through.Rmd Walk-Through</h2>
This script performs a "leave out" analysis, which uses a while-loop that only stops the loop when the resulting OOB error rate is less than 50%. In each loop iteration, the most important KEGG value (judged on MeanDecreaseAccuracy rank) is taken out of the balanced training & testing data frame and concatenated to a data frame called "kegg_left_out". The script returns the list of KEGGs that have the highest value in predicting a yeast specie's possiblity of survival in temperature 37+. 
<br/>
<br/>

```ruby
set.seed(230)
#data frame to store KEGGs
KEGG_left_out <- data.frame()
#data frame to store all oobs (mean)
all_mean_oob <- data.frame()
#set variable to 0
oob <- 0

#Loop until OOB >= 0.5
while (oob < 0.50){
  #run RF
  rf <- randomForest(X37_EVIDENCE~., data=train_balanced, proximity=TRUE, importance=TRUE)
  
  #Relative Importance of all KEGGs (MeanDecreaseAccuracy)
  this_importance<-rf$importance
  sort_indices <- order(this_importance[, "MeanDecreaseAccuracy"], decreasing = TRUE)
  sorted_KEGG <- this_importance[sort_indices,]
  
  #remove the most important KEGG from data frame (use Accuracy) and concatenate it to Kegg_left_out
  KEGG_left_out <- rbind(KEGG_left_out, rownames(sorted_KEGG)[1])
  
  #remove that most importance KEGG from data frame
  train_balanced <- train_balanced[,!colnames(train_balanced) == rownames(sorted_KEGG)[1]]
  test_balanced <- test_balanced[,!colnames(test_balanced) == rownames(sorted_KEGG)[1]]
  #reset variable OOB
  rf_conf<-rf$confusion
  oob <- mean(rf_conf[1,3],rf_conf[2,3])
  #store oob
  all_mean_oob <- rbind(all_mean_oob,oob)
}
```
