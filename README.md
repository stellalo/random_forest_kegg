<h1>Random Forest Kegg Analysis Scripts</h1>

<h2>👩🏻‍💻 Description</h2>
This repository contains scripts that analyze the KEGG gene expression of 3000+ yeast species using ML algorithm Random Forest in R.
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
<p align="center">
  Result:<br/>
  <img src="https://imgur.com/a/dZll4Vi"/>
  <br/>
  

<p align="center">
Launch the utility: <br/>
<img src="https://i.imgur.com/62TgaWL.png" height="80%" width="80%" alt="Disk Sanitization Steps"/>
<br />
