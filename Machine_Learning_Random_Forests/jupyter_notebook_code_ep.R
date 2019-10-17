library(party)
library(rpart)
library(rpart.plot)

epdata <- readRDS("../Machine_Learning_Deep_Learning/data/ep_data.rds")
rownames(epdata) <- c()


set.seed(3)
data_partition <- sample(2, nrow(epdata), replace = TRUE, prob = c(0.7,0.3))
train <- epdata[data_partition==1,]
test <- epdata[data_partition==2,]


cat("We have",nrow(train),"training examples and",nrow(test),"for testing.")


ep_tree <- ctree(class ~ peaks_h3k27ac + peaks_h3k4me3 + peaks_h3k4me2 + peaks_h3k4me1, data= train)
print(ep_tree)

plot(ep_tree)

tree_pred <- predict(ep_tree,test)
tree_pred <- ifelse(tree_pred>0.5, "1", "0")

tree_pred_prob <- predict(ep_tree,test, type="prob")
head(tree_pred_prob)

tab <- table(tree_pred, test$class)
print(tab)

accur <- 1 - sum(diag(tab))/sum(tab)
cat("The accuracy of the model is",accur*100,"%")
#print("\n")
ifelse(accur>90, print("Not bad."), print("Poor."))

ep_tree1 <- rpart(class ~ peaks_h3k27ac + peaks_h3k4me3 + peaks_h3k4me2 + peaks_h3k4me1, data= train)
print(ep_tree1)


rpart.plot(ep_tree1)

tree_pred1 <- predict(ep_tree1,test)
head(tree_pred1)


library(randomForest)
library(caret)


epdata_sample <- readRDS("../Machine_Learning_Deep_Learning/data/ep_data_sample.rds")
rownames(epdata_sample) <- c()
## Data partitioning

set.seed(5)
data_partition <- sample(2, nrow(epdata_sample), replace = TRUE, prob = c(0.7,0.3))
train_sample <- epdata_sample[data_partition==1,]
test_sample <- epdata_sample[data_partition==2,]


rf1 <- randomForest(formula= class ~ peaks_h3k27ac + peaks_h3k4me3 + peaks_h3k4me2 + peaks_h3k4me1, 
                    data= train_sample,
                    proximity = TRUE,
                    ntree=500) # default
print(rf1)


# Prediction and Confusion Matrix
p1<- predict(rf1,train_sample)
head(ifelse(p1>0.5, "Non-Enhancer", "Enhancer"))
p1 <- ifelse(p1>0.5, "1", "0")
table(p1, train_sample$class)

# Alternatively,
confusionMatrix(as.factor(p1),as.factor(train_sample$class))

p2<- predict(rf1,test_sample)
head(ifelse(p2>0.5, "Non-Enhancer", "Enhancer"))
p2 <- ifelse(p2>0.5, "1", "0")
table(p2, test_sample$class)

# Alternatively,
confusionMatrix(as.factor(p2),as.factor(test_sample$class))

plot(rf1)


t <- tuneRF(train_sample[,-5], train_sample[,5],
            stepFactor =2,
            plot = TRUE,
            ntreeTry = 100,
            trace = TRUE,
            improve = 0.005)               

rf2 <- randomForest(class~.,
                    data = train_sample,
                    ntree=250,
                    mtry=1,
                    importance= TRUE,
                    proximity = TRUE)
print(rf2)


rf2 <- randomForest(class~.,
                    data = test_sample,
                    ntree=250,
                    mtry=1,
                    importance= TRUE,
                    proximity = TRUE)
print(rf2) 


p3<- predict(rf2,test_sample)
head(ifelse(p3>0.5, "Non-Enhancer", "Enhancer"))
p3 <- ifelse(p3>0.5, "1", "0")
confusionMatrix(as.factor(p3),as.factor(test_sample$class))


# Number of nodes for the trees
library(ggplot2)
ts <- as.data.frame(treesize(rf2), row.names = c())
ggplot(data = ts, aes(x = `treesize(rf2)`)) + 
  geom_histogram(binwidth = 5, color="purple", fill="lavender", position="identity", alpha=0.7)+
  labs(x="Trees", y="Sizes")+
  ggtitle("Number of Nodes in the Trees")+
  geom_density(alpha=0.6)+
  theme_light()+
  theme(plot.title= element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 12))

# Variable Importance

varImpPlot(rf2)

varImpPlot(rf2, sort = TRUE, n.var = 2, main = "Top 2 variables")

importance(rf2)

varUsed(rf2)

partialPlot(rf2, train_sample, peaks_h3k27ac , "Enhancer")

partialPlot(rf2, train_sample, peaks_h3k4me2, "Non-Enhancer")

# Extract Single tree
getTree(rf2, k=1)

#Multidimensional Scaling (MDS) plot of proximity matrix.

MDSplot(rf2, train_sample$class)




















