# Load 'randomForest' package
install.packages("randomForest")
library(randomForest)

# Execute random forest algorithm
rf <- randomForest(formula= Class~.,data= train)
print(rf)

# Prediction and Confusion Matrix
install.packages("caret")
library(caret)

p1<- predict(rf,test)
confusionMatrix(p1,test$Class)

# Error rate of random forest
plot(rf)

# Tuning mtry

t <- tuneRF(train[,-5], train[,5],
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 120,
            trace = TRUE,
            improve = 0.005)


rf2 <- randomForest(Species~.,
                    data = train,
                    ntree=120,
                    mtry=4,
                    importance= TRUE,
                    proximity = TRUE)
print(rf2)                    

# Number of nodes for the trees
hist(treesize(rf2), col = "green", main = "Number of nodes for the trees")

# Variable Importance

varImpPlot(rf2)

varImpPlot(rf2, sort = TRUE, n.var = 2, main = "Top 2 variables")

importance(rf2)
varUsed(rf2)

# Partial Dependence Plot

partialPlot(rf2, train, Petal.Width, "setosa")
partialPlot(rf2, train, Sepal.Width, "setosa")


# Extract Single tree
getTree(rf2, k=1)

# Multidimensional Scaling (MDS) plot of proximity matrix.

MDSplot(rf2, train$Species)
