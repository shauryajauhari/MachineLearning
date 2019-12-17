## Installing packages and loading libraries.
setwd("C:/Users/rajni/Desktop/Machine_Learning/Machine_Learning_Deep_Learning")
install.packages("keras", 
                 repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
                 dependencies = TRUE)
library(keras)

## The following function helps to setup the "TensorFlow" work environment and 
## points any necessary dependencies to be installed. Note that the directions are
## subjective and may vary from system to system.

install_keras()

# Turn `enhancer prediction data` into a matrix
d <- as.matrix(d)

# Set iris `dimnames` to `NULL`
dimnames(d) <- NULL
str(d)

## Data manipulation


d3 <- as.matrix(d)
dimnames(d3) <- NULL
str(d3)


# Normalize, function is a part of keras package
d3[,1:4] <- normalize(d3[,1:4]) # normalize the features #
d3[,5] <- as.numeric(d3[,5]) # the class variable must be numeric, use of which will be known in a while #



# Data Partitioning

set.seed(1234)
ind <- sample(2, nrow(d3),replace=T, prob=c(0.7,0.3))
training <- d3[ind==1,1:4] 
trainlabels <- d3[ind==1,5] 
testing <- d3[ind==2,1:4] 
testlabels <- d3[ind==2,5]

# One hot encoding

traininglabels <- to_categorical(trainlabels)
testinglabels <- to_categorical(testlabels)
print(testinglabels)

# Create sequential model
dlmodel <- keras_model_sequential()

dlmodel %>% layer_dense(units=8,
                        activation='relu',
                        input_shape=c(4)) %>% 
  layer_dense (units=2,
               activation='softmax')
summary(dlmodel)

# Model compilation

dlmodel %>% compile (loss='binary_crossentropy',
                     optimizer='adam',
                     metrics='accuracy')

# Fitting the model with training data

dlmodeltrain <- dlmodel %>% fit (training, 
                                 traininglabels,
                                 epoch=10, # number of iterations #
                                 batch_size=10, # default=32 # 
                                 validation_split=0.2) # 20% of the data #

plot(dlmodeltrain)


# Evaluating the model with test data

dlmodeltest <- dlmodel %>% evaluate(testing,testinglabels,batch_size=100) 
print(dlmodeltest)

# Prediction and Confusion Matrix from test data

prob <- dlmodel %>% predict_proba(testing)
pred <- dlmodel %>% predict_classes(testing)
table1 <- table(Predicted =pred, Actual= testinglabels)
table1
cbind(prob, pred, testinglabels)


