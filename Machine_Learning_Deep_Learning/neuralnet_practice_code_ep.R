## Install the package and load the library

install.packages("neuralnet", 
                 repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
                 dependencies = TRUE)
library(neuralnet)


## Let us consider our example dataset
d <- readRDS("./data/ep_data_sample.rds")
str(d)

# Classes as factor levels. 
d$class <- as.numeric(as.factor(d$class))-1

# Min-Max normalization
d$reads_h3k27ac <- (d$reads_h3k27ac-min(d$reads_h3k27ac, na.rm = T))/(max(d$reads_h3k27ac, na.rm = T)-min(d$reads_h3k27ac, na.rm = T))
d$reads_h3k4me3 <- (d$reads_h3k4me3-min(d$reads_h3k4me3, na.rm = T))/(max(d$reads_h3k4me3, na.rm = T)-min(d$reads_h3k4me3, na.rm = T))
d$reads_h3k4me2 <- (d$reads_h3k4me2-min(d$reads_h3k4me2, na.rm = T))/(max(d$reads_h3k4me2, na.rm = T)-min(d$reads_h3k4me2, na.rm = T))
d$reads_h3k4me1 <- (d$reads_h3k4me1-min(d$reads_h3k4me1, na.rm = T))/(max(d$reads_h3k4me1, na.rm = T)-min(d$reads_h3k4me1, na.rm = T))

# Data Partition
set.seed(100)
ind <- sample(2, nrow(d), replace = TRUE, prob = c(0.7, 0.3))
training <- d[ind==1,]
testing <- d[ind==2,]

# Neural Networks
set.seed(005)
nn <- neuralnet(class~reads_h3k27ac+reads_h3k4me3+reads_h3k4me2+reads_h3k4me1,
                data = training,
                hidden = 9,
                err.fct = "sse",
                act.fct = "logistic",
                linear.output = FALSE)
plot(nn)

# Prediction
output <- compute(nn, training[,-5])
head(output$net.result)
head(training[1,])

# Node Output Calculations with Sigmoid Activation Function
in5 <- -0.27206 + (0.7777778*-6.9519)+(0.4166667*1.11426)+(0.8305085*1.25498)+(0.8333333*1.02952)
print(in5)
out5 <- (1/1+exp(-in5))
print(out5)

# Confusion Matrix & Misclassification Error - training data
output <- compute(nn, training[,-5])
p1 <- output$net.result
pred1 <- ifelse(p1>0.5, 1, 0)
tab1 <- table(pred1, training$class)
tab1
1-sum(diag(tab1))/sum(tab1)

# Confusion Matrix & Misclassification Error - testing data
output <- compute(n, testing[,-1])
p2 <- output$net.result
pred2 <- ifelse(p2>0.5, 1, 0)
tab2 <- table(pred2, testing$class)
tab2
1-sum(diag(tab2))/sum(tab2)


nn1 <- neuralnet(class~Sepal.Length+Sepal.Width+Petal.Length+Petal.Width,
                 data = training,
                 hidden = c(2,3),
                 err.fct = "ce",
                 act.fct = "logistic",
                 linear.output = FALSE)

nn2 <- neuralnet(class~Sepal.Length+Sepal.Width+Petal.Length+Petal.Width,
                 data = training,
                 hidden = 5,
                 err.fct = "ce",
                 act.fct = "logistic",
                 linear.output = FALSE,
                 lifesign='full',
                 rep=5)