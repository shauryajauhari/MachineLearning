setwd("C:/Users/rajni/Desktop/Machine_Learning/Machine_Learning_Deep_Learning")
install.packages("neuralnet", 
                 repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
                 dependencies = TRUE)
library(neuralnet)
hist(as.numeric(my_consolidated_dataset$Class))
