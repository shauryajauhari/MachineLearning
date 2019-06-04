FROM rocker/r-base:latest
FROM rocker/verse:3.5.1
USER rstudio
WORKDIR /Users/mei/Desktop/Machine_Learning_Logistic_Regression/Logistic_Regression_Lab
RUN export PATH=$PATH:/Applications/RStudio.app/Contents/MacOS/pandoc
COPY ./glmnet_stats_glm_demo.Rmd glmnet_stats_glm_demo.Rmd
RUN R -e "rmarkdown::render('glmnet_stats_glm_demo.Rmd')"
