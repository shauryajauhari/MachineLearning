Decision Trees and Random Forests: A Machine Learning Perspective
========================================================
author: Shaurya Jauhari, Mora Lab, GMU. (Email: shauryajauhari@gzhmu.edu.cn)
date: 2019-06-06
autosize: true

RoadMap
========================================================

- Introduction
- General Concepts
  - **B**ootstrapping
  - **Ag**gregating
  - **Bag**ging
- Classification and Regression Trees (CART)
- Decision Trees
- "Best" Split
  - Information Gain
  - Entropy
- Random Forests
- Packages: *rpart*, *party*, *randomForest*


Decision Trees
========================================================

- **Cleave data** into smaller portions to elicit patterns that aid prediction | Recursive Partitioning
- **Logical structures** hence represented can be construed without any *a-priori knowledge*.
<p float="left">
<img src="./Decision_Tree_Example.jpg" width="900" />
<img src="./Decision_Tree_Example_Table.png" width="900" />
</p>


"Best" Split
========================================================
![Condition Check](./If_Then_Template.png)
- Underperformance noted with nominal/ numerical data
- Which attribute to select for tree partitioning?
- Measures of Purity: **Entropy**, **Information Gain**, **Gain ratio**, **Chi-Squared statistic**, etc.
- $Information Gain = Entropy_{Old} - Entropy_{New}$
- Gini Index
- Only Binary classification problems

Gini Index
========================================================
- Calculated for **each node**
- What's the Root node?
-  First calculate how well each feature classifies data
  - Outlook -> Play Tennis; Temperature -> Play Tennis; Humidity -> Play Tennis; Wind -> Play Tennis
  ![Outlook Gini Index calculation](./Gini_Index_Outlook.jpg)

Gini Index (Contd ...)
========================================================
- Perform such quantitation for all features
- "Impure"
  - None of the leaf nodes denote 100% Yes or 100% No for Play Tennis
- Compare impurities
 - less the better
- $$Gini Impurity_{Sunny}$$ = 1 - (the probability of "yes")^2 - (the probability of "no")^2
  - $1 - (2/2+3)^2 - (3/2+3)^2$
  - **0.52**
- Calculate the same for *Rain* and *Overcast*


Gini Index (Contd ...)
========================================================
- Further, calculate the same for each class/ level of the remaining features (**Temperature**, **Humidity**, and **Wind**)
- Classification Counts could be different for different features. That's OK!
![Temperature Gini Index calculation](./Gini_Index_Temperature.jpg)


Gini Index (Contd ...)
========================================================
- Next is, *Weighted average calculation* for a each feature
  - $$Gini Impurity _{Outlook}$$ = (classifications in the leaf node/total classification in the feature) * $$Gini Impurity_{Sunny}$$ + (classifications in the leaf node/total classification in the feature) * $$Gini Impurity_{Rain}$$ + (classifications in the leaf node/total classification in the feature) * $$Gini Impurity_{Overcast}$$


Gini Index (Contd ...)
========================================================
  - $$Gini Impurity _{Outlook}$$ = (5/5+5+4)* 0.52 + (classifications in the leaf node/total classification in the feature) * $$Gini Impurity_{Rain}$$ + (classifications in the leaf node/total classification in the feature) * $$Gini Impurity_{Overcast}$$



Gini Index (Contd ...)
========================================================
- Eventually impurities for each feature will be compared
  - one with least impurity gets to be the root node/ *splitting attribute*
- What about intermediary nodes?


  
Central Dogma of Random Forests
========================================================
![Central Dogma of Random Forests](./Central_Dogma_of_Random_Forests.jpg)
