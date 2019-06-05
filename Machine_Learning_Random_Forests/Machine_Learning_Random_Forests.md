Decision Trees and Random Forests: A Machine Learning Perspective
========================================================
author: Shaurya Jauhari, Mora Lab, GMU. (Email: shauryajauhari@gzhmu.edu.cn)
date: 2019-05-31
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
![Decision Tree Example](./Decision_Tree_Example.jpg)
![Decision Tree Example Table](./Decision_Tree_Example_Table.png)

"Best" Split
========================================================
![Condition Check](./If_Then_Template.png)
- Underperformance noted with nominal/ numeric data
- Which attribute to select for tree partitioning?
- Measures of Purity: **Entropy** and **Information Gain**
- $Information Gain = Entropy_{Old} - Entropy_{New}$ 

Entropy
========================================================
- 