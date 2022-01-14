## (1)
Run K-means for the Olive Oil data with K = 3 and K = 9 with scaled and unscaled data. Assuming that the macro-areas are the "true" clusters for K = 3, use table to compare the macro-areas with the clustering. Do you think that these
are good clustering results in terms of matching the macro-areas? Why? Does the clustering on scaled or unscaled data look better? Do the same for the regions and the K = 9-clustering.

## (2)
On Moodle you can find the data set Boston.dat. This data set contains information collected by the U.S Census Service concerning housing in the area of Boston Mass. The data was originally published by
> Harrison, D. and Rubinfeld, D.L. 'Hedonic prices and the demand for clean air', J. Environ. Economics & Management, vol.5, 81-102, 1978
 
Every observation refers to a census tract, i.e., a town or district. The data contains the following columns:
  * **crim** per capita crime rate by town.
  * **zn** proportion of residential land zoned for lots over 25,000 sq.ft.
  * **indus** proportion of non-retail business acres per town.
  * **chas** Charles River dummy variable (= 1 if tract bounds river; 0 otherwise).
  * **nox** nitrogen oxides concentration (parts per 10 million).
  * **rm** average number of rooms per dwelling.
  * **age** proportion of owner-occupied units built prior to 1940.
  * **dis** weighted mean of distances to ve Boston employment centres.
  * **rad** index of accessibility to radial highways.
  * **tax** full-value property-tax rate per $10,000.
  * **ptratio** pupil-teacher ratio by town.
  * **black** `1000(Bk-0.63)^2` where Bk is the proportion of blacks by town.
  * **lstat** lower status of the population (percent).
  * **medv** median value of owner-occupied homes in $1000s.
Visualise the data, produce a clustering of this data set that looks reasonable to you, and explain the reasons why you have chosen this and you think it is reasonable.

## (3)
*kmeans++* is the name of a method to initialise the k-means algorithm that has been proposed in the literature (for really big datasets it may be problematic
to run Lloyd's or similar algorithms a lot of times from random starting points, and having just one well picked starting point will be much faster and hopefully not
much worse or even better). Do some research on the internet, find out and explain how this works. It can be run by the following function `kmpp`, where X is a data
matrix and k is the number of clusters.
Run this on one or more data sets on which you have also run kmeans, and compare the achieved values of the objective function with kmeans.
