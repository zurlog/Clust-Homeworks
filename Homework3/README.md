## **(1)**

Consider the following dataset with n = 4 observations and p = 5 variables, the first of which is categorical (for use with the simple matching distance),
the second, third, and fourth are binary (Jaccard distance should be used), and the
fourth is on a continuous scale. "NA" denotes missing values.

x1 = (blue; 1; 1; 0; 12)

x2 = (red; 0; 0; NA; NA)

x3 = (red; 1; 0; NA; 17)

x4 = (green; 1; 0; 0; 21)

What are the Gower (coefficient) dissimilarities between all pairs of observations?
  * Manually compute Gower dissimilarities based on distances for all variables
separately, including variables 2-4.
  * Compute Gower dissimilarities treating variables 2-4 as a single group on which
you compute a Jaccard dissmilarity that has weight 3
(because of 3 variables) in the final Gower dissimilarity. Does this give the
same result as part (a)?
  * Compute the Gower dissimilarites using the daisy-function in R and check
against the manual calculation in (a) and (b).

## **(2)**

Give counterexamples to show that the correlation dissimilarity and the Gower coefficient do not fulull the triangle inequality,
i.e., in each case present three observations of which you show that they violate the triangle inequality.

## **(3)**

Consider two p-dimensional binary observations x1; x2 without missing
values. Prove that the Jaccard distance dJ(x1; x2) is the same as the Gower coefficient dG(x1; x2), computed with wl = 1 for l = 1,...,p and handling of joint absences.

## **(4)**

On Virtuale you can find the data set `covid2021.dat`. This data set
has time series characterising the spread of Covid-19 in 179 countries. The data are
from [Github](https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series).
The time span is 1 April 2020 to 7 October 2021. Data give for each day the number
of additional cases in the previous week (to remove weekday effects) divided by the
country's population (in 1,000).

The data set has 559 variables:
* **country** Name of country

* **continent** The continent to which the country belongs

* **latitude** latitude

* **longitude** longitude

* **X4.1.20-X10.7.21** 555 daily variables giving the counts. Only these should be
used for clustering, the others can be used for interpretation and visualisation.

The task here is to cluster the countries in order to find groups of countries with
similar developments. Try out one or more dissimilarity-based hierarchical clustering
methods together with Euclidean and correlation dissimilarity. You may try to come
up with further ideas for defining a dissimilarity for these data. Choose a number
of clusters, try to understand and interpret the clusters as good as you can, using
the information in the data (using visualisation as you see fit), and built yourself an
opinion which of the tried out clusterings is most appropriate, and how appropriate
they are in general (you may not be happy with any of them, in which case you may
think about what went wrong, and how a better clustering could be achieved)
