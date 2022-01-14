## **(1)**

The package poLCA has an example dataset for latent class clustering with categorical variables that are not just binary called "election".
The categorical variables to be clustered are variables 1-12
The dataset has missing values. There are two different ways of handling them in
the latent class analysis. The first one is to only use the 1311 observations (out of
1785) that do not have missing values.
The second way is to define a new category for the missing values, i.e., replacing
all missing values NA by a category called "NA". Note that all the variables are of
type "factor", and this requires to define a new factor level.
Run the following clusterings and compare them using MDS plots based on the
simple matching distance. Use two different MDS outputs, one for clusterings computed on electioncomplete and one for clusterings computed on electionwithna.
Also compute ARIs for every pair of clusterings. Where you compare one clustering
based on electioncomplete with 1311 observations and one clustering based on
electionwithna with 1785 observations, only use the 1311 observations without
missing values for the ARI computation.

  * Compute a latent class clustering with 3 clusters using poLCA. This will automatically only use
complete cases, so running it as on the help page is equivalent to running it on
electioncomplete.
  * Compute a latent class clustering with 3 clusters using poLCA on the electionwithna
data.
  * Compute a latent class clustering with 3 clusters using flexmxedruns on the
electioncomplete data.
  * Compute a latent class clustering with 3 clusters using flexmixedruns on the
electionwithna data.
  * Compute a distance-based clustering of your choice with 3 clusters based on
the simple matching distance on the electioncomplete data.
  * Compute the same distance-based clustering with 3 clusters on the electionwithna
data.
  * Just another way of handling missing values is to use the simple matching
distance, and to average the distance for every pair of observations only over
those variables with both observations non-missing, i.e., the way missing values
are handled in the Gower coefficient. Compute a dissimilarity-based clustering of your choice with 3 clusters using
this dissimilarity on the election12 data with missing values.
  * Either on electioncomplete or on electionwithna (your choice), compute a
latent class clustering using flexmxedruns with estimated number of clusters.

## **(2)**

For two different latent class clusterings computed in question 1 on the
electionwithna-data produce heatmaps.
Comment on the plots. Do you find the clusters convincing? Why or why not? Is
there evidence against local independence?

## **(3)**

* Assume a situation with 10 categorical variables. Five variables are binary,
three variables have three categories, and two variables have five categories.
What is the number of free parameters for
1. a general categorical model that models all possible probabilities,
2. a latent class mixture model with 4 mixture components?

* Consider a situation where observations X = (X1;X2;X3) have 3 categorical
variables, all binary with values 0 and 1. Consider a latent class mixture with
three components (clusters), and the following parameters:

Mixture proportions  π1 =  π2 =  π3 = 1

**Cluster 1**: η111 = 0.1; η211 = 0.1; η311 = 0.1

**Cluster 2**: η112 = 0.5; η212 = 0:5; η312 = 0.5

**Cluster 3**: η113 = 0.8; η213 = 0.8; η313 = 0.8

Compute the probabilities for all possible observations, i.e., all eight combina-
tions of results 0 and 1 in the three variables resulting from this. Also compute
the degrees of freedom that this model has, and the degrees of freedom of a
general categorical model for this situation.

