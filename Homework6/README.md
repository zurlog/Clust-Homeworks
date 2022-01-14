## **(1)**

On the Virtuale page of the course you find the dataset stars5000.dat.
The data contain information about 5000 celestial objects of which spectra were obtained by mounting a prism in front of a telescope.
The aim is to fiind classes of celestial objects with specific characteristics. There are
six variables containing summary information of the spectra (original spectra con-
tain 300 highly dependent variables; the original database contains about 4 million
objects). These are:

  * **casn** Signal-to-noise ratio of the Calcium break
  * **cacont** Contrast of the Calcium break to smoothed version of the spectrum
  * **kl1** First principal component of smoothed spectrum
  * **kl2** Second principal component of smoothed spectrum
  * **xh1** "Half power point" in upper spectrum
  * **xh2** "Half power point" in lower spectrum
  
Cluster these data using Gaussian mixtures, t-mixtures, skew-normal mixtures, and
skew-t mixtures, and decide which clustering you find most convincing, with reasons.
Although methods with flexible covariance/shape matrices can in principle handle
variables with very different variances, value ranges here are vastly different, and
standardisation may help, maybe in a robust manner (using median and MAD) because of the presence of outliers.

## **(2)**

In a situation with 10 variables and 4 mixture components, what is the number of free parameters for

**(a)**    a "VVV" Gaussian mixture model assuming fully flexible covariance matrices

**(b)** a "VII" Gaussian mixture model assuming spherical covariance matrices with potentially differing volumes

**(c)**  an "EEE" Gaussian mixture model with a flexible covariance matrix assumed equal in all mixture components

**(d)** a fully flexible skew-normal mixture

**(e)** a fully flexible mixture of multivariate t distributions

**(f)** a fully flexible mixture of skew-t distributions 

**(g)** a mixture of skew-t distributions where skewness parameters, degrees of freedom and Σ-matrices are assumed equal in all mixture components?


## **(4)**

A very simple strategy for applying a mixture to a big data set is the following:

**Step** 1 Draw a random data subset of ns observations

**Step 2** Compute the mixture ML estimators using the EM-algorithm on that subset.

**Step 3** Extend the fitted model to all other observations by computing the estimated posterior probabilities pik for observation xi to belong to cluster k (these
can be used as usual by maximisation to find a clustering).

For Gaussian mixtures, this can be done using Mclust in the following way:

**Step 1** as above.
**Step 2** Run `Mclust` on the subset.
**Step 3** Use function `predict.mclust` to extend the fitted model to all observations

Implement the big data method explained above, and apply it to the stars5000
data from Exercise 1 above. Use ns = 1000 (for a data set of size 5000, using
ns = 2000 will not win that much time). Take the time (this can be done using
the `system.time` command in R). Also take the time for running `Mclust` on those
data. Compare the running times and the results of the big data method and of
standard `Mclust` (it would be a good result for the big data method if results are
very similar, but the run time is much faster).
You may want to (but don't have to) write a single function for the big data method
that takes the data as input and optionally some other parameters such as ns or
`Mclust`-input parameters, so that you can use it easily in other situations.

