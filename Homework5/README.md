## **(1)**

On Virtuale under "Data sets" you'll find the data set wdbc.data. This
data set is taken from the [UCI Machine Learning Repository](http://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Diagnostic)).
Data are given about 569 breast cancer patients, and there are the two "true" classes
of benign (357 cases) and malignant (212 cases) tumors. There are ten quantitative
features in the dataset, the documentation of which is in the file wdbc.names, which
I also put on the IOL/Moodle site, however you don't really need to know what's
in there to do this exercise. Actually the features have been recorded for every cell
nucleus in an image and three different statistics of these are in the original dataset
so that the number of variables is 30, but I'm asking you to use only the first ten
variables, otherwise clustering can be very unstable and computationally hard.

Compute different clusterings of the data (use at least two different approaches
including a Gaussian mixture model and try out numbers of clusters up to 10) and
compare them first without using the information about benign vs. malign cancers in
the diagnosis variable wdbcdiag (this also means you should ignore the information
that there are 2 classes). Which clustering do you think is best?
Only after you have made a decision about your favourite clustering, use the ARI
to compare all these clusterings to wdbcdiag.
Discuss how it was possible, without using wdbcdiag, to recognise from the data
whether a clustering would be similar or less similar to the "true" clustering in
wdbcdiag. Consider the clustering that produced the best ARI value with wdbcdiag.
Are there reasons how one could have realised already without access to wdbcdiag
that this is a good and potentially the best clustering? Note also that there may be other legitimate and meaningful clusters in the data
about which we don't have information, so a clustering that has low ARI with
wdbcdiag isn't necessarily bad.

## **(3)**

On Moodle you can find the article 
> "Regularized k-means clustering of high-dimensional data and its asymptotic consistency" by Wei Sun and Junhui Wang, Electronic Journal of Statistics 6 (2012), 148-167.

  * Explain and motivate in your own words how regularized k-means differs from
k-means, and how regularized model-based clustering differs from model-based
clustering (Gaussian mixture models).
  * The authors state that the X-variables should be "centralized". Why is this
important? Do you think it would also be useful to standardize them to unit
variance? Why, or why not, or under what circumstances?
  * As opposed to the use of the Lasso in regression, the tuning constant λ here can-
not be chosen by optimizing a prediction error estimated by cross-validation.
Why not? What do the authors propose instead to choose the λ? Do you think
this could also be done involving the Rand or adjusted Rand index? Why or
why not?
  * What are advantages and good properties of the method according to the
authors? Do you think this is convincing? Can you think of any potential
disadvantages of the method?
