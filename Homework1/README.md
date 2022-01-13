ㅤ

#### **a)**ㅤ *Explain in your own words how trimmed k-means works:*

With trimmed k-means, part of the data is discarded prior to the
application of the (generalized) clustering algorithm, in a way that is
self-determined by the data structure rather then being somehow
arbitrary. Here, arbitrariness is intended in the direction or zones for
removing data, which would require the researcher to do some difficult
choices based on potentially large and complex multivariate samples.
Instead, trimmed k-means only needs to specify an hyperparameter
*α* ∈ \[0;1\] representing the fraction of observations to be ignored in
any case. Such trimming (known as *“impartial trimming”*) is integrated
in the optimization problem:

$\\underset{Y}{\\min} \\underset{t\_{1},...,t\_{k}}{\\min} \\, \\, \\,\\sum\_{X\_{i}\\in Y} \\Phi \\Big (\\underset{1\\leq j \\leq k}{\\inf} \|\|X\_{i} - t\_{j}\|\|\\Big)$

since this searches through the set of trimmed subsets of *X* with
 ≈ \[*n*(1−*α*)\] data points, and leads to a (generalized) k-means
clustering of the subsample with the lowest possible
*Φ* − *p**e**n**a**l**i**z**e**d* variation. In this way both outliers
and noisy observations are handled automatically by the robustified
algorithm, without simply removing extreme observations from the ends of
the feature space (it is also capable of detecting noisy points between
clusters).

The generalized k-means algorithm produces a partition of the
observations in the optimal regions (nontrimmed zone) into a specified
number of groups by minimizing the sum of the penalized euclidean
distances of each *X*<sub>*i*</sub> from its respective cluster center
*t*<sub>*j*</sub>; such optimal centroids no longer correspond, in
general, to group means (as *Φ*(*x*) = *x*<sup>2</sup> guaranteed in the
classical method) but these are rather *Φ* − *m**e**a**n**s* of the
corresponding cluster. The characterization of the optimal region
requires a scale estimate/optimal radius *S*.

ㅤ

#### **b)**ㅤ *How is the definition of k-medoids on the first page related to PAM as introduced in class? What is different, what is the same?*

K-medoids, as defined in the article, is a partitioning clustering
method for euclidean data based on the minimization of the sum of the
(unsquared) euclidean distances of each *X*<sub>*i*</sub> from its
respective cluster center *m*<sub>*c*(*i*)</sub>:

$\\underset{m\_{1},...,m\_{k}}{\\min} \\, \\, \\,\\sum\_{i=1}^n \\underset{1\\leq j \\leq k}{\\inf} \|\|X\_{i} - m\_{j}\|\|$

Such problem is **hard** to solve exactly, and requires local
optimization algorithms which typically only find local minimas of the
objective function. In the article this isn’t an issue since the
proposed example is kind of trivial.

PAM is the most popular algorithm for applying k-medoids (as intended by
Kaufman and Rousseeuw 1990) but it’s dissimilarity based, in the sense
that it just requires any generic dissimilarity matrix as input. Its
objective function is now the sum of pairwise general dissimilarities:

$T(C, m\_{1},...,m\_{k})=\\sum\_{i=1}^n \\, \\, \\,d(X\_{i},m\_{c(i)})$

and is minimized by swapping all the nonmedoid points and medoids
iteratively until convergence. What these methods have, apparently, in
common is the definition of medoid, defined as the data object of the
cluster for which the average dissimilarity to all the objects of the
group is minimal. This is actually employed in PAM but not in the method
discussed here, whose objective function is instead minimized by
geometric medians which are not necessarily points from the original
dataset. In the 1-dimensional case, as in the article’s example, the
geometric median coincides with the median which is again the most
central data object from the sample.

ㅤ

#### **c)**ㅤ *Two desirable robustness characteristics are a bounded influence function and a large breakdown point. Which of these does k-medoids fulfill?*

*K*-medoids method corresponds to a generalized *k*-means with *Φ* = *x*
and bounded *Ψ* = 1,   *ψ* = *s**i**g**n*(*x*). For each *k* medoid
*M*<sub>*i*</sub>   (*i*=1,...,*K*) it has bounded IF, since these are
proportional to *ψ*(*x*−*M*<sub>*i*</sub>)   (*i*=1,...,*K*) on each
cluster *A*<sub>*i*</sub> and decline elsewhere (at points
*x* ∉ *A*<sub>*i*</sub>). Apparently this property could lead us to
consider *K*-medoids as somehow robust, but the IFs just describe the
local effect on estimators given by an infinitesimal contamination at
*x*. Another important tool that partners with the IF is the BP, which
measures in a global way how much contamination the estimator can handle
before *“breaking down”* (i.e. becoming susceptible of arbitrary
deviation from *M*<sub>*i*</sub>(*X*)). Generalized *k*-means does not
achieve a desirable high breakdown point, since it can be proved to be
1/*n* (and converges asymptotically to 0, i.e. breaks down with just one
outlier). This is true regardless of the sample *X* under analysis and
the selected penalty function *Φ*, so this also extends to *k*-medoids.

ㅤ

#### **d)**ㅤ *Explain how and why trimmed k-means in terms of robustness properties is better than k-medoids, and how and why k-means is worse.*

Both *k*-medoids and trimmed *k*-means methods have bounded IF’s, but
this lasts show signs of greater robustness. Infact the IF’s of each
cluster’s center estimator *T*<sub>*i*</sub> are constant on the trimmed
region and again proportional to
*ψ*(*x*−*T*<sub>*i*</sub>)   (*i*=1,...,*K*) on each cluster
*A*<sub>*i*</sub>, but with an effect caused by a contamination at *x*
that is greater on the estimator for the centroid of the cluster
*T*<sub>*i*</sub> to which *x* belongs, and falling significantly for
the others *T*<sub>*j*</sub>, *j* ≠ *i*. This behaviour gets stronger as
the clusters are more separated
(*I**F*(*x*∉*A*<sub>*i*</sub> ; *T*<sub>*i*</sub>,*F*) ≈ 0)), up to the
ideal data structure of well-separated clusters with internal symmetry,
where the centroids estimators are comparable to independent location
*M* estimators applied to the subsamples coming from each distribution
*F*<sub>*i*</sub>.

We should worry less about the choice of *Φ* in trimmed k-means but this
has instead some implications in the untrimmed method: I already
mentioned that bounded *Ψ* is required for bounded IF’s for this (a
clear restriction), but inferential properties may also be compromised.
The population *k*-means (required for consistency and asymptotic
normality) may not exist if the necessary condition
∫ *Φ*(\|*x*\|)  *d**F*(*x*) \< ∞ fails, and this happens particularly
for heavy-tailed distributions that are responsible for outliers in
data. On the contrary, milder regularity conditions are requested for
making inference with the robuster method.

*k*-means (and *k*-medoids) main shortcoming is its low BP, which also
approaches zero asymptotically. For trimmed *k*-means, BP depends
heavily on the data sample structure: if we have well-separated clusters
with moderate contamination and stability, we can get to a BP that is at
most equal to the proportion *α* of discarded observations. The choice
of *α* is crucial as the exclusion of outliers can be missed, leading to
the same issues of the generalized k-means. However, in certain unstable
settings we can still imagine to break down the procedure by taking a
small group of untrimmed observations to infinity; this would create a
new cluster forcing some other two to be joined into only one. Under
these considerations we can also say:
*ε*<sub>*n*</sub><sup>\*</sup> ≤ inf {*α*, *m**i**n*(*n*<sub>*i*</sub>)/*n*}.

Let’s say that trimmed k-means shows its obvious limits in difficult
situations of high contamination and cluster instability, but its
robustness seems to be high when we have a set of well-shaped spherical
clusters. In general, it may be a useful alternative to non-robust
techniques (which break down way more easily).

ㅤ

#### **e)**ㅤ *What limitations do you see in the robustness of trimmed k-means and its treatment in this paper? Are there desirable properties that it does not have, or are not treated?*

The issue I *humbly* see is that it leaves to an optimization process
the decision about what an anomalous observation is. For sure it is an
advantage in many situations but I think it may have some drawbacks: in
general, it is difficult to decide whether a group of scattered
observations should be considered as noise or as an additional small
cluster. The researcher doesn’t have complete control on the impartial
trimming and tuning the method wrongly, by the choice of *α* and *k*,
may prevent him from retrieving the actual group structure. It may be
possible to fit the outlying observations by mixtures of skew and
heavy-tailed distributions which assign groups’ ownership probabilities
rather then discarding them completely, also avoiding the “dilemma”
mentioned earlier. In the article, I would have expected greater
emphasis on the choice of *α* and the effects on both bias and
efficiency, or whether it is best to overtrim when we’re not sure about
the actual proportion of contaminations. Also, given its peculiarity,
I’m not sure its results are comparable at all with clusterings obtained
with other methods.
