# Kernel functions

## Purpose of this document

When first knowing about the `kerntools` package, or perhaps after
downloading it, you may have asked yourself some of these questions:

- What are kernel functions and matrices?
- Why are they useful?
- How to compute them using `kerntools`?

This document is meant to answer exactly that. The first section,
[Introduction](#introduction), deals with the basic concepts regarding
the kernel approach. Subsequent sections review all kernels implemented
in `kerntools` one by one, highlighting their most important traits and
giving some examples of use.

## Introduction

### Informal definition

Kernel methods are all those machine learning methods that use kernel
functions. Let $X = {x_{1},x_{2},...,x_{N}}$ be a dataset of *N* objects
of any type (real-valued, graphs, rankings, matrices, categorical data,
etc.). The kernel function evaluates all possible pairs of objects,
which are then stored in a kernel matrix with dimension *NxN*:

$$\mathbf{K}_{\mathbf{i}\mathbf{j}} = k\left( x_{i},x_{j} \right)$$

Intuitively, we can understand kernels as functions
$k\left( x_{i},x_{j} \right)$ that measure the similarity between
$x_{i}$ and $x_{j}$. Generally speaking, all kernels are similarity
measures; however, not every similarity measure is a kernel: only those
that generate proper kernel matrices. Kernel matrices are characterized
by four properties: they are real-valued, squared, symmetric and PSD
(positive semi-definite: all their eigenvalues are nonnegative). Thus,
the kernel matrix **K**, and not the original dataset *X*, is the input
data for our kernel method. In this approach, choosing a good kernel
function is capital, because it’s the “lens” that the method will use to
“see” the data.

### Advantages

If using kernels implies an “extra step” (computing **K**), why should
we bother at all? Surely, it’s easier to use any other machine learning
method and dealing directly with *X*, isn’t it? Well, the thing is that
the kernel approach has some advantages:

- Most machine learning methods require a numeric input (some of them,
  like Random Forests, also can accept categorical data). Usually, if
  our dataset is different from the typical real-valued matrix **X**
  with *D* features (variables), we will have to recode it somehow. But
  in kernel methods, the objects in *X* are never represented explicitly
  (only via their pairwise similarities), which has an interesting
  consequence: kernels can handle data types that are not standard.
  Instead, the challenge is using (or designing) kernels that extract
  valuable information from specific data types or problems. This
  requires a notion of what is considered “similar” in that given
  context.

- Kernel functions can be combined or “custom-tailored” for a specific
  problem. For instance, previous knowledge about the problem can be
  encoded explicitly into the kernel. This is facilitated by the fact
  that kernel design is modular in nature: kernels can be combined. For
  example, the linear combination of kernel matrices is another kernel
  matrix.

- When using nonlinear kernel functions, methods that were originally
  linear (for instance: Principal Components Analysis or PCA) “become”
  nonlinear. In some sense, we “trick” the method: it *still* works in a
  linear way, but as its input data has been transformed by a nonlinear
  function, the output it produces is nonlinear.

- Curse of dimensionality: contrarily to a lot of other methods, this is
  not a big issue when working with kernels. In fact, when dealing with
  *NxD* datasets where *N \<\< D*, kernels can decrease the
  computational burden if used cleverly, as the method will deal with a
  kernel matrix of dimension *NxN*.

### Feature space

The most basic kernel function is the dot product of two real vectors
(see [Linear kernel](#linear-kernel)). This definition is *very*
important, because all kernels can be written as inner products in some
space:

$$\left. x\rightarrow\varphi(x) \right.$$$$k\left( x_{i},x_{j} \right) = < \varphi\left( x_{i} \right),\varphi\left( x_{j} \right) >$$$\varphi(x)$
maps an object from the original (input) space to the so-called *feature
space*. This feature space typically has a higher dimension than the
original space (other times it has the same dimension) and it’s endowed
with an inner product. However, all the process happens “under the
hood”: the kernel function uses this mapping without telling so (recall
that the objects’ explicit representations are never used). Sometimes,
even knowing the feature space associated to a given kernel can be
tricky. However, this take avoids the computational cost of 1)
explicitly mapping the dataset into feature space (which sometimes can
be *very* high dimensional) + 2) computing the inner product there: the
kernel conflates the two in a single step, which is known as the *kernel
trick*. By the way, this also allows solving nonlinear problems using
linear methods: the kernel takes the original data onto a higher
dimensional space where it can be operated in a linear way with inner
products. (A visual example of this can be found
[here](https://es.m.wikipedia.org/wiki/Archivo:Kernel_trick_idea.svg)).

## Kernels for real vectors

Now we introduce the kernels in `kerntools` that deal with real vectors.
There are *a lot* of kernels that can handle this data, but here we are
going to focus only on those implemented in this package. By real
vectors I refer the typical continuous dataset **X** with *N* samples in
rows and *D* continuous variables (or features) in columns. Thus, for
each sample we have \$\mathbf{x} \in {\rm I\\R}^D\$.

### Linear kernel

#### Definition

The Linear kernel is the dot product of two real vectors:

$$Linear\left( \mathbf{x}_{\mathbf{i}},\mathbf{x}_{\mathbf{j}} \right) = {\mathbf{x}_{\mathbf{i}}}^{T}\mathbf{x}_{\mathbf{j}}$$

It can return positive values, negative values, and even be zero.
Obviously, it has a linear nature. Thus, if we use it in a “kernelized”
machine learning method, the output will be the same that when using the
non-kernelized method (for instance: linear kPCA is equivalent to
standard PCA).

The linear kernel is related to the cosine of the angle between two
vectors, which is maximum when they share the same direction and 0 if
they are perpendicular. Sometimes, is more easy to see the linear kernel
as a “similarity score” if the cosine normalization is applied:

$$Linear_{cosine}\left( \mathbf{x}_{\mathbf{i}},\mathbf{x}_{\mathbf{j}} \right) = \frac{{\mathbf{x}_{\mathbf{i}}}^{T}\mathbf{x}_{\mathbf{j}}}{\sqrt{{\mathbf{x}_{\mathbf{i}}}^{T}\mathbf{x}_{\mathbf{i}}}\sqrt{{\mathbf{x}_{\mathbf{j}}}^{T}\mathbf{x}_{\mathbf{j}}}}$$
This is also known as the “cosine similarity”. A score of 1 implies that
$\mathbf{x}_{\mathbf{i}},\mathbf{x}_{\mathbf{j}}$ are equal, while the
lower the score, the lower the similarity.

#### Feature space

This is the only kernel that fulfills that feature space = input space.

#### Usage

This kernel is suitable for typical numeric, continuous datasets
(especially if you know in advance that your problem has a linear
solution). For instance, if we work with the famous iris dataset:

``` r
Linear(iris[,1:4])
```

In general, however, it is recommended to standardize a dataset before
presenting it to the linear kernel. Imagine that your features are very
different in magnitude: for instance, if the range of feature 1 is
$\lbrack 0.01,0.1\rbrack$ and the range of feature 2 is
$\left\lbrack 10^{2},10^{3} \right\rbrack$, feature 2 will absolutely
dominate the output of the kernel, even if this feature is almost
irrelevant to the problem you want to solve. The reverse case arises
when you want to *manually* prioritize some features over others (for
example, because you have previous knowledge about your problem). In
that case, you can weight them in advance, specifying your desired
coefficients. For example:

``` r
Linear(iris[,1:4], coeff = c(1/2, 1/4, 1/8, 1/8))
```

Cosine-normalization of the linear kernel can be done with:

``` r
Linear(iris[,1:4], cos.norm = TRUE)
```

### RBF kernel

#### Definition

The Gaussian Radial Basis Function (RBF) kernel is regarded as the gold
standard among kernels. Typically is written as:

$$RBF\left( \mathbf{x}_{\mathbf{i}},\mathbf{x}_{\mathbf{j}} \right) = e^{- \gamma{{||}\mathbf{x}_{\mathbf{i}} - \mathbf{x}_{\mathbf{j}}{||}}^{2}}$$
It is immediately obvious that it is nonlinear. The output ranges
between 0 and 1; again, the higher the RBF kernel value, the higher the
similarity between $\mathbf{x}_{\mathbf{i}},\mathbf{x}_{\mathbf{j}}$.
The ${||}\mathbf{x}_{\mathbf{i}} - \mathbf{x}_{\mathbf{j}}{||}$ term
correspond to the Euclidean distance between two real vectors, while
$\gamma > 0$ is a hyperparameter. There is an alternative definition of
the RBF in some books or packages, where it has a “bandwidth”
hyperparameter called $\sigma > 0$. Don’t worry. The relationship
between the two is: $\gamma = \frac{1}{2\sigma^{2}}$.

The RBF kernel is related to the linear kernel. This is because the
inner product induces a L2 norm that, in real vectors, is equivalent to
the Euclidean distance. But contrarily to the linear kernel, the RBF
kernel is an universal approximator: it can be used to model any
decision boundary. In this regard, $\gamma$ governs the “smoothness” of
this kernel when adapting to data. The smaller the value, the more
complex the resulting model, and the risk of overfitting is greater.
Instead, if the value is too small, the decision boundary become very
smooth and “linearized”, to the point that (with very large values) the
RBF kernel starts behaving as the linear kernel.

#### Feature space

With this kernel, data is mapped onto an infinite-dimensional feature
space. Several techniques to approximate it exist. For instance,
$\varphi(\mathbf{x})$ can be computed via Taylor series expansion (more
info in [Explicit Approximations of the Gaussian
Kernel](https://arxiv.org/abs/1109.4603)).

#### Usage

This kernel is perfectly fine for typical numeric, continuous datasets,
especially if the solution to your problem is nonlinear. Here, you can
see how to compute it with `kerntools` using the iris example dataset:

``` r
RBF(iris[,1:4], g= 0.1)
```

(Centering or not the dataset is irrelevant, but it is advisable to
scale the variables by the their standard deviation before presenting
the dataset to the RBF kernel. See [Linear kernel](#linear-kernel) for
more info.)

Here, `RBF(...,g)` is the value of the hyperparameter gamma, which
should be chosen by the user. There are some heuristics that may be
useful for choosing a good value. In that respect, `kerntools` provides
this handy function:

``` r
estimate_gamma(iris[,1:4])
```

(Remember that, in `kerntools`, the hyperparameter of the RBF kernel is
$\gamma$, not $\sigma$.)

Finally, `kerntools`’s implementation currently does not compute the
feature space of the RBF kernel.

### Laplacian kernel

#### Definition

The laplacian kernel is similar to the RBF kernel. However, instead of
computing the L2 norm (Euclidean distance) between two real vectors, it
uses the L1 norm:

$$Laplacian\left( \mathbf{x}_{\mathbf{i}},\mathbf{x}_{\mathbf{j}} \right) = e^{- \gamma{|\mathbf{x}_{\mathbf{i}} - \mathbf{x}_{\mathbf{j}}|}}$$

It ranges between 0 and 1. The higher the laplacian kernel value, the
higher the similarity between
$\mathbf{x}_{\mathbf{i}},\mathbf{x}_{\mathbf{j}}$. The hyperparameter
$\gamma > 0$ has a similar role that in the RBF kernel. However, now its
relationship to $\sigma$ is: $\gamma = \frac{1}{\sigma}$. The
hyperparameter governs how the laplacian kernel adapts to data, though
using the L1 norm makes this kernel more “sharp” than the RBF kernel.

#### Usage

You can use this kernel with a typical numeric, continuous dataset,
especially if the solution to your problem is nonlinear. An example of
use with the famous iris dataset:

``` r
Laplace(iris[,1:4], g= 0.1)
```

Here, `Laplace(...,g)` is the value of the hyperparameter $\gamma$,
which should be chosen by the user. Centering or not the dataset is
irrelevant, but it is advisable to scale the variables by the their
standard deviation before presenting it to the laplacian kernel. See
[Linear kernel](#linear-kernel) for more info.

Remember that in `kerntools`, the hyperparameter of the laplacian kernel
is $\gamma$, not $\sigma$. Moreover, `kerntools`’s implementation right
now does not compute the feature space of the laplace kernel.

## Kernels for real matrices

The following kernel works with data such as each sample is a real
matrix (its entries are real numbers).

### Frobenius kernel

#### Definition

The Frobenius kernel of two matrices
$\mathbf{X}_{\mathbf{i}},\mathbf{X}_{\mathbf{j}}$ with dimension *NxD*
is defined as:

$$Frobenius\left( \mathbf{X}_{\mathbf{i}},\mathbf{X}_{\mathbf{j}} \right) = tr\left( \mathbf{X}_{\mathbf{i}}{\mathbf{X}_{\mathbf{j}}}^{T} \right) = \sum\limits_{n = 1}^{N}\sum\limits_{d = 1}^{D}\mathbf{X}_{ind} \circ \mathbf{X}_{jnd}$$$\circ$
denotes the Hadamard product. It is an inner product and, thus, the
Frobenius kernel is a sort of “linear kernel” for matrices. The result
can be positive, negative or even zero. As always, the higher the kernel
the higher the similarity between the two samples. This is understood
more easily when the cosine normalization is applied:

$$Frobenius_{cosine}\left( \mathbf{X}_{\mathbf{i}},\mathbf{X}_{\mathbf{j}} \right) = \frac{tr\left( \mathbf{X}_{\mathbf{i}}{\mathbf{X}_{\mathbf{j}}}^{T} \right)}{\sqrt{tr\left( \mathbf{X}_{\mathbf{i}}{\mathbf{X}_{\mathbf{i}}}^{T} \right)}\sqrt{tr\left( \mathbf{X}_{\mathbf{j}}{\mathbf{X}_{\mathbf{j}}}^{T} \right)}}$$

A score of 1 implies that the two samples (matrices) are equal.

#### Feature space

The feature space of this kernel is straightforward. It can be obtained
by vectorizing (converting into column vectors) each one of the
matrices.

#### Usage

This kernel can be used to compare matrix-formatted information, for
instance image data. It can be also used to compare the similarity
between two kernel matrices. In `kerntools`, the Frobenius kernel can be
called typing:

``` r
### Let's suppose our samples come from three different matrix tables:
setosa <- iris[iris$Species == "setosa", 1:4]
versicolor <- iris[iris$Species == "versicolor", 1:4]
virginica <- iris[iris$Species == "virginica", 1:4]
matrices <- list(setosa,versicolor,virginica)
Frobenius(matrices)
```

Cosine-normalization of the Frobenius kernel can be done with:

``` r
Frobenius(matrices, cos.norm = TRUE)
```

and the feature space (cosine-normalized or not) can be recovered like
this:

``` r
Frobenius(matrices, feat_space = TRUE)
```

## Kernels for abundances (counts or frequencies)

The following kernels deal with counts and frequency data, be it
absolute or relative frequencies. This kind of data lives in datasets of
dimension *NxD*, which (at first glance) may appear very similar to real
data. However, some nuances should be taken into account. Count data
consists in nonnegative integer vectors of dimension *D* \$\mathbf{x}
\in {\rm I\\R}\_{\>0}^D\$ that quantify the abundance of some trait or
analyzed variable. This data has a non-symmetric distribution, a
potentially large dynamic range (from zero up to millions) and, when the
quantification is done by an instrument, other factors like upper and
lower bounds on detection can be at play. Thus, zeroes may be ambiguous:
they may signal an actual absence, or a quantity below the detection
limit.

Sometimes, counts or absolute frequencies are converted to relative
frequencies. This data is usually called *compositional* and constitutes
a specific category within abundance data. In that case, all samples
(rows) sum to an arbitrary number (for instance 1, or 100). This data
does not live in the positive real space \$\mathbf{x} \in {\rm
I\\R}\_{\>0}^D\$ but in the simplex, which is a *D - 1* subspace of it.
This is important because the variables are *not* independent, but parts
of a whole. Thus, increasing the abundance of a feature within a sample
decreases the abundance of the rest.

The next two kernels are well known in some fields, like ecology
studies, because of their relation to widespread diversity metrics: the
Bray-Curtis dissimilarity and the Ruzicka distance. The latter two (the
compositional linear and Aitchison kernels) are kernels specific for
compositional data.

### Bray-Curtis kernel

#### Definition

The Bray-Curtis kernel is defined as:

$$BCK\left( x_{i},x_{j} \right) = 2\frac{\sum\limits_{d}min\left( x_{id},x_{jd} \right)}{\sum\limits_{d}\left( x_{id} + x_{jd} \right)}$$
where $x_{id},x_{jd}$ are two samples (that is: two rows of the
dataset). It ranges between 0 and 1.

#### Usage

This kernel is $BCK = 1 - BCD$, where BCD is the Bray-Curtis
dissimilarity (it is called “dissimilarity” because it is semimetric
and, therefore, not a distance in the strict sense). It can be applied
over counts and frequencies; however, the “correctness” of using
directly this kernel over compositional data is disputed.

In `kerntools`, it can be computed typing:

``` r
### Our example dataset contains the bacterial abundance of *D* species in *N* soil samples.
data <- soil$abund
BrayCurtis(data)
```

### Ruzicka kernel

#### Definition

The Ruzicka kernel (also known as: quantitative Jaccard, or min-max
kernel) is defined as:

$$Ruzicka\left( x_{i},x_{j} \right) = \frac{\sum\limits_{d}min\left( x_{id},x_{jd} \right)}{\sum\limits_{d}max\left( x_{id} + x_{jd} \right)}$$
where $x_{id},x_{jd}$ are two samples (that is: two rows of the
dataset). It is bounded between 0 and 1.

#### Usage

The Ruzicka kernel is related to the quantitative Jaccard distance, aka
Ruzicka distance, in the same manner as the Bray-Curtis kernel is
related to the Bray-Curtis dissimilarity. The Ruzicka distance is
rank-order similar to BCD, but it fulfills the conditions to be a
metric. It is also related to the Jaccard kernel (which is introduced in
the [Jaccard kernel](#jaccard-kernel) subsection) for sets or
presence/absence data.

It can be applied over counts and frequencies; however, the
“correctness” of using directly this kernel over compositional data is
disputed.

In `kerntools`, it can be computed typing:

``` r
### Our example dataset contains the bacterial abundance of *D* species in *N* soil samples.
data <- soil$abund
Ruzicka(data)
```

### Compositional-linear kernel

#### Definition

We define the compositional-linear kernel as:

$$compLinear\left( \mathbf{x}_{\mathbf{i}},\mathbf{x}_{\mathbf{j}} \right) = \sum\limits_{d = 1}^{D}log\left( \frac{x_{id}}{geo\left( x_{i} \right)} \right) \cdot log\left( \frac{x_{jd}}{geo\left( x_{j} \right)} \right)$$
It is analogous to the linear kernel and it can return positive values,
negative values, and even be zero. We can apply the cosine
normalization, so a score of 1 implies that
$\mathbf{x}_{\mathbf{i}},\mathbf{x}_{\mathbf{j}}$ are equal, while lower
values imply a lower similarity between
$\mathbf{x}_{\mathbf{i}},\mathbf{x}_{\mathbf{j}}$ (see [Linear
kernel](#linear-kernel)).

#### Feature space

The feature space for the compositional-linear kernel is the centered
log-ratio (clr) transformation of the original data.

#### Usage

Unlike the two previous kernels, doing:

``` r
### Our example dataset contains the bacterial abundance of *D* species in *N* soil samples.
data <- soil$abund
cLinear(data)
```

will return an error. This is because compositional kernels require
samples that are strictly positive. This is a problem, because in
real-life, a lot of compositional datasets have zeroes. These zeroes can
have different origins: for example, they may be structural, or they may
be caused by a signal so small that falls below an instrument’s
detection limit, etc. In any case, a zero-replacement strategy should be
applied before computing this kernel. (Other possibilities include
splitting the dataset, if we are sure that the zeroes are structural, or
to sum columns, which is called “amalgamation” in Compositional Data
Analysis. This is not always possible. In the case of the `soil`
dataset, if we wanted to follow the latter route, we could sum all
bacteria of the same species (or genus, or family), aiming for higher
taxonomic ranks, but losing resolution in the process.) The two
compositional kernels present in `kerntools` support a very simple
built-in replacement strategy: replacing all zeroes by a pseudocount
(i.e. a small value):

``` r
cLinear(soil$abund,zeros = "pseudo")
```

(For more sophisticated approaches you may use other R packages like
‘zCompositions’.)

The feature space can be returned typing:

``` r
cLinear(soil$abund,zeros = "pseudo",feat_space=TRUE)
```

And the cosine-normalization of this kernel can be done with:

``` r
cLinear(soil$abund,zeros = "pseudo", cos.norm = TRUE)
```

If both `cos.norm` and `feat_space` are TRUE, the feature space will be
scaled accordingly.

### Aitchison kernel

#### Definition

The Aitchison kernel is analogous to the RBF kernel but for
compositional data:

$$Aitchison\left( \mathbf{x}_{\mathbf{i}},\mathbf{x}_{\mathbf{j}} \right) = exp\left( - \gamma\sum\limits_{d = 1}^{D}{log\left( \frac{x_{id}}{geo\left( x_{i} \right)} \right) - log\left( \frac{x_{jd}}{geo\left( x_{j} \right)} \right)}^{2} \right)$$
where $\mathbf{x}_{\mathbf{i}},\mathbf{x}_{\mathbf{j}}$ represent the
abundances in two different samples, $geo(.)$ denotes the geometric
mean, and $\gamma > 0$ is a hyperparameter that has to be tuned. This
non-linear kernel derives from the Aitchison distance, which is
Euclidean, and thus the previous definition results in a valid kernel.
The output of this kernel ranges between 0 and 1; the higher the
Aitchison kernel value, the higher the similarity between
$\mathbf{x}_{\mathbf{i}},\mathbf{x}_{\mathbf{j}}$. The logarithm terms
can be identified as the compositional clr-transformation of the
original data. The relationship between this kernel and the
compositional-linear kernel is the same than between the RBF and the
linear kernels.

#### Feature space

See: [RBF kernel](#rbf-kernel).

#### Usage

Again, doing:

``` r
### Our example dataset contains the bacterial abundance of *D* species in *N* soil samples.
Aitchison(soil$abund,g=0.1)
```

will return an error, for the very same reasons stated in the “Usage”
subsection of [Compositional-linear
kernel](#compositional-linear-kernel). A possible solution is replacing
all zeroes by a pseudocount (i.e. a small value):

``` r
### Our example dataset contains the bacterial abundance of *D* species in *N* soil samples.
Aitchison(soil$abund,zeros = "pseudo", g=0.1))
```

As always, `Aitchison(...,g)` is the value of the hyperparameter gamma,
which should be chosen by the user. The same heuristics of the RBF
kernel may be useful for choosing a good gamma.

Finally, remember that `kerntools`’s implementation does not support
computing the feature space of the Aitchison kernel.

## Kernels for categorical data

Categorical data consist on qualitative variables that can take a finite
number of values (called classes, modalities, categories, or levels).
Within this package, “categorical data” refers to nominal data whose
categories don’t have a meaningful order, in contrast with ordinal data.
The typical categorical dataset is a matrix *X* with *N* samples and *D*
categorical variables or features. (In R, “data.frames” having columns
of class “factor” can be used for this kind of data).

### Dirac kernel

#### Definition

First we define the Overlap kernel over a single categorical variable
as:

$$\left. x_{id} = x_{jd}\Leftrightarrow Overlap\left( x_{id},x_{jd} \right) = 1 \right.$$$$\left. x_{id} \neq x_{jd}\Leftrightarrow Overlap\left( x_{id},x_{jd} \right) = 0 \right.$$

where $x_{id},x_{jd}$ are two samples and *D = 1*. This is the most
basic categorical kernel, akin in a way to the linear kernel for
continuous data. From there, we can define the Dirac kernel over
multivariate categorical data as the linear combination of the Overlap
kernel’s evaluations:

$$Dirac\left( x_{i},x_{j} \right) = \sum\limits_{d}c_{d}Overlap\left( x_{id},x_{jd} \right)$$
The coefficients $c_{d} > 0$ can be all 1 (and then, the Dirac kernel is
simply the sum of the Overlap kernel evaluations), *1/D* (so the Dirac
kernel is the mean), or may be used to weight differently the variables
(for example, because you have previous knowledge about your problem).
When the coefficients are all set to *1/D*, the values for the Dirac
kernel range between 0 and 1, the 1 meaning that the two samples are
equal, and 0 that they are completely different. This is equivalent to
computing the “sum” Dirac kernel and then cosine-normalizing it. When
custom coefficients are used but $\sum_{d}c_{d} = 1$, this output of
this kernel is a weighted average bounded between 0 and 1 (thus, trying
to cosine-normalize it is useless).

#### Feature space

We can map our data onto the feature space of the Overlap kernel via
one-hot-encoding. That is: the map $\varphi(x)$ recodes a categorical
variable (using the R nomenclature, “factor”) with *L* levels into *L*
dummy variables (each one representing a given level) that can be 0
or 1. For instance, if we have a variable that can take three levels
(“A”, “B” and “C”):

``` r
cat_feat <- data.frame(var=factor(sample(LETTERS[1:3],10,replace = TRUE)))
rownames(cat_feat) <- 1:10
cat_feat
#>    var
#> 1    A
#> 2    C
#> 3    C
#> 4    A
#> 5    B
#> 6    C
#> 7    C
#> 8    C
#> 9    B
#> 10   A
```

the feature space for the Overlap kernel is:

``` r
dummy_data(cat_feat)
#>    var_A var_B var_C
#> 1      1     0     0
#> 2      0     0     1
#> 3      0     0     1
#> 4      1     0     0
#> 5      0     1     0
#> 6      0     0     1
#> 7      0     0     1
#> 8      0     0     1
#> 9      0     1     0
#> 10     1     0     0
```

The feature space for the Dirac kernel results of “concatenating” the
feature spaces for all *D* categorical variables. If the weights
$c_{d} \neq 1$, then data in feature space should be scaled accordingly.
For instance, if $c_{d} = 1/D$ for all weights, which (as stated before)
is the same than cosine-normalize the Dirac kernel, then the feature
space is normalized to unit length.

#### Usage

In `kerntools`, we need a *NxD* dataset of class “matrix” (2-d
dimensional array) that contains only characters, or instead a dataset
of class “data.frame” (with all entries = “character”, or columns =
“factor”). `showdata` is an example dataset that meets this
requirements. The Dirac kernel can be called doing:

``` r
Dirac(showdata)
```

Moreover, if we want to compute the feature space too, we can type:

``` r
Dirac(showdata,feat_space = TRUE)
```

By default `Dirac(...,comp = "sum")`; that is, all coefficients $c_{d}$
are set to 1. If we want the Dirac kernel computed as the mean of the
Overlap kernel evaluations, we can do `Dirac(...,comp = "mean")`.
Instead, to have the average weighting according to some coefficients,
we can do:

``` r
coeffs <- c(1/8,1/4,1/4,1/4,1/8)
Dirac(showdata, comp = "weighted",  coeff=coeffs)
```

The feature space returned by this function will change accordingly to
the coefficients chosen. Of course, `Dirac(...,comp = "mean")` is the
same than doing:

``` r
coeffs <- rep(1/ncol(showdata),length(showdata))
Dirac(showdata, comp = "weighted",  coeff=coeffs)
```

## Kernels for sets

A set *S* is a collection of different objects (characters, numbers,
strings, … even other sets). We call *universe* to the set of all
objects we want to consider in a given problem. The objects of *S* are
called their *elements*, and we denote as $|S|$ the size of the set:
that is, the number of elements it contains.

### Intersect kernel

#### Definition

Imagine that we have a dataset (size *N*) such as each sample is a set:
$S_{1},...,S_{N}$. The most simple kernel for this data is the Intersect
kernel:

$$Intersect\left( S_{i},S_{j} \right) = \left| S_{i} \cap S_{j} \right|$$
that is, the size of the two sets’ shared elements. This is akin to the
linear kernel for continuous data (see the “Feature space” subsection
below). This kernel ranges from 0 (the sets are completely different) to
any positive number. In general higher values indicate higher
similarity; however, some sets may be bigger (that is: they contain more
elements) than others. You can use the cosine-normalization to remove
the size effect and have a result that is bounded between 0 and 1:

$$Intersect_{cosine}\left( S_{i},S_{j} \right) = \frac{\left| S_{i} \cap S_{j} \right|}{\sqrt{\left| S_{i} \right| \cdot \left| S_{j} \right|}}$$

From there, we can consider a multivariate case: each sample has several
(*D*) variables or features that happen to be sets. As we have done with
the Dirac kernel, we can do a linear combination of the Intersect
kernel’s evaluations to produce a “multivariate” kernel:

$$Intersect_{multi}\left( S_{i},S_{j} \right) = \sum\limits_{d = 1}^{D}c_{d}Intersect\left( S_{id},S_{jd} \right)$$

The coefficients $c_{d} > 0$ can be all 1 (and then, this multivariate
kernel is simply the sum of the Intersect kernel evaluations), *1/D* (so
the Intersect-multi kernel is the mean of evaluations), or may be used
to weight differently the variables (for example, because you have
previous knowledge about your problem). In this last case, the
Intersect-multi kernel is a lineal combination of the Intersect kernel’s
evaluations. When the coefficients are all set to *1/D* and we also
perform a cosine normalization over the univariate evaluations, the
values for the Intersect-multi kernel range between 0 and 1, the 1
meaning that the two samples are equal, and 0 that they are completely
different.

#### Feature space

We can map our data onto the feature space of the Intersect kernel via a
sort of one-hot-encoding. In this case, the map $\varphi(S)$ recodes a
set coming from a universe with *L* elements into *L* dummy variables
(each one representing a given level) that can be 1, if a given element
*l* is present in the set *S*, or 0 if it is absent. For example:

**Favorite colors of 3 people:**

``` r
universe <-  c("blue","green","lightblue","orange","purple","red","white","yellow")
person1 <- c(2, 6, 8)
person2 <- c(1, 8)
person3 <- c(2, 3,  6)
list(person1=universe[person1],person2=universe[person2],person3=universe[person3])
#> $person1
#> [1] "green"  "red"    "yellow"
#> 
#> $person2
#> [1] "blue"   "yellow"
#> 
#> $person3
#> [1] "green"     "lightblue" "red"
```

The data mapped onto the feature space for this kernel is:

``` r
colors <- matrix(0,nrow=3,ncol=length(universe))
rownames(colors) <- c("person1","person2","person3")
colnames(colors) <- universe
colors[1,person1] <- 1
colors[2,person2] <- 1
colors[3,person3] <- 1
colors
#>         blue green lightblue orange purple red white yellow
#> person1    0     1         0      0      0   1     0      1
#> person2    1     0         0      0      0   0     0      1
#> person3    0     1         1      0      0   1     0      0
```

and now the Intersect kernel can be computed directly as the dot product
between samples.

If we use the Intersect-cosine, then the feature space is normalized
accordingly to the unit vector:

``` r
cosnormX(colors)
#>              blue     green lightblue orange purple       red white    yellow
#> person1 0.0000000 0.5773503 0.0000000      0      0 0.5773503     0 0.5773503
#> person2 0.7071068 0.0000000 0.0000000      0      0 0.0000000     0 0.7071068
#> person3 0.0000000 0.5773503 0.5773503      0      0 0.5773503     0 0.0000000
```

Moreover, the Intersect-multi kernel’s feature space results of
“concatenating” the feature spaces of all sets associated to a sample,
each one scaled accordingly by the weights $c_{d}$.

#### Usage

`kerntools` implementation requires that the set members (elements) are
character symbols (length=1) combined into a string. Compare this with
the last example:

``` r
universe <- c("b","g","l","o","p","r","w","y")
colors <- matrix("",ncol=1,nrow=3)
rownames(colors) <-  c("person1","person2","person3")
colors[1,] <-paste(universe[person1],collapse="")
colors[2,] <-paste(universe[person2],collapse="")
colors[3,] <-paste(universe[person3],collapse="")
colors
#>         [,1] 
#> person1 "gry"
#> person2 "by" 
#> person3 "glr"
```

That is a *NxD* “matrix” or “data.frame” containing characters. Each
sample is in a different row. Now, you can do:

``` r
Intersect(colors,elements = universe)
```

This function requires that we state the universe’s elements explicitly.
The feature space can be returned too:

``` r
Intersect(colors, elements = universe, feat_space = TRUE)
```

If data is multivariate ($D > 1$ columns, and each one contains a set
feature), the universe should span all elements found in all the *D*
sets. For instance, consider this second feature:

**Best friends of the same 3 people:**

(*Anna = A*; *Bernard = B*; *Cesc = C*; *Diana = D* ; *Elsa = E* ; *Fran
= F*; *Gisele = G*; *Horaci = H*; *Iria = I*; *Jaume = J*)

``` r
universe2 <-  LETTERS[1:10]
person1 <- c(1,2,9,10)
person2 <- c(1,4,6,7)
person3 <- c(2,3,5,9,10)
person1=universe2[person1]
person2=universe2[person2]
person3=universe2[person3]
list(person1=person3,person2=person3,person3=person3)
#> $person1
#> [1] "B" "C" "E" "I" "J"
#> 
#> $person2
#> [1] "B" "C" "E" "I" "J"
#> 
#> $person3
#> [1] "B" "C" "E" "I" "J"
friends <- colors
friends[1,] <- paste(person1,collapse="")
friends[2,] <- paste(person2,collapse="")
friends[3,] <- paste(person3,collapse="")
friends
#>         [,1]   
#> person1 "ABIJ" 
#> person2 "ADFG" 
#> person3 "BCEIJ"
```

So the final dataset is:

``` r
set_data <- data.frame(colors=colors,friends=friends)
set_data
#>         colors friends
#> person1    gry    ABIJ
#> person2     by    ADFG
#> person3    glr   BCEIJ
```

And now we can compute the kernel:

``` r
Intersect(set_data,elements = c(universe,universe2),comp="sum")
```

(We should indicate the elements of the two universes: colors and
friends).

When `Intersect(...,comp = "sum")`, all $c_{d} = 1$ and, in principle,
all “features” (following with our last example: “friends” and “colors”)
have the same weight. In practice you should be careful, because
`kerntools` doesn’t cosine-normalizes each univariate evaluation of the
kernel. Thus, you may inadvertently introduce some bias in favor of
“bigger” sets (“friends” will have, implicitly, more weight than
“colors”). Instead, when `comp = "mean"` or `comp = "weighted"` is
chosen, `kerntools` cosine-normalizes all univariate evaluations, so the
result ranges between 0 and 1, and there are not implicit biases.
`Intersect(...,comp = "mean")` sets $c_{d} = 1/D$ and therefore gives
the same importance to all features. However, if you want to give more
importance to some features over others, you can choose manually the
weights that you want:

``` r
### We will make that "friends" has 3 times more importance than "colors":
coeffs <- c(1/4,3/4)
Intersect(set_data, elements = c(universe,universe2),comp = "weighted",  coeff=coeffs)
```

Here, the result of the kernel (between 0 and 1) is the weighted average
of the features. If you set `Intersect(..., feat_space = TRUE)`, the
feature space returned by this function will be scaled accordingly to
the $c_{d}$ coefficients chosen. The same will occur when doing
`Intersect(...,comp = "mean")`, but in this case the scaling consist in
simply dividing by *D* (the set features).

### Jaccard kernel

#### Definition

The Jaccard index (or Jaccard kernel) of two sets is defined as:

$$Jaccard\left( S_{i},S_{j} \right) = \frac{\left| S_{i} \cap S_{j} \right|}{\left| S_{i} \cup S_{j} \right|}$$
It naturally ranges between 0 (the sets don’t share any element) and 1
(the sets are identical). Thus, the cosine normalization is redundant
for this kernel.

Like the Dirac and the Intersect kernels, we can consider a multivariate
case for the Jaccard kernel. If each sample has several (*D*) variables
that happen to be sets, we can do a linear combination of the kernel’s
evaluations to produce a “multivariate” kernel:

$$Jaccard_{multi}\left( S_{i},S_{j} \right) = \sum\limits_{d = 1}^{D}c_{d}Jaccard\left( S_{id},S_{jd} \right)$$

The coefficients $c_{d} > 0$ can all be 1 (and then, this multivariate
kernel is simply the sum of the Jaccard kernel evaluations), *1/D* (so
the Jaccard-multi kernel is the mean of evaluations), or may be used to
weight variables (sets) differently (for example, because you have
previous knowledge about your problem). In this last case, the
Jaccard-multi kernel is a weighted average of the Jaccard kernel’s
evaluations. When the coefficients are all set to $\sum c_{d} = 1$, the
values given by the Jaccard-multi kernel range between 0 and 1.

#### Usage

The usage is pretty similar to the Intersect kernel:

``` r
Jaccard(set_data,elements = c(universe,universe2),comp="sum")
```

The only nuance here is that the univariate Jaccard kernel always is
bounded between 0 and 1. There are not underlying biases and all
features have the same weight (unless you explicitly do
`Intersect(...,comp = "weighted")` and enter some coefficients). Thus,
this:

``` r
cosNorm(Jaccard(set_data,elements = c(universe,universe2),comp="sum"))
```

or this:

``` r
D <- ncol(set_data)
Jaccard(set_data,elements = c(universe,universe2),comp="sum")/D
```

Is the same than this:

``` r
Jaccard(set_data,elements = c(universe,universe2),comp="mean")
```

Also, take into account that the
[`Jaccard()`](https://elies-ramon.github.io/kerntools/reference/Jaccard.md)
function does not return the data mapped onto the feature space of the
kernel.

## Kernels for ordinal data, ranks, and permutations

Ordinal data consist on qualitative variables that can take a finite
number of values that follow an intrinsic and meaningful order. In that
regard, ordinal data is in-between continuous and categorical (nominal)
data. For instance, in a survey, the typical agree/disagree questions
that have answers like “strongly agree”/“agree”/“neither agree or
disagree”/“disagree”/etc. are ordinal variables. This is related to the
concept of “ranking”, which consists in assigning the ordering labels
“1st”, “2nd”, “3rd”, … to these values.

### Kendall’s $\tau$ kernel

#### Definition

Kendall’s $\tau$ coefficient measures the similarity or rank correlation
between to two ordinal variables. It can be applied to different ordinal
variables or to different rankings of the same variable. It is also a
kernel function, defined as:

$$Kendall_{\tau}(x,y) = \frac{nc - nd}{np}$$ Let
$\left( x_{1},y_{1} \right),...,\left( x_{L},y_{L} \right)$ be two
ordinal variables, or two samples of the same variable. For simplicity,
we shall consider that all values taken by $x_{i}$ and $y_{i}$ are
unique. Observations $\left( x_{i},y_{i} \right)$ and ($x_{j},y_{j})$
with $i < j$ are “concordant” if we either observe $x_{i} < x_{j}$ and
$y_{i} < y_{j}$, or $x_{i} > x_{j}$ and $y_{i} > y_{j}$. Otherwise, they
are “discordant”. (Please note that, as $x_{i}$ and $y_{i}$ are unique,
there are no *ties*: it’s impossible that $x_{i} = x_{j}$ or
$y_{i} = y_{j}$). All possible pairwise evaluations are taken into
account: $np$ is the total number of pairs, given by
$\left( \frac{L}{2} \right)$, $nc$ is the number of total concordant
pairs, and $nd$ is the number of total discordant pairs.

The result of the Kendall’s $\tau$ kernel can range between -1 and 1. 1
means that the two rankings are identical: the two ordinal
samples/variables are equal, -1 that one ranking is the reverse of the
other, and 0 that the two rankings are completely independent.

#### Feature space

If we call $r_{i}$ and $s_{i}$ the ranks of the i-th member in the two
ordinal samples, we can then define two matrices **A** and **B** that:

$$\mathbf{A}_{ij} = sign\left( r_{j} - r_{i} \right)$$$$\mathbf{B}_{ij} = sign\left( s_{j} - s_{i} \right)$$
These matrices fulfill that $\mathbf{A}^{T} = - \mathbf{A}$ and that
$\mathbf{B}^{T} = - \mathbf{B}$. Now, we can observe that:

$$Kendall_{\tau}(x,y) = Frobenius_{cosine}(\mathbf{A},\mathbf{B})$$ From
there, the data in feature space can be found the same way than in the
Frobenius kernel. It can be deduced that, for the Kendall’s $\tau$
kernel, feature space has dimension $\left( \frac{L}{2} \right)$, and
each feature corresponds to a “relative order” ($+ = x_{i} < x_{j}$;
$- = x_{i} > x_{j}$) between all pair of levels within a ranking.

#### Usage

A *NxL* or *LxN* dataset is required. *N* is the sample size and *L* the
number of values (levels) of the ordinal variable. This is an example
that such a dataset:

**Three people are given a list of 10 colors. They rank them from most
(1) to least (10) favorite.**

``` r
color_list <-  c("black","blue","green","grey","lightblue","orange","purple",
"red","white","yellow")
survey1 <- 1:10
survey2 <- 10:1
survey3 <- sample(10)
color <- cbind(survey1,survey2,survey3) # Samples in columns
rownames(color) <- color_list
color
#>           survey1 survey2 survey3
#> black           1      10       5
#> blue            2       9      10
#> green           3       8       7
#> grey            4       7       8
#> lightblue       5       6       9
#> orange          6       5       2
#> purple          7       4       6
#> red             8       3       3
#> white           9       2       1
#> yellow         10       1       4
```

In this case, the ranking is already done. The Kendall’s $\tau$ kernel
is computed like:

``` r
Kendall(color)
```

(If we prefer samples in rows, we can do:)

``` r
color <- t(color)
Kendall(color,samples.in.rows=TRUE)
```

Alternatively, the ranking may be not given explicitly by the user.
Consider this other dataset:

**The same three people are asked the number of times they ate 5
different kinds of food during the last month**

``` r
food <- matrix(c(10, 1,18, 25,30, 7, 5,20, 5, 12, 7,20, 20, 3,22),ncol=3,nrow=5,byrow = TRUE)
colnames(food) <- colnames(color)
rownames(food) <- c("spinach", "chicken", "beef" , "salad","lentils")
food
#>         survey1 survey2 survey3
#> spinach      10       1      18
#> chicken      25      30       7
#> beef          5      20       5
#> salad        12       7      20
#> lentils      20       3      22
```

In that case, the ranking is done internally before computing the
Kendall’s $\tau$.

``` r
Kendall(food)
```

Furthermore, we can combine these two datasets, because the participants
are the same:

``` r
ordinal_data <- list(color=color,food=food) #All samples in columns
Kendall(ordinal_data)
```

By default, the Kendall’s $\tau$ kernel values of the two datasets are
averaged. The alternative is to sum them, doing:
`Kendall(...,comp="sum")`.

Right now,
[`Kendall()`](https://elies-ramon.github.io/kerntools/reference/Kendall.md)
does not return the data mapped onto the feature space.

## Kernels for strings, sequences, or (short) texts

Strings or sequences are vectors of objects (usually represented as
characters) in which repetitions are allowed and order matters. These
characters are drawn from a given alphabet $A$. An example of an
alphabet is:

``` r
LETTERS
#>  [1] "A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M" "N" "O" "P" "Q" "R" "S"
#> [20] "T" "U" "V" "W" "X" "Y" "Z"
```

We denote the size of an alphabet as $|A|$. In the example above, this
is:

``` r
length(LETTERS)
#> [1] 26
```

### Spectrum kernel

#### Feature space

Let $x$ be a string or sequence, and $s$ a substring of length *L* that
consist of characters drawn from an alphabet $A$, so we can count the
times ($t_{s}$) that $s$ occurs in $x$. Then, we do the same for all
possible substrings of length *L* generated from $A$. Each $t_{s}$ is a
feature in feature space, which has dimension $|A|^{L}$.

#### Definition

The Spectrum kernel is the dot product of two strings mapped onto the
feature space we have just defined. Alternatively, it can also be
written as:

$$Spectrum\left( x_{i},x_{j} \right) = \sum\limits_{s \in A^{L}}\left| s \in x_{i} \right| \cdot \left| s \in x_{j} \right|$$
This kernel ranges from 0 (the strings are completely different: they
don’t share any substring) to any positive number. In general, higher
values indicate higher similarity; however, it’s important to remember
that some strings may be longer than others. This is fine if the
strings’ length is relevant to the problem you are studying; but if it
is not, you can cosine-normalize the Spectrum kernel so the result is
bounded between 0 and 1:

$$Spectrum_{cosine}\left( x_{i},x_{j} \right) = \frac{Spectrum\left( x_{i},x_{j} \right)}{\sqrt{Spectrum\left( x_{i},x_{i} \right)Spectrum\left( x_{j},x_{j} \right)}}$$
A result of 0 means complete dissimilarity, 1 means complete similarity,
and the strings’ length is not taken into account.

#### Usage

You need a dataset containing *N* strings or text data. For instance:

``` r
alphabet <- c(letters,"_")
strings <- c("hello_world","hello_word","hola_mon","kaixo_mundua",
"saluton_mondo","ola_mundo", "bonjour_le_monde")
names(strings) <- c("english1","english_typo","catalan","basque",
"esperanto","galician","french")
strings
#>           english1       english_typo            catalan             basque 
#>      "hello_world"       "hello_word"         "hola_mon"     "kaixo_mundua" 
#>          esperanto           galician             french 
#>    "saluton_mondo"        "ola_mundo" "bonjour_le_monde"
alphabet
#>  [1] "a" "b" "c" "d" "e" "f" "g" "h" "i" "j" "k" "l" "m" "n" "o" "p" "q" "r" "s"
#> [20] "t" "u" "v" "w" "x" "y" "z" "_"
```

Then you can compute the Spectrum kernel. For instance, we could
consider substrings of length L = 2:

``` r
Spectrum(strings,l=2,alphabet=alphabet)
```

The feature space can be recovered doing:

``` r
Spectrum(strings,l=2,alphabet=alphabet, feat_space = TRUE)
```

The cosine-normalized version can be computed typing:

``` r
Spectrum(strings,l=2,alphabet=alphabet,cos.norm = TRUE)
```

(if you set `feat_space=TRUE`, the feature space returned by
[`Spectrum()`](https://elies-ramon.github.io/kerntools/reference/Spectrum.md)
will be scaled accordingly).

Right now, the Spectrum kernel implemented by `kerntools` is very basic
and, in large datasets and/or when *L* is high, the computation may be
excessively slow. In that case, there are other R packages that you may
use. For instance, you may refer to
[`kernlab::stringdot()`](https://rdrr.io/pkg/kernlab/man/stringdot.html),
or check the `kebabs` package, which covers complex kernels for strings
like the motif kernel, the mismatch kernel, the gappy pair kernel, etc.

## Kernels for bag-of-words (BoW) or bags-of-visual-words data

Given a document, we define a “bag” as the unordered collection of the
words it contains. A bag-of-words model counts the occurrence of each
word, disregarding word order. Thus, these datasets typically consist of
nonnegative *NxD* tables, where *D* is the total number of different
words present in all documents, and *N* is the number of documents.

### $\chi^{2}$ kernel

#### Definition

The $\chi^{2}$ kernel is defined as:

$$\chi^{2}\left( x_{i},x_{j} \right) = e^{- \gamma{\sum\limits_{d}\frac{{(x_{i}d - x_{j}d)}^{2}}{x_{i}d + x_{j}d}}}$$

This definition is akin to the RBF kernel, as the $\chi^{2}$ kernel is
also the negative exponential of a distance. As in the former, the
output ranges between 0 and 1: the higher the kernel value, the higher
the similarity between
$\mathbf{x}_{\mathbf{i}},\mathbf{x}_{\mathbf{j}}$. $\gamma > 0$ is a
hyperparameter that governs the “smoothness” of this kernel when
adapting to data.

#### Usage

You need a BoW dataset. For instance, consider the following two
documents:

- “John likes to watch movies. Mary likes movies too.”

- “Mary also likes to watch football games.”

If we count the words present in these two documents, we obtain the
following dataset:

``` r
words <-  c("John","likes","to","watch","movies","Mary","too","also","football","games")
documents <- matrix(c(1,2,1,1,2,1,1,0,0,0,0,1,1,1,0,1,0,1,1,1),nrow=2,ncol=length(words),byrow=TRUE)
colnames(documents) <- words
rownames(documents) <- 1:2
documents
#>   John likes to watch movies Mary too also football games
#> 1    1     2  1     1      2    1   1    0        0     0
#> 2    0     1  1     1      0    1   0    1        1     1
```

Now, to compute the $\chi^{2}$ kernel, we can type:

``` r
Chi2(documents,g=0.1)
#>           1         2
#> 1 1.0000000 0.4803053
#> 2 0.4803053 1.0000000
```

`Chi2(...,g)` should be chosen by the user. With $\gamma = 0.1$, the
similarity between documents 1 and 2 is almost 0.50. If
`Chi2(...,g=NULL)`, we obtain instead the LeCam distance (see [On
Measures of Entropy and
Information](https://threeplusone.com/pubs/on_information.pdf).

Finally, `kerntools`’s does not compute the feature space of the
$\chi^{2}$ kernel.
