CGHNormaliter
-------------

In the analysis of aCGH data, normalization is deemed a critical pre-processing step. In general, aCGH normalization approaches are similar to those used for gene expression data, albeit both data-types differ inherently. A particular problem with aCGH data is that imbalanced copy numbers lead to improper normalization using conventional methods.

This issue is addressed by means of an iterative normalization procedure. First, provisory balanced copy numbers are identified and subsequently used for normalization. These two steps are then iterated to refine the normalization.

This component is actually just a python wrapper for the **R** implementation of the algorithm, which can be found at http://www.bioconductor.org/packages/release/bioc/html/CGHnormaliter.html. The package **rpy2** is required in order to include the original library, as well as the **R** code for `CGHNormaliter`.

See [CGHNormaliter]_ for a more detailed description of the algorithm.