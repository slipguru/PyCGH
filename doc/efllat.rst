eFLLAT
------

Copy number variations (CNVs) are alterations of the DNA that result in the cell having an abnormal number of copies of one or more sections of the DNA. Recurrent aberrations across samples may indicate an oncogene or a tumor suppressor gene, but the functional mechanisms that link altered copy numbers to pathogenesis are still to be explained. Array-based Comparative Genomic Hybridization (aCGH) is a modern whole-genome measuring technique that evaluates the occurrence of copy variants across the genome of samples (patients) versus references (controls) on the entire genome, extending the original CGH technology.

A signal measured with an aCGH is made of a piecewise linear (and constant) component plus some noise. The typical analysis on such data is segmentation, that is the automatic detection of altered copy numbers (amplifications or deletions). Differently from other molecular data, as gene expression, with aCGH it is possible to exploit the intrinsic data structure to improve the downstream analysis.

The algorithm uses a Dictionary Learning (DL) approach to identify common structures in the aCGH samples, i.e. locations where alterations (deletions or duplications) are significantly more frequent than in the rest of the genome by minimizing functional

.. math::

    {\bm B}, {\bm \Theta} = \argmin \| {\bm Y} - {\bm B} {\bm \Theta} \|_F^2 + \tau \sum_{s=1}^S \| {\bm \Theta}(:,s) \|_1^2 + \lambda \sum_{j=1}^J \| {\bm B}(:, j) \|_1^2 + \mu \sum_{j=1}^J TV_w ({\bm B}(:, j))

For more details on the algorithm see [masecchia13]_.

.. [masecchia13] _ S. Masecchia, S. Salzo, A. Barla and A. Verri. A dictionary learning based method for aCGH segmentation. *Proceedings of the European Symposium on Artificial Neural Networks, Computational Intelligence and Machine Learning*, 2013.