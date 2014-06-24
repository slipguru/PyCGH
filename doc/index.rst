.. |project_name| replace:: **PyCGH**

===================================================
PyCGH - a Comparative Genomic Hybridization toolkit
===================================================

:Release: |release|
:Homepage: http://slipguru.disi.unige.it/Software/PyCGH
:Repository: https://bitbucket.org/slipguru/pycgh

|project_name| is a Python library for the analysis of aCGH data.
It consists mainly of three components:

 * Synthetic data generation: a script which creates synthetic aCGH data, for testing purposes.
 * CGHNormaliter: A python wrapper for the *CGHNormaliter* algorithm [CGHNormaliter]_, a preprocessing step required to normalizeCGH signals.
 * E-FLLat: an algorithm which uses a dictionary learning approach to discover common patterns in aCGH data.

|project_name| requires the following python libraries:

  * `numpy <http://www.numpy.org/>`_: the fundamental package for scientific computing.
  * `matplotlib <http://matplotlib.org/>`_ : a plotting library.
  * `rpy2 <http://rpy.sourceforge.net/>`_ : a library to include and use **R** code in a python program.

.. note:: Since the function implementing the *CGHNormaliter* algorithm is actually a wrapper for the implementation in the **R** language, an appropriate interpreter for that language is required, as well as the library containing the implementation of the algorithm itself, which can be found at `this link <http://www.bioconductor.org/packages/release/bioc/html/CGHnormaliter.html>`_.

User Documentation
==================

In this Section an overview of each component will be given.


.. toctree::
   :maxdepth: 2

   synthesis
   cghnormaliter
   efllat

PyCGH API
==========

.. toctree::
   :maxdepth: 2

   api.rst

:ref:`genindex`

.. Quick Reference
.. ---------------
.. currentmodule:: pycgh


.. rubric:: References

.. [CGHNormaliter] B\. van Houte, T. Binsl, H. Hettling, W. Pirovano and J. Heringa. CGHnormaliter: an iterative strategy to enhance normalization of array CGH data with imbalanced aberrations. *BMC Genomics*, 2009.