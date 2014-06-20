CGHNormaliter
-------------

In the analysis of aCGH data, normalization is deemed a critical pre-processing step. In general, aCGH normalization approaches are similar to those used for gene expression data, albeit both data-types differ inherently. A particular problem with aCGH data is that imbalanced copy numbers lead to improper normalization using conventional methods.

This issue is addressed by means of an iterative normalization procedure. First, provisory balanced copy numbers are identified and subsequently used for normalization. These two steps are then iterated to refine the normalization.

This component is actually just a python wrapper for the **R** implementation of the algorithm, which can be found at http://www.bioconductor.org/packages/release/bioc/html/CGHnormaliter.html. The package **rpy2** is required in order to include the original library, as well as the **R** code for `CGHNormaliter`.

To use this component it is sufficient to import in the main script the module containing the CGHNormaliter wrapper, and then call function ``cghnormaliter`` passing a :class:`pycgh.datatypes.ArrayCGH` object as parameter: the following example shows a script which first creates a aCGH object using the module for the generation of synthetic data and then starts the normalization, returning the normalized aCGH as a result [#f1]_.

.. code-block:: python

    from pycgh.analysis.cghnormaliter import cghnormaliter
    from pycgh import synth, datatypes, readers
    
    def main():
        # Read chip definition (file downloaded from UCSC Genome Browser)
        chip_design = readers.ucsc_mapping('agilentCgh4x44k.txt.gz',filter_valid=True)
        
        # Read cytobands information (file downloaded from UCSC Genome Browser)
        cs = datatypes.CytoStructure('cytoBandIdeo.txt.gz')
        
        # Create aCGH synthetizer
        source = synth.ArrayCGHSynth(geometry=(430, 103),
            design=chip_design,
            cytostructure=cs,
            alterations = {'17q': [(3, .9)]}) # Gain
        
        # Generate a male sample
        acgh = source.draw('male')
        
        normalized_acgh = cghnormaliter(acgh)

See [CGHNormaliter]_ for a more detailed description of the algorithm.

.. rubric:: Footnotes

.. [#f1] The two files ``agilentCgh4x44k.txt.gz`` and ``cytoBandIdeo.txt.gz`` must be in the same folder as the script.
