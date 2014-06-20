Synthetic data generation
-------------------------

This component consists of a module, :mod:`pycgh.synth`, which allows to easily generate synthetic Array CGH data which resemble real data.
The core class is :class:`pycgh.synth.ArrayCGHSynth`, which must be instantiated passing all parameters required to generate the data; part of these arguments describe the structure of the data, for instance which probes (i.e. which DNA fragments) will be present on the simulated chip, as well as the structure of chromosomes, that is their division in regions (bands), other parameters control different characteristics such as the presence of alterations, noise or spatial bias.

The design of the chip is defined in an external file which must be provided to the script and is elaborated by auxiliary function :func:`pycgh.readers.ucsc_mapping`.

Once instantiated, the method :meth:`pycgh.synth.ArrayCGHSynth.draw` can be used to generate a :class:`pycgh.datatypes.ArrayCGH` object representing the Array CGH with the required characteristics.

The following is an example script which generates synthetic aCGH data (variable ``acgh``) and plots its profile, the relative MA plot and a representation of the virtual chip. The required files ``agilentCgh4x44k.txt.gz`` and ``cytoBandIdeo.txt.gz`` must be in the same folder as the script.

.. code-block:: python

    import pylab as pl # standard plotting library
    
    from pycgh import synth, datatypes, readers, plots

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
        
        # Profile plot (top)
        pl.subplot2grid((2,2),(0, 0), colspan=2)
        plots.profile(acgh)
        
        # Spatial plot (bottom-left)
        pl.subplot2grid((2,2), (1, 0))
        plots.spatial(acgh)
        
        # MA-plot (bottom-right)
        pl.subplot2grid((2,2), (1, 1))
        plots.MA(acgh)
        
        # Save to aCGH to file
        acgh.save('acgh_file.npz')
        pl.show()

    if __name__ == '__main__':
        main()