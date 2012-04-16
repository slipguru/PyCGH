from setuptools import setup, find_packages

setup(
    name = 'PyCGH',
    use_hg_version = True,
    packages = find_packages(),

    # Dependencies
    setup_requires=['hgtools'],
    #install_requires = [],
    requires = ['numpy (>=1.4.0)'], # More??

    # Data Files
    #package_data = { }

    # Metadata
    #url='http://slipguru.disi.unige.it/Software/L1L2Py',
    #description=description,
    #long_description=long_description,
    #keywords='feature selection, regularization, regression, classification,'
    #         'l1-norm, l2-norm',
    #author='L1L2Py developers - SlipGURU',
    #author_email='slipguru@disi.unige.it',
    #maintainer='Salvatore Masecchia',
    #maintainer_email='salvatore.masecchia@disi.unige.it',
    #license='GNU GPL version 3',
    #download_url = 'http://slipguru.disi.unige.it/Software/L1L2Py/L1L2Py-%s.tar.gz' % version,
    #classifiers=[
    #    'Development Status :: 5 - Production/Stable',
    #    'Intended Audience :: Science/Research',
    #    'License :: OSI Approved :: GNU General Public License (GPL)',
    #    'Natural Language :: English',
    #    'Operating System :: POSIX :: Linux',
    #    'Operating System :: MacOS',
    #    'Operating System :: Microsoft :: Windows',
    #    'Programming Language :: Python',
    #    'Topic :: Scientific/Engineering :: Artificial Intelligence',
    #],
)
