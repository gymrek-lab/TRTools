import os
from setuptools import setup, find_packages

DESCRIPTION = "Toolkit for genome-wide analysis of STRs"
LONG_DESCRIPTION = DESCRIPTION
NAME = "trtools"
AUTHOR = "Melissa Gymrek"
AUTHOR_EMAIL = "mgymrek@ucsd.edu"
MAINTAINER = "Melissa Gymrek"
MAINTAINER_EMAIL = "mgymrek@ucsd.edu"
DOWNLOAD_URL = 'http://github.com/gymreklab/TRTools'
LICENSE = 'MIT'

# version-keeping code based on pybedtools
curdir = os.path.abspath(os.path.dirname(__file__))
MAJ = 5
MIN = 0
REV = 1
VERSION = '%d.%d.%d' % (MAJ, MIN, REV)
with open(os.path.join(curdir, 'trtools/version.py'), 'w') as fout:
        fout.write(
            "\n".join(["",
                       "# THIS FILE IS GENERATED FROM SETUP.PY",
                       "version = '{version}'",
                       "__version__ = version"]).format(version=VERSION)
        )

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=DOWNLOAD_URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      python_requires='>=3.5',
      packages=find_packages(),
      include_package_data=True,
      license_file="LICENSE.txt",
      scripts=[
          "trtools/testsupport/test_trtools.sh",
          'scripts/trtools_prep_beagle_vcf.sh'
      ],
      entry_points={
          'console_scripts': [
              'dumpSTR=trtools.dumpSTR:run',
              'mergeSTR=trtools.mergeSTR:run',
              'statSTR=trtools.statSTR:run',
              'compareSTR=trtools.compareSTR:run',
              'qcSTR=trtools.qcSTR:run',
              'associaTR=trtools.associaTR:run'
          ],
      },
      install_requires=['cyvcf2',
                        'matplotlib',
                        'numpy',
                        'pandas',
                        'pybedtools',
                        'pysam',
                        'scikit-learn',
                        'scipy',
                        'statsmodels'],
      classifiers=['Development Status :: 4 - Beta',\
                       'Programming Language :: Python :: 3.5',\
                       'License :: OSI Approved :: MIT License',\
                       'Operating System :: OS Independent',\
                       'Intended Audience :: Science/Research',\
                       'Topic :: Scientific/Engineering :: Bio-Informatics']
     )
