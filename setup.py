from setuptools import setup, find_packages
import os

DESCRIPTION = "Toolkit for genome-wide analysis of STRs"
LONG_DESCRIPTION = DESCRIPTION
NAME = "trtools"
AUTHOR = "Melissa Gymrek"
AUTHOR_EMAIL = "mgymrek@ucsd.edu"
MAINTAINER = "Melissa Gymrek"
MAINTAINER_EMAIL = "mgymrek@ucsd.edu"
DOWNLOAD_URL = 'http://github.com/gymreklab/TRTools'
LICENSE = 'MIT'

##
# version-keeping code based on pybedtools
curdir = os.path.abspath(os.path.dirname(__file__))
MAJ = 2
MIN = 0
REV = 5
VERSION = '%d.%d.%d' % (MAJ, MIN, REV)
with open(os.path.join(curdir, 'trtools/version.py'), 'w') as fout:
        fout.write(
            "\n".join(["",
                       "# THIS FILE IS GENERATED FROM SETUP.PY",
                       "version = '{version}'",
                       "__version__ = version"]).format(version=VERSION)
        )
##

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
      packages=find_packages(),
      entry_points={
          'console_scripts': [
              'dumpSTR=dumpSTR.dumpSTR:run',
              'mergeSTR=mergeSTR.mergeSTR:run',
              'statSTR=statSTR.statSTR:run',
              'compareSTR=compareSTR.compareSTR:run',
              'qcSTR=qcSTR.qcSTR:run'
          ],
      },
      install_requires=['argparse',
                        'matplotlib',
                        'numpy',
                        'pandas',
                        'pybedtools',
                        'pyvcf',
                        'scipy',
                        'pysam'],
      classifiers=['Development Status :: 4 - Beta',\
                       'Programming Language :: Python :: 3.2',\
                       'License :: OSI Approved :: MIT License',\
                       'Operating System :: OS Independent',\
                       'Intended Audience :: Science/Research',\
                       'Topic :: Scientific/Engineering :: Bio-Informatics']
     )
