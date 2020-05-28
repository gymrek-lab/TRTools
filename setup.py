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

#version_file = open('_version')
VERSION = "2.0.16test" #version_file.read().strip()

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
      scripts=["trtools/testsupport/test_trtools.sh"],
      entry_points={
          'console_scripts': [
              'dumpSTR=trtools.dumpSTR:run',
              'mergeSTR=trtools.mergeSTR:run',
              'statSTR=trtools.statSTR:run',
              'compareSTR=trtools.compareSTR:run',
              'qcSTR=trtools.qcSTR:run'
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
                       'Programming Language :: Python :: 3.5',\
                       'License :: OSI Approved :: MIT License',\
                       'Operating System :: OS Independent',\
                       'Intended Audience :: Science/Research',\
                       'Topic :: Scientific/Engineering :: Bio-Informatics']
     )
