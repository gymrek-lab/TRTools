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
VERSION = "2.0.13" #version_file.read().strip()

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
                       'Programming Language :: Python :: 3.5',\
                       'License :: OSI Approved :: MIT License',\
                       'Operating System :: OS Independent',\
                       'Intended Audience :: Science/Research',\
                       'Topic :: Scientific/Engineering :: Bio-Informatics'],
	  data_files=[
		('', ['LICENSE.txt'])
	  ]
     )
