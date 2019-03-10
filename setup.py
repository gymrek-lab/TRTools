from setuptools import setup, find_packages

DESCRIPTION = "Toolkit for genome-wide analysis of STRs"
LONG_DESCRIPTION = DESCRIPTION
NAME = "strtools"
AUTHOR = "Melissa Gymrek"
AUTHOR_EMAIL = "mgymrek@ucsd.edu"
MAINTAINER = "Melissa Gymrek"
MAINTAINER_EMAIL = "mgymrek@ucsd.edu"
DOWNLOAD_URL = 'http://github.com/gymreklab/STRTools'
LICENSE = 'MIT'

VERSION = '1.0.0'

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
              'plinkSTR=plinkSTR.plinkSTR:main',
              'dumpSTR=dumpSTR.dumpSTR:main',
              'mergeSTR=mergeSTR.mergeSTR:main',
              'statSTR=statSTR.statSTR:main'
          ],
      },
      install_requires=['argparse',
                        'numpy',
                        'pybedtools',
                        'pyvcf'],
      classifiers=['Development Status :: 4 - Beta',\
                       'Programming Language :: Python :: 3.2',\
                       'License :: OSI Approved :: MIT License',\
                       'Operating System :: OS Independent',\
                       'Intended Audience :: Science/Research',\
                       'Topic :: Scientific/Engineering :: Bio-Informatics']
     )
