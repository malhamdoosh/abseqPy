from setuptools import setup, find_packages

import subprocess
import sys


def readme():
    with open("README.md") as f:
        return f.read()


# weird bug where setup_requires numpy is not obeyed
setup_requires = ['numpy>=1.11.3']
for pack in setup_requires:
    # pip.main(['install', pack]) no longer supported in pip >= 10
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', pack])

setup(name="abseqPy",
      version="1.1.15",
      description="Quality control pipeline for antibody libraries",
      license="placeholder",
      long_description=readme(),
      author="CSL",
      author_email="placeholder",
      maintainer="CSL",
      maintainer_email="placeholder",
      # pandas requires numpy installed, it's a known bug in setuptools - put in both setup and install requires
      # UPDATE Wed Feb 21 13:15:43 AEDT 2018 - moved into pre-installation stage
      setup_requires=['numpy>=1.11.3', 'pandas>=0.20.1', 'biopython>=1.66', 'weblogo>=3.4', 'matplotlib>=1.5.1',
                      'tables>=3.2.3.1', 'psutil', 'matplotlib-venn', 'pyyaml'],
      install_requires=['numpy>=1.11.3', 'pandas>=0.20.1', 'biopython>=1.66', 'weblogo>=3.4', 'matplotlib>=1.5.1',
                        'tables>=3.2.3.1', 'psutil', 'matplotlib-venn', 'pyyaml', 'scipy>=0.19.0'],
      packages=find_packages(),
      include_package_data=True,
      entry_points={
          'console_scripts': ['abseq=abseqPy.abseqQC:main'],
      })
