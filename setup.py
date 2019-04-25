from setuptools import setup, find_packages, Extension

import subprocess
import sys
import platform


def readme():
    with open("README.md") as f:
        return f.read()


def windows_filter(package):
    """
    TAMO not supported in windows:

    MDsupport.cxx
    c:\users\hello\appdata\local\temp\pip-dgi1lr-build\tamo\md\mdsupport_source\MDsupport.h(15) :
    fatal error C1083: Cannot open include file: 'regex.h': No such file or directory

    :param package: package name
    :return: True if this package should be installed in the host OS (True for all packages in non-windows machines)
    Windows machines will not install TAMO
    """
    if platform.system() == "Windows":
        return 'TAMO' not in package.upper()
    return True


# weird bug where setup_requires numpy is not obeyed
setup_requires = ['numpy>=1.11.3']
for pack in setup_requires:
    # pip.main(['install', pack]) no longer supported in pip >= 10
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', pack])

setup(name="abseqPy",
      version="0.99.2",
      description="Quality control pipeline for antibody libraries",
      long_description=readme(),
      long_description_content_type="text/markdown",
      author="Monther Alhamdoosh",
      author_email="m.hamdoosh@gmail.com",
      url="https://github.com/malhamdoosh/abseqPy",
      maintainer="Jia Hong Fong",
      maintainer_email="jiahfong@gmail.com",
      # pandas requires numpy installed, it's a known bug in setuptools - put in both setup and install requires
      # UPDATE Wed Feb 21 13:15:43 AEDT 2018 - moved into pre-installation stage
      setup_requires=['numpy>=1.11.3', 'pandas>=0.20.1', 'biopython>=1.66', 'matplotlib>=1.5.1, <3.0',
                      'tables>=3.2.3.1', 'psutil', 'matplotlib-venn', 'pyyaml'] +
                     ['weblogo>=3.4'] if platform.system() != 'Windows' else [],
      install_requires=['numpy>=1.11.3', 'pandas>=0.20.1', 'biopython>=1.66', 'matplotlib>=1.5.1, <3.0',
                        'tables>=3.2.3.1', 'psutil', 'matplotlib-venn', 'pyyaml', 'scipy>=0.19.0'] +
                       ['weblogo>=3.4'] if platform.system() != 'Windows' else [],
      packages=filter(windows_filter, find_packages()),
      ext_modules=([
                       Extension('TAMO.MD._MDsupport',
                                 ['TAMO/MD/MDsupport_source/MDsupport.cxx',
                                  'TAMO/MD/MDsupport_source/MDsupport_wrap.cxx']),
                       Extension('TAMO.util._swilk',
                                 ['TAMO/util/swilk_source/swilk.cxx',
                                  'TAMO/util/swilk_source/swilk_wrap.cxx'])
                   ] if platform.system() != "Windows" else []),
      include_package_data=True,
      entry_points={
          'console_scripts': ['abseq=abseqPy.abseqQC:main'],
      },
      classifiers=[
            "Programming Language :: Python :: 2.7",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Operating System :: MacOS",
            "Operating System :: Microsoft",
            "Operating System :: POSIX",
            "Operating System :: Unix",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Intended Audience :: Science/Research",
            "Environment :: Console"
      ]
)
