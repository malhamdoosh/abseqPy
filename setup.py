import os
import sys
from setuptools import setup, find_packages, Command



# https://stackoverflow.com/questions/3779915/why-does-python-setup-py-sdist-create-unwanted-project-egg-info-in-project-r
class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        #os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')
        os.system('rm -vrf ./build ./dist ./*.egg-info')

setup(name="AbSeq",
    version="1.1.1",
    description="Quality control pipeline for antibody libraries",
    author="CSL",
    setup_requires=['numpy>=1.11.3', 'pytz', 'python-dateutil'],
    install_requires=['pandas>=0.20.1', 'biopython==1.66', 'weblogo>=3.4'] + (['psutil'] if "darwin" in sys.platform else []),
    packages=find_packages(),
	cmdclass={
		'clean':CleanCommand,
	},
)
