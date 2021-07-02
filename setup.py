from setuptools import setup, find_packages
from who2tbp import __version__ as VERSION

setup(
    name='who2tbp',
    version=VERSION,
    packages=find_packages(),
    url='github.com/MDU-PHL/who2tbp',
    license='GPLv3',
    author='Anders Goncalves da Silva',
    author_email='andersgs@gmail.com',
    description='A script to convert WHO Excel mutations to TBProfiler database format',
    install_requires=[
        'biopython>=1.79',
        'openpyxl>=3.0.7',
        'setuptools>=49.6.0',
        'tqdm>=4.61.1',
        'gffutils>=0.10.1'
    ],
    entry_points={
        "console_scripts":
            [
                'who2tbp=who2tbp.bin.main:main'
            ]
    },
    package_data={'who2tbp': ['data/genome.gff']},
    include_package_data=True,
)
