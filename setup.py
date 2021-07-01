from setuptools import setup

setup(
    name='who2tbp',
    version='0.1.0',
    packages=['who2tbp'],
    url='github.com/MDU-PHL/who2tbp',
    license='GPLv3',
    author='Anders Goncalves da Silva',
    author_email='andersgs@gmail.com',
    description='A script to convert WHO Excel mutations to TBProfiler database format',
    entry_points={
        "console_scripts":
            [
                'who2tbp=who2tbp.bin.main:main'
            ]
    }
)
