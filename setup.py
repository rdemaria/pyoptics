from setuptools import setup, find_packages

setup(
        name='pyoptics',
        version='0.0.1',
        description='Particle Beam Optics Tools',
        author='Riccardo De Maria',
        author_email='riccardo.de.maria@cern.ch',
        url='https://github.com/rdemaria/pyoptics',
        packages=find_packages(),
        install_requires=['numpy','matplotlib','cpymad'],
)
