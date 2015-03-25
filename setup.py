try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name='EcoliME',
      version='0.01',
      description='ME model data files specific to EColi K12-MG1655',
      author='Ali Ebrahmim and Colton Lloyd',
      author_email='minime_dev@googlegroups.com',
      url='https://github.com/SBRG/ecolime',
      install_requires=['Biopython'],
      packages=['ecolime'],
      )
