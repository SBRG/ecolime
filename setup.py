from setuptools import setup, find_packages


setup(name='EcoliME',
      version='0.0.5',
      description='ME model data files specific to E. coli K-12 MG1655',
      author='Ali Ebrahim and Colton Lloyd',
      author_email='minime_dev@googlegroups.com',
      url='https://github.com/SBRG/ecolime',
      install_requires=['Biopython', 'setuptools', 'cobra<=5.11'],
      packages=find_packages(),
      package_data={'': ['building_data/*', 'me_models/*']},
      license='MIT')
