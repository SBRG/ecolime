from setuptools import setup, find_packages


setup(name='ECOLIme',
      version='0.0.9',
      description='ME-model data files specific to E. coli K-12 MG1655',
      author='Colton Lloyd and Ali Ebrahim',
      url='https://github.com/SBRG/ecolime',
      install_requires=['Biopython', 'setuptools', 'cobra<=0.5.11', 'xlrd',
                        'pandas', 'six', 'scipy', 'numpy'],
      packages=find_packages(),
      package_data={'': ['building_data/*', 'me_models/*']},
      license='MIT')
