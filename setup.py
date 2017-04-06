from setuptools import setup, find_packages


setup(name='EcoliME',
      version='0.0.5',
      description='ME model data files specific to E. Coli K12-MG1655',
      author='Ali Ebrahmim and Colton Lloyd',
      author_email='minime_dev@googlegroups.com',
      url='https://github.com/SBRG/ecolime',
      install_requires=['Biopython', 'setuptools'],
      packages=find_packages(),
      package_data={'': ['building_data/*']},
      license='MIT')
