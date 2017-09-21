# ECOLIme

ECOLIme contains all the information and building scripts required to construct
a ME-model for E. coli K-12 MG1655, *i*LE1678-ME using
[COBRAme](https://github.com/sbrg/cobrame). See the COBRAme 
[github](https://github.com/sbrg/cobrame) for further instructions on how to 
build and solve ME-models.

**Note**: The current version of *i*LE1678-ME is provided prior to peer-reviewed
publication. The model will be edited and updated until it is finalized for 
publication.


## Installation

1. clone the repository
2. run ```python setup.py develop --user```
3. Run build_ME_model.py or build_ME_model.ipynb to construct and save *i*LE1678-ME
4. The saved model can be solved using solve_demo.ipynb


codon table source:
http://openwetware.org/wiki/Escherichia_coli/Codon_usage

