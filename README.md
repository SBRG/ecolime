# ECOLIme

ECOLIme contains all the information and building scripts required to construct
a ME-model for E. coli K-12 MG1655, *i*JL1678b-ME using
[COBRAme](https://github.com/sbrg/cobrame). See the COBRAme 
[github](https://github.com/sbrg/cobrame) and 
[documentation](https://cobrame.readthedocs.io) for further instructions on 
how to build and solve ME-models.

**Note**: The current version of *i*JL1678b-ME is provided prior to peer-reviewed
publication. The model may be edited and updated until it is finalized for 
publication.


## Installation

1. clone the repository
2. run ```python setup.py develop --user```
3. Run build_me_model.py or build_me_model.ipynb to construct and save *i*JL1678b-ME
4. The saved model can be solved using solve_demo.ipynb


codon table source:
http://openwetware.org/wiki/Escherichia_coli/Codon_usage

