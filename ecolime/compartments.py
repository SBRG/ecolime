import cobrame
from ecolime.util.helper_functions import get_base_complex_data
from collections import defaultdict


def _return_compartments_of_complexes(model, cplx):
    try:
        data = model.process_data.get_by_id(cplx.id)
    except KeyError:
        data = get_base_complex_data(model, cplx.id)

    mem_dict = defaultdict(int)
    for s in data.stoichiometry:
        if '_Inner_Membrane' in s:
            mem_dict['im'] += 1
        elif '_Outer_Membrane' in s:
            mem_dict['om'] += 1
        elif '_Periplasm' in s:
            mem_dict['p'] += 1

    # if no membrane associated with membrane subunit, assume complex is
    # cytosolic
    if len(mem_dict) == 0:
        return 'c'
    # if only one membrane is represented in protein subunits, use this
    # membrane for the compartment
    elif len(mem_dict) == 1:
        return mem_dict.popitem()[0]
    # if multiple membrane compartements are represented, use generic "m" for
    # "membrane" for now
    else:
        return 'm'


def add_compartments_to_model(model):
    """Firsts adds compartments based on suffix of metabolite ID. If metabolite
    is a complex, the protein subunit stoichiometry is used to infer
    compartment. All remaining metabolites without a compartment suffix (RNAs,
    generic metabolites, nonmembrane proteins, etc.) are assumed to be
    cytosolic"""

    for m in model.metabolites:
        if m.compartment:
            continue
        if isinstance(m, cobrame.Constraint):
            m.compartment = 'mc'
        elif '_Inner_Membrane' in m.id:
            m.compartment = 'im'
        elif '_Outer_Membrane' in m.id:
            m.compartment = 'om'
        elif '_Periplasm' in m.id or m.id.endswith('_p'):
            m.compartment = 'p'
        elif m.id.endswith('_e'):
            m.compartment = 'e'
        elif m.id.endswith('_c'):
            m.compartment = 'c'
        elif isinstance(m, cobrame.Complex):
            m.compartment = _return_compartments_of_complexes(model, m)
        else:
            m.compartment = 'c'
