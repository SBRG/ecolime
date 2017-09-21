from collections import Counter

import cobrame
from cobrame.util import massbalance


def get_remaining_complex_elements(model, complex):

    tmp_met = cobrame.Metabolite('tmp_met')
    mets = model.metabolites
    components = complex.id.split('_mod_')
    base_complex = components[0]
    elements = Counter()

    # If a the completely unmodified complex is present in the model
    # has a formula, initialize the elements dictionary with that
    if base_complex in mets and mets.get_by_id(base_complex).formula:
        elements.update(mets.get_by_id(base_complex).elements)
    for component in components[1:]:
        new_elements = elements.copy()
        new_complex = '_mod_'.join([base_complex, component])
        if new_complex in mets and mets.get_by_id(new_complex).formula:

            # default to new_complex elements if both new and old exist
            if base_complex in mets and mets.get_by_id(base_complex).formula:
                new_elements = Counter()
            formula = mets.get_by_id(new_complex).formula
            tmp_met.formula = formula
            new_elements.update(tmp_met.elements)

        # Net effect of an SH modification is adding a Sulfur to elements
        elif ':SH' in component:
            new_elements['S'] += 1

        # modifies O- to SH
        elif component == 'cosh':
            new_elements['O'] -= 1
            new_elements['S'] += 1
            new_elements['H'] += 1

        elif component in model.global_info['modification_formulas']:
            formula = model.global_info['modification_formulas'][component]
            tmp_met.formula = formula
            new_elements.update(tmp_met.elements)

        elif ':' in component:
            value, component = component.split(':')
            if component in model.global_info['modification_formulas']:
                formula = \
                    model.global_info['modification_formulas'][component]['formula']
            elif component + '_c' in mets:
                formula = mets.get_by_id(component + '_c').formula
            else:
                raise UserWarning('No formula found for modification (%s)'
                                  % component)
            tmp_met.formula = formula

            for e, v in tmp_met.elements.items():
                new_elements[e] += v * float(value)

        elif 'Oxidized' in component and 'FLAVODOXIN' not in base_complex:
            new_elements.update({'H': -2})

        if elements == new_elements and 'FLAVODOXIN' not in base_complex:
            print(complex.id, base_complex, component)
        base_complex = '_mod_'.join([base_complex, component])
        elements = new_elements.copy()

    return elements


def add_remaining_complex_formulas(model):
    """
        Add formula to complexes that are not formed from a complex formation
        reaction (ie. complexes involved in metabolic reactions)
    """
    element_dict = {}

    # Reset all formulas first
    complex_list = []
    for c in model.metabolites:
        if not isinstance(c, cobrame.Complex) or c.id in model.complex_data:
            continue
        for r in c.reactions:
            if hasattr(r, 'update'):
                r.update()

        c.formula = ''
        c.elements = {}
        complex_list.append(c)

    # Get formulas only for complexes without complex formation reaction
    for c in complex_list:
        element_dict[c] = get_remaining_complex_elements(model, c)

    # Adding elements for complexes dynamically can change function output
    # Update all of them after
    for c, elements in element_dict.items():
        massbalance.elements_to_formula(c, elements)

