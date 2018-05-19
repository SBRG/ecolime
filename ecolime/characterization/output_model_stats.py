from collections import defaultdict

import pandas as pd

from cobrame import MetabolicReaction


def output_model_reactions_stats(model, output_file, solution=None):
    if solution:
        model.solution = solution

    """Return general statistics about the reactions in the ME-model"""
    df = pd.DataFrame(columns=['id', 'genes', 'products_stoichiometry',
                               'products', 'reactants_stoichiometry',
                               'reactants', 'reaction_type', 'reaction',
                               'upper_bound', 'lower_bound'],
                      index=["index"])
    writer = pd.ExcelWriter(output_file)
    rxn_types = set()
    for r in model.reactions:
        r_index = model.reactions.index(r.id)
        products = []
        reactants = []
        product_stoich = []
        reactant_stoich = []
        for met, stoich in r.metabolites.items():
            try:
                if stoich < 0:
                    reactants.append(met.id)
                    reactant_stoich.append(stoich)
                else:
                    products.append(met.id)
                    product_stoich.append(stoich)
            except TypeError:
                reactants.append(met.id)
                reactant_stoich.append(stoich)

        df.loc[r_index, 'products_stoichiometry'] = product_stoich
        df.loc[r_index, 'products'] = products
        df.loc[r_index, 'reactants_stoichiometry'] = reactant_stoich
        df.loc[r_index, 'reactants'] = reactants
        df.loc[r_index, 'id'] = r.id
        df.loc[r_index, 'lower_bound'] = str(r.lower_bound)
        df.loc[r_index, 'upper_bound'] = str(r.upper_bound)

        # Add flux if model has been solved
        if model.solution.x_dict:
            df.loc[r_index, 'flux'] = model.solution.x_dict[r.id]
        try:
            df.loc[r_index, 'reaction'] = r.reaction
        except TypeError:
            ub, lb = r.upper_bound, r.lower_bound
            r.upper_bound = 1000.
            r.lower_bound = 0.
            df.loc[r_index, 'reaction'] = r.reaction
            print(r, r.reaction)
            r.upper_bound = ub
            r.lower_bound = lb
            print(r.upper_bound)

        gene_list = []
        if isinstance(r, MetabolicReaction) and r.complex_data:
            for met in r.complex_data.stoichiometry:
                gene_list.append(met.replace('protein_', ''))

        df.loc[r_index, 'genes'] = gene_list
        rxn_type = r.__class__.__name__
        df.loc[r_index, 'reaction_type'] = rxn_type
        rxn_types.add(rxn_type)

    for rxn_type in rxn_types:
        df_filtered = df[df.rxn_type == rxn_type]
        df_filtered.to_excel(writer, sheet_name=rxn_type)

    writer.save()


def output_model_component_stats(model, output_file):
    """Return general statistics about the components in the ME-model"""
    df = pd.DataFrame(columns=['met_id', 'met_name', 'formula', 'charge',
                               'met_type'],
                      index=["index"])
    writer = pd.ExcelWriter(output_file)
    met_types = set()
    for m in model.metabolites:
        m_index = model.metabolites.index(m.id)
        for key, value in m.__dict__.items():
            if not key.startswith('_'):
                df.loc[m_index, key] = str(value)

        met_type = m.__class__.__name__
        df.loc[m_index, 'met_type'] = met_type
        met_types.add(met_type)

        rxn_involvement_set = set()
        for r in m.reactions:
            rxn_involvement_set.add(r.__class__.__name__)
        df.loc[m_index, 'reaction_involvement'] = str(rxn_involvement_set)

    for met_type in met_types:
        df_filtered = df[df.met_type == met_type]
        df_filtered = df_filtered.T.dropna(how='all').T
        df_filtered.set_index('id')
        df_filtered.to_excel(writer, sheet_name=met_type)

    writer.save()


def output_model_process_data_stats(model, output_file):
    """Return general statistics about the components in the ME-model"""
    df = pd.DataFrame()
    writer = pd.ExcelWriter(output_file)
    data_types = set()
    for d in model.process_data:
        d_index = model.process_data.index(d.id)
        for key, value in d.__dict__.items():
            if isinstance(value, defaultdict):
                df.loc[d_index, key] = str(dict(value))
            else:
                df.loc[d_index, key] = str(value)
        data_type = d.__class__.__name__
        df.loc[d_index, 'data_type'] = data_type
        data_types.add(data_type)

    for data_type in data_types:
        df_filtered = df[df.data_type == data_type]
        df_filtered = df_filtered.T.dropna(how='all').T
        df_filtered.set_index('id')
        df_filtered.to_excel(writer, sheet_name=data_type)

    writer.save()
