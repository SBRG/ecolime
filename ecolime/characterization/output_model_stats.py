import pandas as pd
from cobrame import MetabolicReaction


def output_model_reactions_stats(model, output_file):
    """Return general statistics about the reactions in the ME-model"""
    df = pd.DataFrame(columns=['genes', 'products_stoichiometry', 'products',
                               'reactants_stoichiometry', 'reactants', 'rxn_id',
                               'reaction_type', 'reaction'],
                      index=["index"])
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
            except:
                print r.id

        df.loc[r_index, 'products_stoichiometry'] = product_stoich
        df.loc[r_index, 'products'] = products
        df.loc[r_index, 'reactants_stoichiometry'] = reactant_stoich
        df.loc[r_index, 'reactants'] = reactants
        df.loc[r_index, 'rxn_id'] = r.id
        df.loc[r_index, 'reaction'] = r.reaction

        gene_list = []
        if isinstance(r, MetabolicReaction) and r.complex_data:
            for met in r.complex_data.stoichiometry:
                gene_list.append(met.replace('protein_', ''))

        df.loc[r_index, 'genes'] = gene_list
        df.loc[r_index, 'reaction_type'] = r.__class__.__name__

    df.to_csv(output_file)


def output_model_component_stats(model, output_file):
    """Return general statistics about the components in the ME-model"""
    df = pd.DataFrame(columns=['met_id', 'met_name', 'formula', 'charge',
                               'met_type'],
                      index=["index"])
    for m in model.metabolites:
        m_index = model.metabolites.index(m.id)

        df.loc[m_index, 'met_id'] = m.id
        df.loc[m_index, 'met_name'] = m.name
        df.loc[m_index, 'formula'] = m.formula
        df.loc[m_index, 'charge'] = m.charge
        df.loc[m_index, 'met_type'] = m.__class__.__name__

    df.to_csv(output_file)
