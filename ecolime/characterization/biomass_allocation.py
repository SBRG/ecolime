from collections import defaultdict
import pandas as pd
from sympy import Basic, Symbol
from os.path import dirname, join, abspath

import cobra.test
from cobrame import dogma

ecoli_data_files_dir = dirname(abspath(__file__))

del dirname, abspath


def fixpath(filename):
    return join(ecoli_data_files_dir, filename)


def get_biomass_composition(model, solution=None):
    if solution is None:
        solution = model.solution
    biomass_composition = defaultdict(float)
    for met, stoich in model._protein_biomass_dilution.metabolites.items():
        if abs(stoich) > 1:
            protein_stoich = stoich
    biomass_composition['Protein'] = \
        solution.x_dict['protein_biomass_dilution'] * protein_stoich
    biomass_composition['tRNA'] = \
        solution.x_dict['tRNA_biomass_dilution']
    biomass_composition['mRNA'] = \
        solution.x_dict['mRNA_biomass_dilution']
    biomass_composition['ncRNA'] = \
        solution.x_dict['ncRNA_biomass_dilution']
    biomass_composition['rRNA'] = \
        solution.x_dict['rRNA_biomass_dilution']
    biomass_composition['Other'] = \
        solution.x_dict['biomass_component_dilution']

    return biomass_composition


def RNA_to_protein_ratio(model, solution=None):
    composition = get_biomass_composition(solution=solution)
    RNA_to_protein = (composition['mRNA'] + composition['tRNA'] +
                      composition['rRNA'] + composition['ncRNA']) / \
                     (composition['Protein'] +
                      composition['Unmodeled Protein'])
    return RNA_to_protein


def get_RNA_fractions_dict(model, solution=None):
    RNA_fractions = {}
    composition = get_biomass_composition(solution=solution)

    tRNA_to_RNA = (composition['tRNA']) / (
        composition['mRNA'] + composition['tRNA'] + composition['rRNA'] +
        composition['ncRNA'])
    RNA_fractions['tRNA'] = tRNA_to_RNA

    rRNA_to_RNA = (composition['rRNA']) / (
        composition['mRNA'] + composition['tRNA'] + composition['rRNA'] +
        composition['ncRNA'])
    RNA_fractions['rRNA'] = rRNA_to_RNA

    mRNA_to_RNA = (composition['mRNA']) / (
        composition['mRNA'] + composition['tRNA'] + composition['rRNA'] +
        composition['ncRNA'])
    RNA_fractions['mRNA'] = mRNA_to_RNA

    ncRNA_to_RNA = (composition['ncRNA']) / (
        composition['mRNA'] + composition['tRNA'] + composition['rRNA'] +
        composition['ncRNA'])
    RNA_fractions['ncRNA'] = ncRNA_to_RNA

    return RNA_fractions


def get_membrane_composition(model, membrane='Inner_Membrane',
                             solution=None):
    if not solution:
        solution = model.solution

    composition = {}

    composition['Dummy_protein'] = solution.x_dict[
        'SA_demand_dummy_protein_' + membrane]
    composition['Modeled_protein'] = solution.x_dict[
        'SA_demand_protein_' + membrane]
    composition['Lipid'] = solution.x_dict[
        'SA_demand_lipid_' + membrane]

    if membrane == 'Outer_Membrane':
        composition['LPS'] = solution.x_dict[
            'SA_demand_lps_' + membrane]
        composition['Lipoprotein'] = solution.x_dict[
            'SA_demand_lipoprotein_' + membrane]
    return composition


def get_membrane_protein(model, membrane='Inner_Membrane', solution=None):
    if not solution:
        solution = model.solution

    composition = defaultdict(float)

    membrane_SA = model.metabolites.get_by_id('SA_protein_' + membrane)
    for rxn in membrane_SA.reactions:
        if membrane_SA in rxn.products:
            coeff = rxn._metabolites[membrane_SA]
            flux = coeff * solution.x_dict[rxn.id]
            if flux != 0:
                composition[
                    rxn.id.replace('translocation_', '')] += flux
        elif 'lipid_modification' in rxn.id:
            coeff = rxn._metabolites[membrane_SA]
            flux = coeff * solution.x_dict[rxn.id]
            if flux != 0:
                composition[rxn.id.split('_')[0]] += flux
    composition['dummy'] = solution.x_dict[
        'SA_demand_dummy_protein_' + membrane]
    return composition


def make_composition_piechart(model, kind='Biomass', solution=None):
    try:
        import brewer2mpl
    except ImportError:
        color_map = None
    else:
        color_map = brewer2mpl.wesanderson.Zissou.mpl_colormap

    try:
        import pandas
    except ImportError:
        raise Exception("Pandas must be installed to get biomass piechart")

    if solution is None:
        solution = model.solution

    summary = {}
    if kind == 'Biomass':
        summary['Biomass composition'] = \
            get_biomass_composition(solution=solution)
        frame = pandas.DataFrame.from_dict(summary) / solution.f

    elif kind == 'Inner_Membrane':
        SA = solution.x_dict['SA_demand_Inner_Membrane']
        summary['Inner_Membrane'] = \
            model.get_membrane_composition(solution=solution, membrane=kind)
        frame = pandas.DataFrame.from_dict(summary) / SA
    elif kind == 'Outer_Membrane':
        SA = solution.x_dict['SA_demand_Outer_Membrane']
        summary['Outer_Membrane'] = \
            model.get_membrane_composition(solution=solution, membrane=kind)
        frame = pandas.DataFrame.from_dict(summary) / SA
    elif kind == 'Protein_Inner_Membrane':
        SA = solution.x_dict['SA_demand_protein_Inner_Membrane'] + \
             solution.x_dict['SA_demand_dummy_protein_Inner_Membrane']
        summary['Inner_Membrane_Protein'] = \
            model.get_membrane_protein(solution=solution,
                                      membrane='Inner_Membrane')
        frame = pandas.DataFrame.from_dict(summary) / SA
    elif kind == 'Protein_Outer_Membrane':
        SA = solution.x_dict['SA_demand_protein_Outer_Membrane'] + \
             solution.x_dict['SA_demand_dummy_protein_Outer_Membrane']
        summary['Outer_Membrane_Protein'] = \
            model.get_membrane_protein(solution=solution,
                                      membrane='Outer_Membrane')
        frame = pandas.DataFrame.from_dict(summary) / SA

    else:
        raise ('%s not a valid composition kind' % kind)

    print('Component sum =', frame.sum().values[0])
    frame = frame[frame > 1e-4].dropna()
    return frame.plot(kind='pie', subplots=True, legend=None,
                      colormap=color_map)


def compare_to_ijo_biomass(model, kind='amino_acid', solution=None):
    ijo = cobra.test.create_test_model('ecoli')
    biomass_rxn = ijo.reactions.Ec_biomass_iJO1366_core_53p95M
    me_demand = defaultdict(float)

    # These are reactions that incorporate metabolites into biomass
    skip_list = ['SummaryVariable', 'ComplexFormation',
                 'TranscriptionReaction', 'TranslationReaction']

    growth_rate = model.solution.x_dict['biomass_dilution']
    mu = Symbol('mu')

    if solution:
        model.solution = solution

    if kind == 'amino_acid':
        for r in model.reactions.query('translation_'):
            for aa_letter, aa in dogma.amino_acids.items():
                me_demand[aa] += r.translation_data.amino_acid_count[aa] * r.x
    elif kind == 'cofactors':
        for d in model.complex_data:
            for mod, num in d.modifications.items():
                me_demand[mod.replace('mod_', '')] += d.formation.x * num
    else:
        for met_id in biomass_rxn.metabolites:
            met = model.metabolites.get_by_id(met_id.id)
            for r in met.reactions:
                if r.__class__.__name__ not in skip_list:
                    stoich = r._metabolites[met]
                    if isinstance(stoich, Basic):
                        stoich = stoich.subs(mu, growth_rate)
                    me_demand[met_id.id] += r.x * stoich


    compare ={}
    compare['ME_gr_%.2f' % growth_rate] = me_demand.copy()
    compare['Measured'] = {}
    for met in me_demand:

        if met in ijo.metabolites:
            ijo_met = ijo.metabolites.get_by_id(met)
        else:
            continue

        if ijo_met in biomass_rxn.metabolites and biomass_rxn._metabolites[ijo_met] < 0:
            compare['Measured'][met] = \
                abs(model.solution.x_dict['biomass_dilution'] *
                    biomass_rxn._metabolites[ijo_met])

    return pd.DataFrame.from_dict(compare).dropna(how='any')


def get_protein_distribution(model, solution=None, groupby='COG'):
    """
    Return the synthesis flux for each protein.

    First implementation-can be improved
    """
    protein_dict = defaultdict(float)
    if not solution:
        solution = model.solution

    transcription = model.get_translation_flux(solution=solution)

    if groupby == 'COG':
        COG_df = pd.read_csv(fixpath('data/cogs_ecoli_mg1655.csv'))
        COG_df = COG_df.set_index('locus')
        for protein, flux in transcription.items():
            protein_mass = model.translation_data.get_by_id(protein).mass
            if protein != 'dummy':
                if protein in COG_df. index:
                    try:
                        COG = COG_df.loc[protein, 'COG description'].values[0]
                    except:
                        COG = COG_df.loc[protein, 'COG description']
                else:
                    COG = 'No COG'

                protein_dict[COG] += protein_mass * flux
            else:
                protein_dict['dummy'] += protein_mass * flux

    elif groupby == 'Metabolic_Subsystem':
        for protein_id, flux in transcription.items():
            if flux <= 0:
                continue
            subsystem = set()
            protein = model.metabolites.get_by_id('protein_' + protein_id)
            protein_mass = protein.mass
            for complex in protein.complexes:
                for metabolic_reaction in complex.metabolic_reactions:
                    subsystem.add(metabolic_reaction.subsystem)
            if len(subsystem) > 0:
                protein_dict['_'.join(subsystem)] = protein_mass * flux
    else:
        raise('groupby flag is not "COG" or "Metabolic_Subsystem')

    return protein_dict
