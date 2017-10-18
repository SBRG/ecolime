# coding: utf-8

# # Build *i*LE1678-ME

# We will build an ME model from the NC_000913.2 Genbank file, the iJO1366
# M model, and the complex reconstruction from iJL1678-ME

# In[ ]:

# python imports

import re
from os.path import join
from collections import defaultdict
import pickle

# third party imports
import pandas
import cobra

# ECOLIme
import ecolime
from ecolime import (transcription, translation, flat_files, generics,
                     formulas)

# COBRAme
import cobrame
from cobrame.util import building, mu, me_model_interface
from cobrame.util.mass import dna_mw_no_ppi

# ## Part 1: Create minimum solveable ME-model
# This will include the bare minimum representations of
# - Transcription Reactions
# - Translation Reactions
# - Complex Formation Reactions
# - Metabolic Reactions
#
# that still produce a working ME-model
#
# ### 1) Create Model Object and populate its global info
# This includes important parameters that are used to calculate coupling
# constraints as well as organism-specific information such as peptide
# processing types

# In[ ]:


def return_me_model():
    # Define Models
    ijo_directory = join(flat_files.ecoli_files_dir, 'iJO1366.json')
    ijo = cobra.io.load_json_model(ijo_directory)
    me = cobrame.MEModel('iJL1678b-ME')

    # ME-models require special OptLang interface if cobrapy version >= 0.6.0
    # If cannot import SymbolicParameter, assume using cobrapy
    # versions <= 0.5.11
    try:
        from optlang.interface import SymbolicParameter
    except:
        pass
    else:
        me.solver = me_model_interface

    # "Translational capacity" of organism
    me.global_info['kt'] = 4.5  # (in h-1)scott 2010, RNA-to-protein curve fit
    me.global_info['r0'] = 0.087  # scott 2010, RNA-to-protein curve fit
    me.global_info['k_deg'] = 1.0 / 5. * 60.0  # 1/5 1/min 60 min/h # h-1

    # Molecular mass of RNA component of ribosome
    me.global_info['m_rr'] = 1453.  # in kDa

    # Average molecular mass of an amino acid
    me.global_info['m_aa'] = 109. / 1000.  # in kDa

    # Proportion of RNA that is rRNA
    me.global_info['f_rRNA'] = .86
    me.global_info['m_nt'] = 324. / 1000.  # in kDa
    me.global_info['f_mRNA'] = .02

    # tRNA associated global information
    me.global_info['m_tRNA'] = 25000. / 1000.  # in kDA
    me.global_info['f_tRNA'] = .12

    # Folding Properties
    me.global_info['temperature'] = 37
    me.global_info['propensity_scaling'] = .45

    # DNA Replication Parameters
    me.global_info['GC_fraction'] = 0.507896997096

    # Define the types of biomass that will be synthesized in the model
    me.add_biomass_constraints_to_model(
        ["protein_biomass", "mRNA_biomass", "tRNA_biomass", "rRNA_biomass",
         "ncRNA_biomass", "DNA_biomass", "lipid_biomass",
         "constituent_biomass", "prosthetic_group_biomass",
         "peptidoglycan_biomass"])

    # ### 2) Load metabolites and build Metabolic reactions
    # The below reads in:
    # - Required
    #  * **reaction_matrix.txt** (reaction matrix w/ reactions unlumped, metabolites renamed etc.)
    #  * **metabolites.txt** (metabolite properties)
    #  * **reactions.txt** (info on reversiblity, whether enzyme catalyzed etc.)
    #  * **m_to_me_mets.csv** (mapping of enzymes/complexes used in M-model to their ME-model compatible ID)
    #  * **protein_complexes.txt** (protein subunit stoichiometry of all complexes, used to identify metabolites as such)
    #
    # It creates a new e coli M-model from this info then incorporates it into the ME-model using *add_m_model_content*. metabolites are added directly reactions are added as StoichiometricData
    #
    # Metabolite types have different properties in an ME-model so enzyme complexes need added to the model as Complexes not Metabolites. Components in the E. coli M-model that are actually Complexes are compiled in *complex_list*

    # In[ ]:

    # m_model = flat_files.get_m_model()
    m_model = flat_files.process_m_model(ijo, 'metabolites.txt',
                                         'm_to_me_mets.csv', 'reactions.txt',
                                         'reaction_matrix.txt',
                                         'protein_complexes.txt',
                                         defer_to_rxn_matrix={'GLUTRR', 'PAPSR2'})
    m_model.reactions.EX_glc_e.id = 'EX_glc__D_e'
    m_model.repair()
    # some of the "metabolites" in iJO1366 "M" model are actually complexes. We pass those in
    # so they get created as complexes, not metabolites.
    complexes = flat_files.get_complex_subunit_stoichiometry(
        'protein_complexes.txt').keys()
    complex_set = set(
        [i.id for i in m_model.metabolites if i.id.split('_mod_')[0] in complexes])
    building.add_m_model_content(me, m_model, complex_metabolite_ids=complex_set)

    # In[ ]:

    # This adds exchange reactions for metabolites not contained in iJO1366
    # Some of these cannot be produced by the model so they are added here
    exchange_list = [
        'LI_c',
        'pqq_e',
        'cs_e',
        'tl_c',
        'RNase_m5', 'RNase_m16', 'RNase_m23']  # RNAses are gaps in model

    for met_id in exchange_list:
        r = cobrame.MEReaction("EX_" + met_id)
        me.add_reaction(r)
        r.reaction = met_id + " <=> "

    # ### 3) Add Transcription and Translation
    # The below reads in:
    # - Required
    #  * **NC_000913.2.gb** (Genbank sequence annotation)
    #  * **ecolime/translation.py** (codon to tRNA mapping)
    # - Optional
    #  * **TUs_from_ecocyc.txt** (TU definitions, start/stop positions, etc.)
    #  * **ecolime/translation.py** (dictionary of gene to frameshift mutation)
    #
    # To construct the bare minimimum components of a transcription and translation reactions. For example, transcription reactions at this point include nucleotides and the synthesized RNAs.

    # In[ ]:

    gb_filename = join(flat_files.ecoli_files_dir, 'NC_000913.2.gb')
    tu_df = flat_files.get_tu_dataframe('TUs_from_ecocyc.txt')

    building.build_reactions_from_genbank(me, gb_filename, tu_df, verbose=False,
                                          frameshift_dict=translation.frameshift_dict,
                                          trna_to_codon=translation.trna_to_codon)

    # ### 4) Add in complex Formation without modifications (for now)
    #
    # The below reads in:
    # - Required
    #  * **protein_complexes.txt** (Metabolic complexes' protein subunit stoichiometries)
    #  * **protein_modification.txt** (Type and number of modifications for each protein)

    # In[ ]:

    # complex_stoichiometry_dict is a dict of {'complex_id': [{'bnum' : count}]}
    rna_components = {"b3123",
                      "b0455"}  # component id should have 'RNA_ instead' of 'protein_'
    # get Metabolic Complex composition from ECOLIme
    complex_stoichiometry_dict = flat_files.get_complex_subunit_stoichiometry(
        'protein_complexes.txt', rna_components)
    # add complexes to model
    complex_modification_dict = flat_files.get_complex_modifications(
        'protein_modification.txt', 'protein_complexes.txt')
    building.add_model_complexes(me, complex_stoichiometry_dict,
                                 complex_modification_dict, verbose=False)

    # remove modifications. they will be added back in later
    for data in me.complex_data:
        data.subreactions = {}

    # add formation reactions for each of the ComplexDatas
    for cplx_data in me.complex_data:
        formation = cplx_data.formation
        if formation:
            formation.update()
        else:
            cplx_data.create_complex_formation()

    # ### 5) Add dummy reaction to model and unmodeled_protein_fraction
    #
    # Includes the transcription, translation, complex_formation, and metabolic reaction. Sequence based on prevelance of each codon found in *E. coli*.
    #  - Required
    #      * [**codon_usage.csv**](http://openwetware.org/wiki/Escherichia_coli/Codon_usage) (codon prevelance)

    # In[ ]:

    seq = "ATG"
    codons = pandas.read_csv(join(flat_files.ecoli_files_dir, "codon_usage.csv"),
                             index_col=0)
    for codon, row in codons.iterrows():
        if row.amino_acid == "Stop":
            continue
        seq += codon * int(row.per_1000 // 3)  # want roughly 300 aa
    # get the most used stop codon
    seq += codons[codons.amino_acid == "Stop"].sort_values("per_1000").index[-1]
    building.add_dummy_reactions(me, seq, update=True)

    rxn = cobrame.SummaryVariable('dummy_protein_to_mass')
    me.add_reactions([rxn])
    mass = me.metabolites.protein_dummy.mass
    rxn.add_metabolites({'protein_biomass': -mass, 'protein_dummy': -1,
                         cobrame.Constraint('unmodeled_protein_biomass'): mass})

    # ### 6) Assocated Complexes and build Metabolic Reactions
    # - Required
    #     * **enzyme_reaction_association.txt**
    #     * **reactions.txt** (gives reaction name, reversiblity, source and whether reaction is spontaneous)
    #

    # In[ ]:

    # associate reaction id with the old ME complex id (including modifications)
    rxn_to_cplx_dict = flat_files.get_reaction_to_complex(m_model)
    rxn_info = flat_files.get_reaction_info_frame('reactions.txt')

    # Required to add dummy reaction as spontaneous reaction
    rxn_info = rxn_info.append(
        pandas.Series({'description': 'dummy reaction', 'is_reversible': 0,
                       'is_spontaneous': 1}, name='dummy_reaction'))

    building.add_reactions_from_stoichiometric_data(me, rxn_to_cplx_dict, rxn_info,
                                                    update=True)

    # ### 7) Incorporate remaining biomass constituents
    # There are leftover components from the *i*JO1366 biomass equation that either:
    # 1. have no mechanistic function in the model (glycogen)
    # 2. are cofactors that are regenerated (nad)
    #
    # Applies demands and coefficients from the *i*JO1366 biomass objective function

    # In[ ]:

    me.ngam = 1
    me.gam = 34.98
    me.unmodeled_protein_fraction = .36

    biomass_constituents = {
        "glycogen_c": -.023 / (me.metabolites.glycogen_c.formula_weight / 1000.),
        "2ohph_c": -0.000223,
        "nad_c": -.001831,
        "udcpdp_c": -5.5e-05,
        "coa_c": -0.000576,
        "ribflv_c": -0.000223,
        "nadp_c": -0.000447,
        "mlthf_c": -0.000223,
        "thf_c": -0.000223,
        "10fthf_c": -0.000223
    }

    rxn = cobrame.SummaryVariable('biomass_constituent_demand')
    me.add_reactions([rxn])
    rxn.add_metabolites(biomass_constituents)
    constituent_mass = sum(
        me.metabolites.get_by_id(c).formula_weight / 1000. * -v
        for c, v in biomass_constituents.items())
    rxn.lower_bound = mu
    rxn.upper_bound = mu
    rxn.add_metabolites({me.metabolites.constituent_biomass: constituent_mass})

    #  #### Lipid components
    #  Metabolites and coefficients from *i*JO1366 biomass objective function

    # In[ ]:

    # Find biomass constituents with 3 numbers followed by a compartment in the BOF
    lipid = re.compile('\d{3}_.')
    lipid_demand = {}
    for key, value in ijo.reactions.Ec_biomass_iJO1366_core_53p95M.metabolites.items():
        if lipid.search(key.id):
            lipid_demand[key.id] = abs(value)

    for met, requirement in lipid_demand.items():
        component_mass = me.metabolites.get_by_id(met).formula_weight / 1000.
        rxn = cobrame.SummaryVariable('Demand_' + met)
        me.add_reactions([rxn])
        rxn.add_metabolites({met: -1 * requirement,
                             'lipid_biomass': component_mass * requirement})
        rxn.lower_bound = mu
        rxn.upper_bound = 1000.

    # Kdo2lipid4
    requirement = 0.01945  # in mmol/gDW
    met = me.metabolites.get_by_id('kdo2lipid4_e')
    component_mass = met.formula_weight / 1000.
    rxn = cobrame.SummaryVariable('Demand_' + met.id)
    me.add_reactions([rxn])

    rxn.add_metabolites({met.id: -1. * requirement,
                         'lipid_biomass': component_mass * requirement})
    rxn.lower_bound = mu
    rxn.upper_bound = mu

    # #### DNA Demand Requirements
    # Added based on growth rate dependent DNA levels as in [O'brien EJ et al 2013](https://www.ncbi.nlm.nih.gov/pubmed/24084808)

    # In[ ]:

    dna_demand_stoich, dna_demand_bound = ecolime.dna_replication.return_gr_dependent_dna_demand(
        me.global_info['GC_fraction'])

    dna_replication = cobrame.SummaryVariable("DNA_replication")
    me.add_reaction(dna_replication)
    dna_replication.add_metabolites(dna_demand_stoich)
    dna_biomass = cobrame.Constraint("DNA_biomass")
    dna_biomass.elements = {e: abs(v) for e, v in
                            dna_replication.check_mass_balance().items()}

    dna_mw = 0
    for met, value in me.reactions.DNA_replication.metabolites.items():
        if met.id != 'ppi_c':
            dna_mw -= value * dna_mw_no_ppi[met.id.replace('_c', '')] / 1000.

    dna_replication.add_metabolites({dna_biomass: dna_mw})
    dna_replication.lower_bound = dna_demand_bound
    dna_replication.upper_bound = dna_demand_bound

    # **Note**: From this point forward, executing every codeblock should result in a solveable ME-model
    #
    # ------
    #
    # ## Part 2: Add metastructures to solving ME-model
    # This includes:
    # 1. ribosome
    # 2. RNA polymerase
    # 3. charged_tRNAs

    # Sometimes multiple entities can perform the same role. To prevent a combinatorial explosion of possibilities, we can create  "generic" version, where any of those entities can fill in.

    # In[ ]:

    for generic, components in generics.generic_dict.items():
        cobrame.GenericData(generic, me, components).create_reactions()

    # ### 1) Add ribosome
    # This uses the ribosome composition definition in **ecolime/ribosome.py**

    # In[ ]:

    ecolime.ribosome.add_ribosome(me, verbose=False)

    # ### 2) Add charged tRNA reactions

    # The tRNA charging reactions were automatically added when loading the genome from the genbank file. However, the charging reactions still need to be made aware of the tRNA synthetases which are responsible.
    #
    # Uses **trna_charging.py**

    # In[ ]:

    aa_synthetase_dict = ecolime.trna_charging.amino_acid_trna_synthetase
    for data in me.tRNA_data:
        data.synthetase = str(aa_synthetase_dict[data.amino_acid])

    # Generic charged tRNAs are added to translation reactions via SubreactionData below.
    #
    # All new data added in this block contained in **ecolime/translation.py**

    # In[ ]:

    ecolime.translation.add_charged_trna_subreactions(me)
    for data in me.translation_data:
        data.add_initiation_subreactions(
            start_codons=translation.translation_start_codons,
            start_subreactions=set(translation.initiation_subreactions.keys()))

        data.add_elongation_subreactions(elongation_subreactions=set(
            translation.elongation_subreactions.keys()))

        data.add_termination_subreactions(
            translation_terminator_dict=translation.translation_stop_dict)

    # ### 3) Add Transcription Metacomplexes
    # #### RNA Polymerase
    #
    # Data for RNA_polymerase composition fround in **ecolime/transcription**
    #
    # Uses *tu_df* from **TUs_from_ecocyc.txt**, above

    # In[ ]:

    for met in transcription.rna_polymerases:
        rnap_obj = cobrame.RNAP(met)
        me.add_metabolites(rnap_obj)
    transcription.add_rna_polymerase_complexes(me, verbose=False)

    # associate the correct RNA_polymerase and factors to TUs
    sigma_to_rnap_dict = transcription.sigma_factor_complex_to_rna_polymerase_dict
    for tu_id in tu_df.index:
        transcription_data = me.process_data.get_by_id(tu_id)
        sigma = tu_df.sigma[tu_id]
        rna_polymerase = sigma_to_rnap_dict[sigma]
        transcription_data.RNA_polymerase = rna_polymerase

    # #### Degradosome (both for RNA degradation and RNA splicing)
    #
    # All new data contained in **transcription.py**

    # In[ ]:

    me.add_metabolites([cobrame.Complex('RNA_degradosome')])
    data = cobrame.ComplexData('RNA_degradosome', me)
    for subunit, value in transcription.rna_degradosome.items():
        data.stoichiometry[subunit] = value
    data.create_complex_formation(verbose=False)

    # Used for RNA splicing
    data = cobrame.SubreactionData('RNA_degradation_machine', me)
    data.enzyme = 'RNA_degradosome'

    data = cobrame.SubreactionData('RNA_degradation_atp_requirement', me)
    # .25 water equivaltent for atp hydrolysis per nucleotide
    data.stoichiometry = {'atp_c': -.25, 'h2o_c': -.25, 'adp_c': .25,
                          'pi_c': .25, 'h_c': .25}

    transcription.add_rna_splicing(me)

    # ------
    # ## Part 3: Add remaining modifications
    # rRNA modifications handled in *add_ribosome*
    #
    # ### 1) Add complex modifications
    # *complex_modification_dict* from **protein_modification.text**, above
    #
    # The rest of the new data contained in **ecolime/modifications.py**

    # In[ ]:

    for complex_id, info in complex_modification_dict.items():
        modifications = {}
        for mod, value in info['modifications'].items():
            # stoichiometry of modification determined in
            # subreaction_data.stoichiometry
            modifications['mod_' + mod] = abs(value)
        me.process_data.get_by_id(complex_id).subreactions = modifications

    # Adds modification data for more complicated enzyme modifications
    # (ie, iron sulfur cluster modification)
    ecolime.modifications.add_modification_procedures(me)

    # add formation reactions for each of the ComplexDatas
    for cplx_data in me.complex_data:
        formation = cplx_data.formation
        if formation:
            formation.update()
        else:
            cplx_data.create_complex_formation()

    # ### 2) Add tRNA mods and asocciate them with tRNA charging reactions
    # New data from:
    # 1. **ecolime/trna_charging.py** (read via *add_trna_modification_procedures()*)
    # 2. **post_transcriptional_modification_of_tRNA.txt** (modification types per tRNA)
    #

    # In[ ]:

    # Add tRNA modifications to ME-model
    ecolime.trna_charging.add_trna_modification_procedures(me)

    # trna_modifications = {tRNA_id: {modifications: count}}
    trna_modifications = flat_files.get_trna_modification_targets()
    for trna in trna_modifications:
        for data in me.process_data.query(trna):
            data.subreactions = trna_modifications[trna]

    # ---
    # ## Part 4: Add remaining subreactions
    # ### 1) Add translation related subreactions
    # All new data from **ecolime/translation.py**

    # In[ ]:

    # add the translation subreaction data objects to model
    translation.add_translation_subreactions_to_model(me)

    # add translation subreaction data to reactions
    methionine_cleaved = translation.methionine_cleaved
    folding_dict = translation.folding_dict

    for data in me.translation_data:
        data.term_enzyme = \
            translation.translation_stop_dict.get(data.last_codon)

        locus_id = data.id
        if locus_id in methionine_cleaved:
            data.subreactions['N_terminal_methionine_cleavage'] = 1

        for folding_type in folding_dict:
            if locus_id in folding_dict[folding_type]:
                data.subreactions[folding_type] = 1

        # This block was ran above, but should be ran again to
        # incorporate any subreactions not added previously
        data.add_initiation_subreactions(
            start_codons=translation.translation_start_codons,
            start_subreactions=set(translation.initiation_subreactions.keys()))
        data.add_elongation_subreactions(elongation_subreactions=set(
            translation.elongation_subreactions.keys()))
        data.add_termination_subreactions(
            translation_terminator_dict=translation.translation_stop_dict)

        # add organism specific subreactions associated with peptide processing
        for subrxn in translation.peptide_processing_subreactions:
            data.subreactions[subrxn] = 1

    # ### 2) Add transcription related subreactions
    # All new data from **ecolime/transcription.py**

    # In[ ]:

    for subreaction in transcription.transcription_subreactions:
        subreaction_data = cobrame.SubreactionData(subreaction, me)
        enzymes = transcription.transcription_subreactions[subreaction]['enzymes']
        subreaction_data.stoichiometry = \
        transcription.transcription_subreactions[subreaction]['stoich']
        subreaction_data.enzyme = enzymes

    for transcription_data in me.transcription_data:
        # Assume false if not in tu_df\n",
        rho_dependent = tu_df.rho_dependent.get(transcription_data.id, False)
        rho = 'dependent' if rho_dependent else 'independent'
        stable = 'stable' if transcription_data.codes_stable_rna else 'normal'
        transcription_data.subreactions['Transcription_%s_rho_%s' % (stable,
                                                                     rho)] = 1

    # ----
    # ## Part 5: Add in translocation
    #
    # New data from:
    # 1. **peptide_compartment_and_pathways.txt** (Protein compartment and translocation pathway for each membrane complex)
    # 2. **ecolime/translocation.py** (definitions of each translocation pathway)

    # In[ ]:

    # Add TranslocationData
    transloc = pandas.read_csv(
        join(flat_files.ecoli_files_dir, "peptide_compartment_and_pathways2.txt"),
        sep='\t', comment="#")
    for pathway, info in ecolime.translocation.pathway.items():
        if 'alt' not in pathway:
            transloc_data = cobrame.TranslocationData(pathway + '_translocation',
                                                      me)
        else:
            transloc_data = cobrame.TranslocationData(
                pathway.replace('_alt', '_translocation_alt'), me)
        transloc_data.enzyme_dict = info['enzymes']
        transloc_data.keff = info['keff']
        transloc_data.length_dependent_energy = info['length_dependent_energy']
        transloc_data.stoichiometry = info['stoichiometry']

    # Associate data and add translocation reactions
    ecolime.translocation.add_translocation_pathways(me, transloc,
                                                     membrane_constraints=False)

    # Update stoichiometry of membrane complexes
    # new_stoich = {complex_id: protein_w_compartment}
    new_stoich = defaultdict(dict)
    for cplx, row in transloc.set_index('Complex').iterrows():
        if cplx == 'EG10544-MONOMER':
            continue
        protein = row.Protein.split('(')[0] + '_' + row.Protein_compartment
        value = row.Protein.split('(')[1][:-1].split(':')[0]
        new_stoich[cplx]['protein_' + protein] = float(value)

    for cplx, stoich in new_stoich.items():
        complex_data = me.process_data.get_by_id(cplx)
        for met, value in stoich.items():
            complex_data.stoichiometry.pop(met[0:13])
            complex_data.stoichiometry[met] = value
            complex_data.formation.update()
        # Complex ids in protein compartment file doesn't include mods
        # Some have multiple alternative modifications so must loop through these
        for complex_data in me.process_data.query(cplx + '_mod_'):
            for met, value in stoich.items():
                complex_data.stoichiometry.pop(met[0:13])
                complex_data.stoichiometry[met] = value
                complex_data.formation.update()

    # ---
    # ## Part 6: Add Cell Wall Components
    # All new data from **ecolime/translocation.py**

    # In[ ]:

    compartment_dict = {}
    for prot, compartment in transloc.set_index(
            'Protein').Protein_compartment.to_dict().items():
        compartment_dict[prot.split('(')[0]] = compartment

    # #### Add lipid modification SubreactionData

    # In[ ]:

    lipid_modifications = ecolime.translocation.lipid_modifications

    for lipid in lipid_modifications:
        data = cobrame.SubreactionData('mod_' + lipid, me)
        data.stoichiometry = {lipid: -1, 'g3p_c': 1}
        data.enzyme = ['Lgt_MONOMER', 'LspA_MONOMER']
        # The element contribution is based on the lipid involved in the
        # modification, so calculate based on the metabolite formula
        data._element_contribution = data.calculate_element_contribution()

    data = cobrame.SubreactionData('mod2_pg160_p', me)
    data.stoichiometry = {'pg160_p': -1, '2agpg160_p': 1}
    data.enzyme = 'EG10168-MONOMER'
    data._element_contribution = data.calculate_element_contribution()

    data = cobrame.SubreactionData('mod2_pe160_p', me)
    data.stoichiometry = {'pe160_p': -1, '2agpe160_p': 1}
    data.enzyme = 'EG10168-MONOMER'
    data._element_contribution = data.calculate_element_contribution()

    ecolime.translocation.add_lipoprotein_formation(me, compartment_dict,
                                                    membrane_constraints=False)

    # #### Correct complex formation IDs if they contain lipoproteins

    # In[ ]:

    for gene in ecolime.translocation.lipoprotein_precursors.values():
        compartment = compartment_dict.get(gene)
        for rxn in me.metabolites.get_by_id(
                                        'protein_' + gene + '_' + compartment).reactions:
            if isinstance(rxn, cobrame.ComplexFormation):
                data = me.process_data.get_by_id(rxn.complex_data_id)
                value = data.stoichiometry.pop(
                    'protein_' + gene + '_' + compartment)
                data.stoichiometry[
                    'protein_' + gene + '_lipoprotein' + '_' + compartment] = value
                rxn.update()

    # #### Braun's lipoprotein demand
    # Metabolites and coefficients as defined in [Liu et al 2014](http://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-014-0110-6)

    # In[ ]:

    rxn = cobrame.SummaryVariable('core_structural_demand_brauns')
    met1 = me.metabolites.get_by_id('murein5px4p_p')
    met1_mass = met1.formula_weight / 1000.
    met2 = me.metabolites.get_by_id('protein_b1677_lipoprotein_Outer_Membrane')
    me.add_reactions([rxn])
    # biomass of lipoprotein accounted for in translation and lipip_modification
    rxn.add_metabolites({met1: -0.013894, met2: -0.003597,
                         'peptidoglycan_biomass': (0.013894 * met1_mass)},
                         combine=False)
    rxn.lower_bound = mu
    rxn.upper_bound = mu

    # -----
    # ## Part 7: Set keffs
    #
    # Either entirely based on SASA or using fit keffs from [Ebrahim et al 2016](https://www.ncbi.nlm.nih.gov/pubmed/27782110?dopt=Abstract)
    # # Set keffs to sasa fluxes centered around 65.
    # me.set_SASA_keffs(65)
    # In[ ]:

    keff_list = []
    keffs = flat_files.get_reaction_keffs(me, verbose=True)
    for reaction_id, keff in keffs.items():
        if keff > 3000:
            keff = 3000.
        elif keff < .01:
            keff = .01
        keff_list.append(keff)
        me.reactions.get_by_id(reaction_id).keff = keff
        me.reactions.get_by_id(reaction_id).update()

    # Keffs that were not set in the above block
    me.process_data.N_terminal_methionine_cleavage.keff = 1339.4233102860871
    me.process_data.peptide_deformylase_processing.keff = 1019.5963333345715
    me.reactions.get_by_id(
        'GLUTRR_FWD_CPLX0-3741').keff = 3000  # 3269.0108007383374
    me.process_data.fmet_addition_at_START.keff = 1540.4356849968603
    me.process_data.ribosome_recycler.keff = 1059.6910912619182
    me.process_data.UAG_PrfA_mono_mediated_termination.keff = 1721.7910609284945
    me.process_data.UGA_PrfB_mono_mediated_termination.keff = 1700.2966587695353
    me.process_data.UAA_generic_RF_mediated_termination.keff = 1753.4238515034572

    # -----
    # ## Part 8: Model updates and corrections
    # # Add NDH flux split constraint
    # for rxn in me.metabolites.get_by_id('NADH-DHII-MONOMER_mod_mg2_mod_cu_mod_fad').metabolic_reactions:
    #     rxn.stoichiometric_data.stoichiometry['ndh2_constraint'] = 1
    #     rxn.update()
    # for rxn in me.metabolites.get_by_id('NADH-DHI-CPLX_mod_2fe2s_mod_4fe4s_mod_fmn').metabolic_reactions:
    #     rxn.stoichiometric_data.stoichiometry['ndh1_constraint'] = 1
    #     rxn.update()
    # rxn = cobra.Reaction('ndh_flux_split_constraint')
    # me.add_reaction(rxn)
    # rxn.reaction = 'ndh1_constraint + ndh2_constraint ->'
    # In[ ]:

    # Add reaction subsystems from iJO to model
    for rxn in ijo.reactions:
        if rxn.id in me.process_data:
            data = me.process_data.get_by_id(rxn.id)
        else:
            continue
        for r in data.parent_reactions:
            r.subsystem = rxn.subsystem

    # #### Corrections and final updates

    # In[ ]:

    ecolime.corrections.correct_reaction_stoichiometries(me, join(
        flat_files.ecoli_files_dir,
        'iLE1678_model_changes.xlsx'))
    # RNA_dummy, TU_b3247, TU_b3705 do not have RNAP, this is set as the most common RNAP
    for data in me.transcription_data:
        if len(data.RNA_polymerase) == 0:
            data.RNA_polymerase = 'RNAP70-CPLX'

    # If lower_bound open, model feeds G6P into EDD
    me.reactions.EX_pqq_e.lower_bound = 0
    me.reactions.EX_pqq_e.upper_bound = 0

    # cobalamin is not in glucose M9 media
    me.reactions.EX_cbl1_e.lower_bound = 0

    me.process_data.PPKr.lower_bound = 0.
    me.process_data.PPKr.update_parent_reactions()

    # This enyzme is involved in catalyzing this reaction
    sub = cobrame.SubreactionData('EG12450-MONOMER_activity', me)
    sub.enzyme = 'EG12450-MONOMER'
    me.process_data.NHFRBO.subreactions['EG12450-MONOMER_activity'] = 1

    # #### Add enzymatic coupling for "carriers"
    # These are enzyme complexes that act as metabolites in a metabolic reaction (i.e. are metabolites in iJO1366)

    # In[ ]:

    for data in me.stoichiometric_data:
        if data.id == 'dummy_reaction':
            continue

        for met, value in data.stoichiometry.items():
            if not isinstance(me.metabolites.get_by_id(met),
                              cobrame.Complex) or value > 0:
                continue

            subreaction_id = met + '_carrier_activity'
            if subreaction_id not in me.process_data:
                sub = cobrame.SubreactionData(met + '_carrier_activity', me)
                sub.enzyme = met
            data.subreactions[subreaction_id] = abs(value)

    # ----
    # ## Part 9: Update and save

    # In[ ]:

    me.reactions.dummy_reaction_FWD_SPONT.objective_coefficient = 1.
    me.reactions.EX_glc__D_e.lower_bound = -1000
    me.reactions.EX_o2_e.lower_bound = -1000.
    me.ngam = 1.
    me.gam = 34.98
    me.unmodeled_protein_fraction = .36

    # In[ ]:

    me.update()
    me.prune()

    # ### Add remaining metabolite formulas to model

    # In[ ]:

    # Update a second time to incorporate all of the metabolite formulas corectly
    for r in me.reactions.query('formation_'):
        r.update()

    # Update complex formulas
    df = \
        pandas.read_table(join(flat_files.ecoli_files_dir, 'modification.txt'),
                          names=['mod', 'formula', 'na'])
    df = df.drop('na', axis=1).set_index('mod').dropna(how='any')
    modification_formulas = df.T.to_dict()
    formulas.add_remaining_complex_formulas(me, modification_formulas)

    me.metabolites.get_by_id(
        'CPLX0-782_mod_1:2fe2s_mod_1:4fe4s').formula = 'C3164Fe6H5090N920O920S50'
    me.metabolites.get_by_id(
        'EG50003-MONOMER_mod_pan4p_mod_lipo').formula = 'C387H606N95O142PS4'
    # Update reactions affected by formula update
    for r in me.reactions.query('_mod_lipo'):
        r.update()
    for r in me.reactions.query('_mod_glycl'):
        r.update()

    # In[ ]:

    n_genes = len(me.metabolites.query(re.compile('RNA_b[0-9]')))
    print("number of genes in the model %d (%.2f%%)" % (
    n_genes, n_genes * 100. / (1678)))

    return me

if __name__ == '__main__':

    me = return_me_model()
    with open("./me_models/iLE1678.pickle", "wb") as outfile:
        pickle.dump(me, outfile)
