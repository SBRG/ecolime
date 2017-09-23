from __future__ import division, absolute_import, print_function

from six import iteritems

from cobrame import ComplexData, ModificationData
from ecolime import generics

transcription_subreactions = {
    'Transcription_normal_rho_independent':
        {'enzymes': ['Mfd_mono_mod_1:mg2', 'NusA_mono', 'NusG_mono',
                     'GreA_mono', 'GreB_mono', 'RpoZ_mono_mod_1:mg2'],
         'stoich': {}},
    'Transcription_normal_rho_dependent':
        {'enzymes': ['Mfd_mono_mod_1:mg2', 'NusA_mono', 'NusG_mono',
                     'GreA_mono', 'GreB_mono', 'RpoZ_mono_mod_1:mg2',
                     'Rho_hexa_mod_3:mg2'],
         'stoich': {'atp_c': -3,
                    'h2o_c': -3,
                    'adp_c': 3,
                    'pi_c': 3,
                    'h_c': 3}},
    'Transcription_stable_rho_independent':
        {'enzymes': ['Mfd_mono_mod_1:mg2', 'NusA_mono', 'NusG_mono',
                     'GreA_mono', 'GreB_mono', 'RpoZ_mono_mod_1:mg2',
                     'RpsJ_mono',  'RpsD_mono', 'RplC_mono', 'RplD_mono',
                     'RplM_mono', 'NusB_mono'],
         'stoich': {}},
    'Transcription_stable_rho_dependent':
        {'enzymes': ['Mfd_mono_mod_1:mg2', 'NusA_mono', 'NusG_mono',
                     'GreA_mono', 'GreB_mono', 'RpoZ_mono_mod_1:mg2',
                     'Rho_hexa_mod_3:mg2',  'RpsJ_mono',  'RpsD_mono',
                     'RplC_mono', 'RplD_mono', 'RplM_mono', 'NusB_mono'],
         'stoich': {'atp_c': -3,
                    'h2o_c': -3,
                    'adp_c': 3,
                    'pi_c': 3,
                    'h_c': 3}}
}

sigma_factor_complex_to_rna_polymerase_dict = {
    # this is notation used in Tu_from_ecocyc, protein_complexes (and below)
    'RpoD_mono': 'RNAP70-CPLX',
    'RpoH_mono': 'RNAP32-CPLX',
    'RpoE_mono': 'RNAPE-CPLX',
    'RpoS_mono': 'RNAPS-CPLX',
    'RpoN_mono': 'RNAP54-CPLX',
    'FliA_mono': 'CPLX0-222',
    'FecI_mono': 'CPLX0-221'}

rna_polymerase_sigma_factor_components = {
    'CPLX0-221': {'sigma_factor': 'FecI_mono',
                  'polymerase': 'hRNAP_mod_1:zn2_mod_2:mg2'},
    'RNAPE-CPLX': {'sigma_factor': 'RpoE_mono',
                   'polymerase': 'hRNAP_mod_1:zn2_mod_2:mg2'},
    'CPLX0-222': {'sigma_factor': 'FliA_mono',
                  'polymerase': 'hRNAP_mod_1:zn2_mod_2:mg2'},
    'RNAP32-CPLX': {'sigma_factor': 'RpoH_mono',
                    'polymerase': 'hRNAP_mod_1:zn2_mod_2:mg2'},
    'RNAP54-CPLX': {'sigma_factor': 'RpoN_mono',
                    'polymerase': 'hRNAP_mod_1:zn2_mod_2:mg2'},
    'RNAP70-CPLX': {'sigma_factor': 'RpoD_mono',
                    'polymerase': 'hRNAP_mod_1:zn2_mod_2:mg2'},
    'RNAPS-CPLX': {'sigma_factor': 'RpoS_mono',
                   'polymerase': 'hRNAP_mod_1:zn2_mod_2:mg2'}}

RNA_polymerases = {'CPLX0-221', 'RNAPE-CPLX', 'CPLX0-222',
                   'RNAP32-CPLX', 'RNAP54-CPLX',
                   'RNAP70-CPLX', 'RNAPS-CPLX'}

RNA_degradosome = {'Eno_dim_mod_4:mg2': 1, 'Pnp_trim': 1,
                   'RNase_E_tetra_mod_2:zn2': 1, 'RhlB_dim': 1,
                   'Orn_dim_mod_2:mg2': 1}

excision_machinery = {
    'rRNA_containing': ['RNase_E_tetra_mod_2:zn2', 'RNase_P_cplx_mod_2:mg2',
                        'generic_RNase', 'RNase_m5', 'RNase_m16', 'RNase_m23',
                        'RNase_III_dim_mod_2:mg2', 'RNase_G_dim',
                        'RNase_T_dim_mod_4:mg2'],
    'monocistronic': ['RNase_E_tetra_mod_2:zn2', 'RNase_P_cplx_mod_2:mg2',
                      'generic_RNase'],
    'polycistronic_wout_rRNA': ['RNase_E_tetra_mod_2:zn2',
                                'RNase_P_cplx_mod_2:mg2', 'generic_RNase',
                                'RNase_III_dim', 'RNase_G_dim',
                                'RNase_T_dim_mod_4:mg2']}


def add_RNA_polymerase_complexes(me_model, verbose=True):

    for complex, components in iteritems(rna_polymerase_sigma_factor_components):
        rnap_complex = ComplexData(complex, me_model)
        rnap_components = rnap_complex.stoichiometry
        sigma_factor = components['sigma_factor']
        polymerase = components['polymerase']

        rnap_components[sigma_factor] = 1
        rnap_components[polymerase] = 1

        rnap_complex.create_complex_formation(verbose=verbose)


def add_RNA_splicing(me_model):

    # Ecoli has three alternatie mechanisms for splicing RNA, depending
    # on what RNA types the TU contains
    excision_types = ['rRNA_containing', 'monocistronic',
                      'polycistronic_wout_rRNA']

    for excision_type in excision_types:
        complex_data = ComplexData(excision_type + "_excision_machinery",
                                   me_model)

        for machine in excision_machinery[excision_type]:
            complex_data.stoichiometry[machine] = 1

        complex_data.create_complex_formation()
        modification = ModificationData(excision_type + "_excision", me_model)
        modification.enzyme = complex_data.id

    # Loop through transcription reactions and add appropriate splicing
    # machinery based on RNA types and number of splices required
    for t in me_model.transcription_data:
        n_excised = sum(t.excised_bases.values())
        n_cuts = len(t.RNA_products) * 2
        if n_excised == 0 or n_cuts == 0:
            continue
        RNA_types = list(t.RNA_types)
        n_tRNA = RNA_types.count("tRNA")

        if "rRNA" in set(RNA_types):
            t.modifications["rRNA_containing_excision"] = n_cuts
        elif n_tRNA == 1:
            t.modifications["monocistronic_excision"] = n_cuts
        elif n_tRNA > 1:
            t.modifications["polycistronic_wout_rRNA_excision"] = n_cuts
        else:  # only applies to rnpB
            t.modifications["monocistronic_excision"] = n_cuts

        # The non functional RNA segments need degraded back to nucleotides
        # TODO check if RNA_degradation requirement is per nucleotide
        t.modifications["RNA_degradation_machine"] = n_cuts
        t.modifications["RNA_degradation_atp_requirement"] = n_excised
