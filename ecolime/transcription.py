from minime import ComplexData, TranscribedGene, ModificationData
from minime.util.building import add_modification_data

import cobra


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
    'sigma_19': 'CPLX0-221',
    # this is standard notation
    'sigma_28': 'CPLX0-222',
    # this is standard notation
    'sigma_32': 'RNAP32-CPLX',
    # this is standard notation
    'sigma_54': 'RNAP54-CPLX',
    # this is standard notation
    'sigma_70': 'RNAP70-CPLX',
    # this is standard notation
    'sigma_24': 'RNAPE-CPLX',
    # this is standard notation
    'sigma_38': 'RNAPS-CPLX',
    # this is standard notation
    'RPOD-MONOMER': 'RNAP70-CPLX',
    # this is the notation used in promoters.dat
    'RPOH-MONOMER': 'RNAP32-CPLX',
    # this is the notation used in promoters.dat
    'RPOE-MONOMER': 'RNAPE-CPLX',
    # this is the notation used in promoters.dat
    'RPOS-MONOMER': 'RNAPS-CPLX',
    # this is the notation used in promoters.dat
    'RPON-MONOMER': 'RNAP54-CPLX',
    # this is the notation used in promoters.dat
    'EG11355-MONOMER': 'CPLX0-222',
    # this is the notation used in promoters.dat
    'PD00440': 'CPLX0-221',
    # this is the notation used in promoters.dat
    'RpoD_mono': 'RNAP70-CPLX',
    # this is notation used in Tu_from_ecocyc, protein_complexes (and below)
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


def add_RNA_polymerase_complexes(me_model, verbose=True):

    for complex, components in rna_polymerase_sigma_factor_components.iteritems():
        rnap_complex = ComplexData(complex, me_model)
        rnap_components = rnap_complex.stoichiometry
        sigma_factor = components['sigma_factor']
        polymerase = components['polymerase']

        rnap_components[sigma_factor] = 1
        rnap_components[polymerase] = 1

        rnap_complex.create_complex_formation()
