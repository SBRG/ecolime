from cobrame import (ComplexData, Complex, GenericData, ModificationData,
                     StoichiometricData)
from cobrame.util import building

# these guys can transfer assembled iron sulfur clusters to the various enzymes
fes_transfer = {"erpA": "CPLX0-7617", "iscA": "IscA_tetra",
                "sufA": "CPLX0-7824"}
fes_synthesizing_complexes = {
    "sufBC2DES_pathway_complex": {"CPLX0-1341": 1,
                                  "CPLX0-246_CPLX0-1342_mod_pydx5p": 1},

    "iscUS_cyaY_pathway_complex": {"IscU": 1,
                                   "IscS_mod_2:pydx5p": 1,
                                   "EG11653-MONOMER": 1}}

fes_transfer_reactions = {
    "suf_2fe2s_formation": {"enzyme": "sufBC2DES_pathway_complex",
                            "stoich": {'2fe2s_c': 1.0,
                                       'ala__L_c': 2.0,
                                       'cys__L_c': -2.0,
                                       'fad_c': 1.0,
                                       'fadh2_c': -1.0,
                                       'fe3_c': -2.0,
                                       'h_c': 6.0}},
    "suf_4fe4s_formation": {"enzyme": "sufBC2DES_pathway_complex",
                            "stoich": {'4fe4s_c': 1.0,
                                       'ala__L_c': 4.0,
                                       'cys__L_c': -4.0,
                                       'fad_c': 3.0,
                                       'fadh2_c': -3.0,
                                       'fe3_c': -4.0,
                                       'h_c': 6.0}},
    'isc_2fe2s_formation': {"enzyme": "iscUS_cyaY_pathway_complex",
                            "stoich": {'2fe2s_c': 1.0, 'ala__L_c': 2.0,
                                       'cys__L_c': -2.0, 'fe2_c': -2.0,
                                       'h_c': 4.0}},
    'isc_4fe4s_formation': {"enzyme": "iscUS_cyaY_pathway_complex",
                            "stoich": {'4fe4s_c': 1.0,
                                       'ala__L_c': 4.0,
                                       'cys__L_c': -4.0,
                                       'fad_c': 1.0,
                                       'fadh2_c': -1.0,
                                       'fe2_c': -4.0,
                                       'h_c': 6.0}}}

# Add known specific chaperone
fes_chaperones = {'CPLX0-1762': 'G6712-MONOMER'}  # FE-S modification

# complexes that can transfer an iron sulfur cluster to an enzyme target
generic_fes_transfer_complexes = ['CPLX0-7617', 'CPLX0-7824', 'IscA_tetra']

lipoate_modifications = {
    "mod_lipo_c": {"enzyme": 'EG11796-MONOMER',
                   "stoich": {"lipoamp_c": -1,
                              "amp_c": 1}},

    "mod_lipo_c_alt": {"enzyme": 'EG11591-MONOMER',
                       "stoich": {'EG50003-MONOMER_mod_pan4p_mod_lipo': -1,
                                  'EG50003-MONOMER_mod_pan4p': 1,
                                  'h_c': -1}}

                         }

bmocogdp_chaperones = {'TMAOREDUCTI-CPLX': 'EG12195-MONOMER',
                       'DIMESULFREDUCT-CPLX': 'G6849-MONOMER',
                       'NITRATREDUCTA-CPLX': 'NARJ-MONOMER',
                       'NITRATREDUCTZ-CPLX': 'NARW-MONOMER',
                       'NAP-CPLX': 'NAPD-MONOMER',
                       'NAPAB-CPLX_NAPC-MONOMER': 'NAPD-MONOMER'}


def add_iron_sulfur_modifications(me_model):

    for i in generic_fes_transfer_complexes:
        me_model.add_metabolites([Complex(i)])

    generic_fes_transfer = GenericData("generic_fes_transfer", me_model,
                                       generic_fes_transfer_complexes)
    generic_fes_transfer.create_reactions()

    # add fes transfer enzymes to proper modification data
    me_model.modification_data.mod_2fe2s_c.enzyme = generic_fes_transfer.id
    me_model.modification_data.mod_2fe2s_c.keff = 65.
    me_model.modification_data.mod_4fe4s_c.enzyme = generic_fes_transfer.id
    me_model.modification_data.mod_4fe4s_c.keff = 65.

    for chaperone in set(fes_chaperones.values()):
        new_mod = ModificationData('mod_2fe2s_c_' + chaperone, me_model)
        new_mod.enzyme = chaperone
        new_mod.stoichiometry = {'2fe2s_c': -1}

    for cplx_data in me_model.modification_data.get_by_id(
            'mod_2fe2s_c').get_complex_data():
        cplx_id = cplx_data.id.split('_mod')[0]
        if cplx_id in fes_chaperones:
            cplx_data.modifications['mod_2fe2s_c_' + fes_chaperones[
                cplx_id]] = \
                cplx_data.modifications.pop('mod_2fe2s_c')


def add_iron_sulfur_reactions_and_complexes(me_model):
    for fes_complex, stoich in fes_synthesizing_complexes.items():
        complex_data = ComplexData(fes_complex, me_model)
        complex_data.stoichiometry = stoich
        complex_data.create_complex_formation()

    for key, values in fes_transfer_reactions.items():
        # Define stoichiometric data
        stoich_data = StoichiometricData(key, me_model)
        stoich_data._stoichiometry = values['stoich']
        stoich_data.lower_bound = 0.
        stoich_data.upper_bound = 1000.

        # Create MetabolicReaction and associate stoichiometric data
        complex_id = values['enzyme']
        building.add_metabolic_reaction_to_model(me_model, stoich_data.id,
                                                 'forward',
                                                 complex_id=complex_id,
                                                 update=True)


def add_lipoate_modifications(me_model):
    # two different reactions can add a lipoate modification.
    # We create a separate ModificationData for each one

    for mod, info in lipoate_modifications.items():
        if mod in me_model.modification_data:
            mod_data = me_model.modification_data.get_by_id(mod)
        else:
            mod_data = ModificationData(mod, me_model)

        mod_data.stoichiometry = info["stoich"]
        mod_data.enzyme = info["enzyme"]

    lipo = me_model.modification_data.get_by_id('mod_lipo_c')
    alt_lipo = me_model.modification_data.get_by_id('mod_lipo_c_alt')
    for cplx_data in lipo.get_complex_data():
        alt_cplx_data = ComplexData(cplx_data.id + "alt", me_model)
        alt_cplx_data.complex_id = cplx_data.complex_id
        alt_cplx_data.stoichiometry = cplx_data.stoichiometry
        alt_cplx_data.chaperones = cplx_data.chaperones
        alt_cplx_data.modifications = cplx_data.modifications.copy()
        alt_cplx_data.modifications[alt_lipo.id] = \
            alt_cplx_data.modifications.pop(lipo.id)
        alt_cplx_data.create_complex_formation()


def add_bmocogdp_modifications(me_model):
    for chaperone in set(bmocogdp_chaperones.values()):
        new_mod = ModificationData('mod_bmocogdp_c_' + chaperone, me_model)
        new_mod.enzyme = chaperone
        new_mod.stoichiometry = {'bmocogdp_c': -1}

    for cplx_data in me_model.modification_data.get_by_id(
            'mod_bmocogdp_c').get_complex_data():
        cplx_id = cplx_data.id.split('_mod')[0]
        if cplx_id in bmocogdp_chaperones:
            cplx_data.modifications[
                'mod_bmocogdp_c_' + bmocogdp_chaperones[cplx_id]] = \
                cplx_data.modifications.pop('mod_bmocogdp_c')
