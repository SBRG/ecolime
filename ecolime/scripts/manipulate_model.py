def change_substrate(model, subrate_rxn, source='C', value=1000):
    """value should be positive"""
    if source == 'C':
        model.reactions.EX_glc__D_e.lower_bound = 0
    model.reactions.get_by_id(subrate_rxn).lower_bound = -value


def change_PO_ratio(model, value):
    rxn = model.reactions.ndh_flux_split_constraint
    rxn.add_metabolites({'ndh1_constraint': -value,
                         'ndh2_constraint': -(1-value)}, combine=False)

def modify_model_parameter(model, parameter, mult_value, default):
    if parameter == 'Q':
        Q = default
        model.unmodeled_protein_fraction = Q * mult_value
    elif parameter == 'ngam':
        LB = default
        model.reactions.ATPM_FWD_SPONT.lower_bound = LB * mult_value
    elif parameter == 'Inner_Membrane':
        membrane = default
        value = default * mult_value
        rxn_prefix = 'SA_components_to_SA_'
        rxn = model.reactions.get_by_id(rxn_prefix + membrane)
        rxn.add_metabolites({'SA_total_protein_' + membrane: -value,
                             'SA_nonprotein_' + membrane: -(1-value)},
                            combine=False)
    elif parameter == 'gam':
        gam_dict = {'atp_c': -(default * mult_value),
                    'h2o_c': -(default * mult_value),
                    'adp_c': (default * mult_value),
                    'h_c': (default * mult_value),
                    'pi_c': (default * mult_value)}
        model.reactions.biomass_component_demand.add_metabolites(gam_dict,
                                                      combine=False)
    else:
        raise('%s not a valid parameter type' % parameter)