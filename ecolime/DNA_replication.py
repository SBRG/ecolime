import math
from scipy.optimize import leastsq
from minime import mu
from minime.util.mass import dna_mw_no_ppi

# Experimental Data
gr_data_doublings_per_hour = [0, 0.6, 1.0, 1.5, 2.0, 2.5]
gr_data = [m * math.log(2) for m in gr_data_doublings_per_hour]
#first point is mass 1 genome in 0.4 um^3 at density
percent_dna_data = [0.0592, 0.0512, 0.0330, 0.0252, 0.0222, 0.0208]

DNA_polymerase_stoichiometry = {
    "b3702": 1,  # dnaA (initiator)
    "b4052": 6,  # dnaB (helicase)
    "b4361": 1,  # dnaC
    "b3066": 3,  # dnaG (primase)
    # CPLXO-3803 (dna polymerase iii, holoenzyme)
    # core
    "b0184": 3,  # dnaE
    "b0215": 3,  # dnaQ
    "b1842": 3,  # holE
    # holo
    "b0470": 5,  # dnaX
    "b1099": 1,  # holB
    "b0640": 1,  # holA
    "b3701": 4,  # dnaN
    "b4259": 4,  # holC
    "b4372": 4,  # holD

    # CPLX0-2425 (dna gyrase)
    "b2231": 2,  # gyrA
    "b3699": 2,  # gyrB
}


def percent_dna_template_function(params, gr):
    [g_p_gDW_0, g_per_gDW_inf, b, d] = params
    c = g_per_gDW_inf
    a = g_p_gDW_0 - g_per_gDW_inf
    g_p_gDW = (-a * gr ** d) / (b + gr ** d) + a + c
    return g_p_gDW


def optimize_dna_function(gr, percent_dna):
    params = [0.9, 0.3, 0.2, 1]

    def minimization_function(params, gr, percent_dna):
        return percent_dna_template_function(params, gr) - percent_dna
    a = leastsq(minimization_function, params, args=(gr, percent_dna))
    return a[0]


def return_gr_dependent_dna_demand(GC_fraction):
    fit_params = optimize_dna_function(gr_data, percent_dna_data)
    dna_g_per_g = percent_dna_template_function(fit_params, mu)

    # average dinucleotide molecular weight
    dna_mw = dna_mw_no_ppi
    dNTP_mw = (GC_fraction * (dna_mw['dctp'] + dna_mw['dgtp']) +
               (1 - GC_fraction) * (dna_mw['datp'] + dna_mw['dttp'])) / 2

    mmol_dNTPs_per_gram_dry_weight = dna_g_per_g / dNTP_mw * 1000

    dna_demand_stoich = {
        'datp_c': -mmol_dNTPs_per_gram_dry_weight * ((1 - GC_fraction) / 2),
        'dctp_c': -mmol_dNTPs_per_gram_dry_weight * (GC_fraction / 2),
        'dgtp_c': -mmol_dNTPs_per_gram_dry_weight * (GC_fraction / 2),
        'dttp_c': -mmol_dNTPs_per_gram_dry_weight * ((1 - GC_fraction) / 2),
        'ppi_c': mmol_dNTPs_per_gram_dry_weight}

    return dna_demand_stoich
