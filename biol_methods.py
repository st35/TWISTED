from model_setup import *

from math import tanh, sqrt

def get_RNAP_velocity(model: Model, gene_index: int, left_segment_length: float, right_segment_length: float, left_torque: float, right_torque: float) -> float: # Get the velocity of an RNAP based on the torques on the DNA segments ahead and behind it
	v0 = model.model_setup.v0
	tau_c = model.model_setup.tau_c

	tau_f = 0.0
	tau_b = 0.0
	segment_ahead_length = 0.0

	if model.genomic_setup.gene_directions[gene_index] == 1:
		segment_ahead_length = right_segment_length
		tau_f = right_torque
		tau_b = left_torque
	else:
		segment_ahead_length = left_segment_length
		tau_f = left_torque
		tau_b = right_torque
	
	if segment_ahead_length < model.model_setup.between_RNAPs_steric_effect_cutoff:
		return 0.0
	
	return (v0 / 2.0)*(1.0 - tanh((tau_f - tau_b) / tau_c))

def get_RNAP_angular_velocity(model: Model, gene_index: int, x: float, dx_dt: float, left_torque: float, right_torque: float) -> float: # Get the angular velocity of an RNAP based on its position, linear velocity, and the torques on the DNA segments ahead and behind it
    w0 = model.model_setup.w0
    chi = model.model_setup.chi
    eta = model.model_setup.eta
    alpha = model.model_setup.alpha

    x = abs(x - model.genomic_setup.TSSes[gene_index])
    tau_f = 0.0
    tau_b = 0.0
    if model.genomic_setup.gene_directions[gene_index] == 1:
        tau_f = right_torque
        tau_b = left_torque
    else:
        tau_f = left_torque
        tau_b = right_torque
    
    denom = (chi + eta*(x**alpha))
    return ((w0*dx_dt)*(chi / denom)) + ((tau_f - tau_b) / denom)

def get_segment_Lk_dynamics(model: Model, dtheta_dt_front: float, dtheta_dt_back: float) -> float: # Get the rate of change of linking number for a DNA segment based on the angular velocities at its front and back
    return (1.0 / (2.0*3.14))*(dtheta_dt_front - dtheta_dt_back)

def get_RNAP_recruitment_rate(model: Model, TSS_index: int, promoter_status: int, TSS_sigma: float) -> float: # Get the RNAP recruitment rate at a given TSS based on promoter status and local supercoiling
    if promoter_status == 0: # Promoter is OFF
        return 0.0
    
    return model.genomic_setup.RNAP_on_rates[TSS_index]

def get_TOP1_events_rates(model: Model, segments_lengths: list[float], segments_sigmas: list[float]) -> list[float]: # Get the rates of TOP1 binding and unbinding
    return [0.0, 0.0]

def get_TOP2_events_rates(model: Model, segments_lengths: list[float], segments_sigmas: list[float]) -> list[float]: # Get the rates of TOP2 binding and unbinding
    return [0.0, 0.0]

def get_prokaryotic_torque(w0: float, force: float, kBT: float, segment_length: float, sigma: float, finite_size_effect_flag: int, finite_size_effect_length: float) -> tuple[float, int, float]: # Get the torque, DNA state, and writhe fraction for prokaryotic DNA based on supercoiling density and force
    A = 50.0
    C = 95.0
    P = 24.0

    A_m = 4.0
    C_m = 1.75
    e_m = 6.0*kBT
    sigma_0 = -1.0

    c = kBT*C*w0*w0
    p = kBT*P*w0*w0

    cs = c*(1.0 - ((C / (4.0*A))*sqrt(kBT / (A*force))))
    g = force - sqrt((kBT*force / A))

    factor = sqrt((2.0*p*g) / (1.0 - (p / cs)))

    sigma_s = factor / cs
    if finite_size_effect_flag == 1:
        sigma_s = sigma_s*(1.0 + pow((finite_size_effect_length / segment_length), 2.0))
    sigma_p = factor / p

    g_m = 1.2*(force - sqrt(kBT*force / A_m))
    c_m = kBT*C_m*w0*w0

    sigma_sm = (c_m / (cs - c_m))*(-sigma_0 - sqrt(sigma_0*sigma_0 + (2.0*(cs - c_m) / (cs*c_m))*(g + e_m - g_m)))
    sigma_m = sigma_0 + (cs / (cs - c_m))*(-sigma_0 - sqrt(sigma_0*sigma_0 + (2.0*(cs - c_m) / (cs*c_m))*(g + e_m - g_m)))

    torque = 0.0
    dna_state = -1
    writhe_frac = 0.0

    if sigma <= sigma_m:
        torque = (c_m / w0)*(sigma - sigma_0)
        dna_state = 0 # Twisted, melted DNA
    elif sigma > sigma_m and sigma <= sigma_sm:
        torque = (c_m / w0)*(sigma_m - sigma_0)
        dna_state = 1 # Melted DNA
    elif sigma > sigma_sm and sigma <= sigma_s:
        torque = (cs / w0)*sigma
        dna_state = 2 # Twisted DNA
    elif sigma > sigma_s and sigma <= sigma_p:
        torque = factor / w0
        dna_state = 5 # Positive plectoneme
        writhe_frac = (sigma - sigma_s) / (sigma_p - sigma_s)
    elif sigma > sigma_p:
        torque = (p / w0)*sigma
        if torque > 40.0:
            torque = 40.0
        dna_state = 6 # Twisted plectoneme
        writhe_frac = 1.0

    return (torque, dna_state, writhe_frac)