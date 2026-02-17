from model_setup import *

from math import tanh, sqrt, exp

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
    
    if model.genomic_setup.gene_directions[gene_index] == 1:
        return (v0 / 2.0)*(1.0 - tanh((tau_f - tau_b) / tau_c))
    return (-v0 / 2.0)*(1.0 - tanh((tau_f - tau_b) / tau_c))

def get_RNAP_angular_velocity(model: Model, gene_index: int, x: float, dx_dt: float, left_torque: float, right_torque: float) -> float: # Get the angular velocity of an RNAP based on its position, linear velocity, and the torques on the DNA segments ahead and behind it
    w0 = model.model_setup.w0
    chi = model.model_setup.chi
    eta = model.model_setup.eta
    alpha = model.model_setup.alpha

    x = abs(x - model.genomic_setup.TSSes[gene_index])
    tau_f = right_torque
    tau_b = left_torque
    
    denom = (chi + eta*(x**alpha))
    return ((w0*dx_dt)*(chi / denom)) + ((tau_f - tau_b) / denom)

def get_segment_Lk_dynamics(model: Model, dtheta_dt_front: float, dtheta_dt_back: float) -> float: # Get the rate of change of linking number for a DNA segment based on the angular velocities at its front and back
    return (1.0 / (2.0*3.14))*(dtheta_dt_front - dtheta_dt_back)

def get_RNAP_recruitment_rate(model: Model, TSS_index: int, promoter_status: int, TSS_sigma: float) -> float: # Get the RNAP recruitment rate at a given TSS based on promoter status and local supercoiling
    if promoter_status == 0: # Promoter is OFF
        return 0.0
    
    return model.genomic_setup.RNAP_on_rates[TSS_index]

def get_per_TOP1_binding_rate_for_each_segment(model: Model, segments_lengths: float, segments_sigmas: float) -> list[float]: # Get the TOP1 binding rate for a DNA segment based on its length and supercoiling density
    TOP1_on_rate = model.model_setup.topoisomerase_on_off_rates[0][0] # Note that this is the binding rate for the overall genomic segment, not the per unit length rate

    TOP1_on_rates_per_segment = []
    for i in range(len(segments_lengths)):
        segment_length = segments_lengths[i]
        segment_sigma = segments_sigmas[i]

        TOP1_on_rates_per_segment.append(TOP1_on_rate*(segment_length / (model.genomic_setup.clamp_right - model.genomic_setup.clamp_left)))

    return TOP1_on_rates_per_segment

def get_per_TOP2_binding_rate_for_each_segment(model: Model, segments_lengths: float, segments_sigmas: float) -> list[float]: # Get the TOP2 binding rate for a DNA segment based on its length and supercoiling density
    TOP2_on_rate = model.model_setup.topoisomerase_on_off_rates[1][0] # Note that this is the binding rate for the overall genomic segment, not the per unit length rate

    TOP2_on_rates_per_segment = []
    for i in range(len(segments_lengths)):
        segment_length = segments_lengths[i]
        segment_sigma = segments_sigmas[i]

        TOP2_on_rates_per_segment.append(TOP2_on_rate*(segment_length / (model.genomic_setup.clamp_right - model.genomic_setup.clamp_left)))

    return TOP2_on_rates_per_segment

def get_TOP1_effect_on_Lk_dynamics(model: Model, segment_length: float, segment_sigma: float, segment_torque: float, segment_writhe_frac: float, bound_TOP1_count: int) -> float: # Get the effect of bound TOP1 enzymes on the rate of change of linking number for a DNA segment
    if segment_writhe_frac > 0.0: # Non-zero writhe; TOP1 cannot act
        return 0.0
    k0 = model.model_setup.TOP1_k0
    beta = 1.0 / model.model_setup.kBT
    theta = model.model_setup.TOP1_theta
    tau = abs(segment_torque)

    if segment_sigma > 0.0: # Positively supercoiled DNA
        try:
            return -bound_TOP1_count*k0*exp(theta*beta*tau)*(1.0 - exp(-beta*2*3.14*tau))
        except OverflowError as e:
            print('Overflow error in calculating TOP1 effect on Lk dynamics:' , theta, beta, tau, theta*beta*tau, exp(theta*beta*tau))
            raise e
    else: # Negatively supercoiled DNA
        try:
            return bound_TOP1_count*k0*exp(-theta*beta*tau)*(1.0 - exp(-beta*2*3.14*tau))
        except OverflowError as e:
            print('Overflow error in calculating TOP1 effect on Lk dynamics:' , theta, beta, tau, theta*beta*tau, exp(theta*beta*tau))
            raise e

def get_TOP2_effect_on_Lk_dynamics(model: Model, segment_length: float, segment_sigma: float, segment_torque: float, segment_writhe_frac: float, bound_TOP2_count: int) -> float: # Get the effect of bound TOP2 enzymes on the rate of change of linking number for a DNA segment
    if segment_writhe_frac > 0.0: # Non-zero writhe; TOP2 can act
        Wr = abs(segment_sigma*segment_writhe_frac)
        V0 = model.model_setup.TOP2_V0
        k12 = model.model_setup.TOP2_k12

        return -bound_TOP2_count*V0*(Wr / (Wr + k12)) if segment_sigma > 0.0 else bound_TOP2_count*V0*(Wr / (Wr + k12))
    
    return 0.0

def get_prokaryotic_torque(w0: float, force: float, kBT: float, segment_length: float, sigma: float, finite_size_effect_flag: int, finite_size_effect_length: float) -> tuple[float, int, float, float]: # Get the torque, DNA state, writhe fraction, and sigma_s (threshold beyond which plectonemes form) for prokaryotic DNA based on supercoiling density and force
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

    return (torque, dna_state, writhe_frac, sigma_s)