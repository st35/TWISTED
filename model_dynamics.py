from model_setup import *
from biol_methods import *
from utilities import *

from scipy.integrate import solve_ivp
from math import log
import numpy as np

def get_state_vectors_from_dicts(model: Model) -> tuple[list[int], list[float]]: # Get the RNAP_gene_index and state_vector from the model's x_dict, theta_dict, and Lk RNAPs are ordered from right to left on the DNA
	RNAP_gene_index = []
	x_vector = []
	theta_vector = []
	for i in range(len(model.genomic_setup.gene_names)):
		RNAP_gene_index = RNAP_gene_index + [i for _ in model.x_dict[i]]
		if model.genomic_setup.gene_directions[i] == 1:
			x_vector = x_vector + model.x_dict[i]
			theta_vector = theta_vector + model.theta_dict[i]
		else:
			x_vector = x_vector + model.x_dict[i][::-1]
			theta_vector = theta_vector + model.theta_dict[i][::-1]
	
	rates_vector, _ = get_events_rates(model, RNAP_gene_index, x_vector + theta_vector + model.Lk + [0.0])

	return (RNAP_gene_index, x_vector + theta_vector + model.Lk + [sum(rates_vector)])

def update_dicts_from_state_vector(model: Model, RNAP_gene_index: list[int], state_vector: list[float]) -> None: # Update the model's x_dict, theta_dict, and Lk from the given RNAP_gene_index and state_vector
	RNAP_count = len(RNAP_gene_index)
	x_vector = state_vector[0:RNAP_count]
	theta_vector = state_vector[RNAP_count:2*RNAP_count]
	Lk_vector = state_vector[2*RNAP_count:-1]

	for i in range(len(model.genomic_setup.gene_names)):
		model.x_dict[i] = []
		model.theta_dict[i] = []

	for i in range(RNAP_count):
		gene_index = RNAP_gene_index[i]
		model.x_dict[gene_index].append(x_vector[i])
		model.theta_dict[gene_index].append(theta_vector[i])
	for i in range(len(model.genomic_setup.gene_names)):
		if model.genomic_setup.gene_directions[i] == -1: # Reverse order for - strand genes
			model.x_dict[i] = model.x_dict[i][::-1]
			model.theta_dict[i] = model.theta_dict[i][::-1]
	
	model.Lk = [val for val in Lk_vector]

def calculate_segments_attributes(model: Model, RNAP_gene_index: list[int], state_vector: list[float]) -> tuple[list[float], list[float], list[float], list[int], list[float]]: # Calculate and return the segments_lengths, segments_sigmas, segments_torques, segments_dna_states, and segments_writhe_fractions based on the current state_vector; segments are ordered from right to left on the DNA
	RNAP_count = len(RNAP_gene_index)
	x_vector = state_vector[0:RNAP_count]
	theta_vector = state_vector[RNAP_count:2*RNAP_count]
	Lk_vector = state_vector[2*RNAP_count:-1]

	segments_lengths = []
	if RNAP_count == 0: # No RNAPs on the DNA; single segment between clamps
		segments_lengths.append(model.genomic_setup.clamp_right - model.genomic_setup.clamp_left)
	else:
		for i in range(RNAP_count):
			if i == 0:
				segments_lengths.append(model.genomic_setup.clamp_right - x_vector[i])
			else:
				segments_lengths.append(x_vector[i - 1] - x_vector[i])
		segments_lengths.append(x_vector[-1] - model.genomic_setup.clamp_left)
	
	segments_LK0 = []
	for length in segments_lengths:
		segments_LK0.append(length / model.model_setup.h_dna)
	
	segments_sigmas = []
	for Lk, Lk0 in zip(Lk_vector, segments_LK0):
		segments_sigmas.append((Lk - Lk0) / Lk0)
	
	segments_torques = []
	for segment_length, segment_sigma in zip(segments_lengths, segments_sigmas):
		if model.genomic_setup.chromatin_type == 'prokaryotic':
			segments_torques.append(get_prokaryotic_torque(model.model_setup.w0, model.model_setup.force, model.model_setup.kBT, segment_length, segment_sigma, model.model_setup.finite_size_effect_flag, model.model_setup.finite_size_effect_length))
		else:
			raise NotImplementedError('Chromatin type "eukaryotic" not yet implemented.')
	
	segments_writhe_fractions = [val[2] for val in segments_torques]
	segments_dna_states = [val[1] for val in segments_torques]
	segments_torques = [val[0] for val in segments_torques]
	
	return (segments_lengths, segments_sigmas, segments_torques, segments_dna_states, segments_writhe_fractions)

def get_RNAP_velocities(model: Model, RNAP_gene_index: list[int], segments_lengths: list[float], segments_torques: list[float]) -> list[float]: # Get the velocities of all RNAPs based on the segments lengths and torques
	RNAP_count = len(RNAP_gene_index)
	if RNAP_count == 0: # No RNAPs on the DNA; return empty list
		return []
	
	RNAP_velocities = []
	left_segment_length = 0.0
	right_segment_length = 0.0
	left_torque = 0.0
	right_torque = 0.0
	for i in range(RNAP_count):
		left_segment_length = segments_lengths[i + 1]
		right_segment_length = segments_lengths[i]
		left_torque = segments_torques[i + 1]
		right_torque = segments_torques[i]

		RNAP_velocities.append(get_RNAP_velocity(model, RNAP_gene_index[i], left_segment_length, right_segment_length, left_torque, right_torque))
	
	return RNAP_velocities

def get_RNAP_angular_velocities(model: Model, RNAP_gene_index: list[int], state_vector: list[float], segments_lengths: list[float], segments_torques: list[float], RNAP_velocities: list[float]) -> list[float]: # Get the angular velocities of all RNAPs based on the state_vector, segments lengths and torques, and RNAP linear velocities
	RNAP_count = len(RNAP_gene_index)
	if RNAP_count == 0: # No RNAPs on the DNA; return empty list
		return []
	
	x_vector = state_vector[0:RNAP_count]
	
	RNAP_angular_velocities = []
	left_torque = 0.0
	right_torque = 0.0
	for i in range(RNAP_count):
		left_torque = segments_torques[i + 1]
		right_torque = segments_torques[i]

		RNAP_angular_velocities.append(get_RNAP_angular_velocity(model, RNAP_gene_index[i], x_vector[i], RNAP_velocities[i], left_torque, right_torque))
	
	return RNAP_angular_velocities

def get_segments_Lk_dynamics(model: Model, dtheta_dt: list[float]) -> list[float]: # Get the rates of change of linking number for all DNA segments based on the RNAP angular velocities
	RNAP_count = len(dtheta_dt)
	if RNAP_count == 0: # No RNAPs on the DNA; there is only one segment between clamps and there is no change in linking number
		return [0.0]
	
	dLk_dt = []
	dtheta_dt_front = 0.0
	dtheta_dt_back = 0.0
	for i in range(RNAP_count + 1):
		if i == 0: # Rightmost segment
			dtheta_dt_front = 0.0
			dtheta_dt_back = dtheta_dt[i]
		elif i == RNAP_count: # Leftmost segment
			dtheta_dt_front = dtheta_dt[i - 1]
			dtheta_dt_back = 0.0
		else:
			dtheta_dt_front = dtheta_dt[i - 1]
			dtheta_dt_back = dtheta_dt[i]
		dLk_dt.append(get_segment_Lk_dynamics(model, dtheta_dt_front, dtheta_dt_back))
	
	return dLk_dt

def get_RNAP_recruitment_rates(model: Model, RNAP_gene_index: list[int], state_vector: list[float], segments_lengths: list[float], segments_sigmas: list[float]) -> list[float]: # Get the RNAP recruitment rates at all TSSes based on the current state_vector and segments attributes
	RNAP_count = len(RNAP_gene_index)
	x_vector = state_vector[0:RNAP_count]

	TSS_segments_indices = [get_spot_segment_index(model.genomic_setup.TSSes[i], segments_lengths) for i in range(len(model.genomic_setup.gene_names))] # Get the segment indices for all TSSes
	is_TSS_blocked = [get_TSS_steric_hindrance_status(model.genomic_setup.TSSes[i], RNAP_gene_index, state_vector, model.model_setup.between_RNAPs_steric_effect_cutoff) for i in range(len(model.genomic_setup.gene_names))] # Get the steric hindrance status for all TSSes

	RNAP_recruitment_rates = []
	for i in range(len(model.genomic_setup.gene_names)):
		if False: # Implement additional conditions for blocking RNAP recruitment if needed
			RNAP_recruitment_rates.append(0.0)
		else:
			RNAP_recruitment_rates.append(get_RNAP_recruitment_rate(model, i, model.promoter_status[i], segments_sigmas[TSS_segments_indices[i]]))

	return RNAP_recruitment_rates

def are_RNAPs_alive(model: Model, RNAP_gene_index: list[int], state_vector: list[float]) -> list[int]: # Determine whether each RNAP is still alive (1) or has finished transcription (0) based on the current state_vector
	RNAP_count = len(RNAP_gene_index)
	x_vector = state_vector[0:RNAP_count]

	RNAPs_alive_status = []
	for i in range(RNAP_count):
		gene_direction = model.genomic_setup.gene_directions[RNAP_gene_index[i]]
		TSS_position = model.genomic_setup.TSSes[RNAP_gene_index[i]]
		gene_length = model.genomic_setup.gene_lengths[RNAP_gene_index[i]]
		if gene_direction == 1:
			if x_vector[i] > TSS_position + gene_length:
				RNAPs_alive_status.append(0)
			else:
				RNAPs_alive_status.append(1)
		else:
			if x_vector[i] < TSS_position - gene_length:
				RNAPs_alive_status.append(0)
			else:
				RNAPs_alive_status.append(1)
	
	return RNAPs_alive_status

def update_Lk_vector_after_RNAP_recruitment(model: Model, TSS_index: int, RNAP_gene_index: list[int], state_vector: list[float], segments_lengths: list[float], segments_sigmas: list[float]) -> None: # Update the model's Lk vector after an RNAP is recruited at the given TSS_index; Lk for the segments on either side of the TSS are calculated such that the supercoiling density in the segments is the same as in the segment spanning the TSS before recruitment
	TSS_segment_index = get_spot_segment_index(model.genomic_setup.TSSes[TSS_index], segments_lengths)
	TSS_sigma = segments_sigmas[TSS_segment_index]

	left_segment_boundary = model.genomic_setup.clamp_left
	right_segment_boundary = model.genomic_setup.clamp_right

	if len(segments_lengths) > 1:
		if TSS_segment_index == 0:
			left_segment_boundary = state_vector[0]
		elif TSS_segment_index == len(segments_lengths) - 1:
			right_segment_boundary = state_vector[len(RNAP_gene_index) - 1]
		else:
			left_segment_boundary = state_vector[TSS_segment_index]
			right_segment_boundary = state_vector[TSS_segment_index - 1]
	
	left_segment_length = model.genomic_setup.TSSes[TSS_index] - left_segment_boundary
	right_segment_length = right_segment_boundary - model.genomic_setup.TSSes[TSS_index]
	assert left_segment_length >= 0.0 and right_segment_length >= 0.0, 'Error in calculating segment lengths to the left and right of the TSS during RNAP recruitment.'

	left_Lk0 = left_segment_length / model.model_setup.h_dna
	right_Lk0 = right_segment_length / model.model_setup.h_dna

	# Based on sigma = (Lk - Lk0) / Lk0 => Lk = Lk0 * (1 + sigma)
	left_Lk = left_Lk0*(1.0 + TSS_sigma)
	right_Lk = right_Lk0*(1.0 + TSS_sigma)

	model.Lk = model.Lk[:TSS_segment_index] + [right_Lk, left_Lk] + model.Lk[TSS_segment_index + 1:]

def update_state_vector_to_remove_dead_RNAPs(model: Model, RNAP_gene_index: list[int], t: float, state_vector: list[float], simulation_setup_and_state: SimulationSetupAndState) -> None: # Update the state_vector and RNAP_gene_index to remove RNAPs that have finished transcription; also update the simulation_setup_and_state with transcription completion data
	RNAP_count = len(RNAP_gene_index)
	if RNAP_count == 0: # No RNAPs on the DNA; nothing to update
		return

	x_vector = state_vector[0:RNAP_count]
	theta_vector = state_vector[RNAP_count:2*RNAP_count]
	Lk_vector = state_vector[2*RNAP_count:-1]
	
	RNAPs_alive_status = are_RNAPs_alive(model, RNAP_gene_index, state_vector)
	for i in range(RNAP_count):
		if RNAPs_alive_status[i] == 0: # RNAP has finished transcription; update simulation_setup_and_state accordingly
			simulation_setup_and_state.RNAPs_finished_transcription[RNAP_gene_index[i]] += 1
			simulation_setup_and_state.RNAPs_exit_positions[RNAP_gene_index[i]].append(x_vector[i])
			simulation_setup_and_state.RNAP_exit_times[RNAP_gene_index[i]].append(t)

	new_RNAP_gene_index = [RNAP_gene_index[i] for i in range(RNAP_count) if RNAPs_alive_status[i] == 1]
	new_x_vector = [x_vector[i] for i in range(RNAP_count) if RNAPs_alive_status[i] == 1]
	new_theta_vector = [theta_vector[i] for i in range(RNAP_count) if RNAPs_alive_status[i] == 1]

	# When an RNAP finishes transcription, the two segments on either side of it merge into one segment; thus, we need to update the Lk_vector accordingly. Lk of the merged segment is the sum of the Lk of the two segments.
	new_Lk_vector = []
	RNAP_index = 0
	while RNAP_index < RNAP_count:
		Lk_front_segment = Lk_vector[RNAP_index]
		Lk_back_segment = Lk_vector[RNAP_index + 1]
		if RNAPs_alive_status[RNAP_index] == 1:
			new_Lk_vector.append(Lk_front_segment)
			RNAP_index += 1
		else:
			new_Lk_vector.append(Lk_front_segment + Lk_back_segment)
			RNAP_index += 2
	if RNAPs_alive_status[-1] == 1:
		new_Lk_vector.append(Lk_vector[-1])
	
	state_vector[:] = new_x_vector + new_theta_vector + new_Lk_vector + [state_vector[-1]]
	RNAP_gene_index[:] = new_RNAP_gene_index

def get_events_rates(model: Model, RNAP_gene_index: list[int], state_vector: list[float]) -> tuple[list[float], list[int]]: # Get the rates of all possible events and the indices that separate different event types in the rates_vector
	segments_lengths, segments_sigmas, segments_torques, segments_dna_states, segments_writhe_fractions = calculate_segments_attributes(model, RNAP_gene_index, state_vector)

	RNAP_recruitment_rates = get_RNAP_recruitment_rates(model, RNAP_gene_index, state_vector, segments_lengths, segments_torques)
	model_observation_event_rate = [model.model_setup.model_observation_event_rate]
	global_supercoiling_relaxation_rate = [model.model_setup.global_supercoiling_relaxation_rate]
	local_supercoiling_relaxation_rates = model.model_setup.local_supercoiling_relaxation_rates
	TOP1_rates = get_TOP1_events_rates(model, segments_lengths, segments_sigmas)
	TOP2_rates = get_TOP2_events_rates(model, segments_lengths, segments_sigmas)

	rates_vector = RNAP_recruitment_rates + model_observation_event_rate + global_supercoiling_relaxation_rate + local_supercoiling_relaxation_rates + TOP1_rates + TOP2_rates
	assert all(rate >= 0.0 for rate in rates_vector), 'Negative rate encountered in calculating events rates.'

	events_indices = []
	events_indices.append(len(RNAP_recruitment_rates))
	events_indices.append(events_indices[-1] + len(model_observation_event_rate))
	events_indices.append(events_indices[-1] + len(global_supercoiling_relaxation_rate))
	events_indices.append(events_indices[-1] + len(local_supercoiling_relaxation_rates))
	events_indices.append(events_indices[-1] + len(TOP1_rates))
	events_indices.append(events_indices[-1] + len(TOP2_rates))

	return rates_vector, events_indices

def model_dynamics(t: float, state_vector: list[float], RNAP_gene_index: list[int], model: Model, simulation_setup_and_state: SimulationSetupAndState) -> list[float]: # Calculate the time derivatives of the state_vector and the cumulative propensity for event selection
	RNAP_count = len(RNAP_gene_index)

	segments_lengths, segments_sigmas, segments_torques, segments_dna_states, segments_writhe_fractions = calculate_segments_attributes(model, RNAP_gene_index, state_vector)

	rates_vector, _ = get_events_rates(model, RNAP_gene_index, state_vector)

	dx_dt = get_RNAP_velocities(model, RNAP_gene_index, segments_lengths, segments_torques)
	dtheta_dt = get_RNAP_angular_velocities(model, RNAP_gene_index, state_vector, segments_lengths, segments_torques, dx_dt)
	dLk_dt = get_segments_Lk_dynamics(model, dtheta_dt)

	return dx_dt + dtheta_dt + dLk_dt + [sum(rates_vector)]

def integrate(model: Model, simulation_setup_and_state: SimulationSetupAndState, t_start: float, state_vector: list[float], RNAP_gene_index: list[int], p0: float) -> None: # Integrate the model dynamics until the next event occurs; return the time step dt and cumulative propensity a0 at event time
	dt = simulation_setup_and_state.RNAP_alive_status_check_interval
	t = t_start
	a0_ini = state_vector[-1]
	t_event = -1.0
	t_event_index = -1

	while True:
		t_eval = []
		t_first = t
		while t_first < t + dt:
			t_eval.append(t_first)
			t_first += simulation_setup_and_state.integration_time_resolution
		t_eval.append(t + dt)
		sol = solve_ivp(lambda t, y: model_dynamics(t, y, RNAP_gene_index, model, simulation_setup_and_state), (t, t + dt), state_vector, t_eval = t_eval, method = 'RK45')
		assert sol.success, 'ODE integration failed during simulation.'
		assert np.all(np.diff(sol.y[-1, :]) >= 0.0), 'Cumulative propensity decreased during integration.' # Cumulative propensity should never decrease over time
		for t_index, t_val in enumerate(sol.t): # Check for event occurrence at each time point in the solution
			a0 = sol.y[-1, t_index] - a0_ini
			assert a0 >= 0.0, 'Cumulative propensity increment is negative during integration.'
			if a0 > log(1.0 / p0):
				t_event = t_val
				t_event_index = t_index
				dt_current = t_event - t_start
				break
		if t_event > 0.0: # Event has occurred within this integration interval
			state_vector[:] = sol.y[:, t_event_index]
			update_state_vector_to_remove_dead_RNAPs(model, RNAP_gene_index, t + dt_current, state_vector, simulation_setup_and_state)
			return (dt_current, a0)
		# No event has occurred; update state_vector and continue integration
		state_vector[:] = sol.y[:, -1]
		update_state_vector_to_remove_dead_RNAPs(model, RNAP_gene_index, t + dt, state_vector, simulation_setup_and_state)
		t += dt