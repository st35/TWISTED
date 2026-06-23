from model_setup import *
from biol_methods import *
from utilities import *

import sys
from scipy.integrate import solve_ivp
from math import log
import numpy as np
from typing import Callable

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

def calculate_segments_attributes(model: Model, RNAP_gene_index: list[int], state_vector: list[float]) -> tuple[list[float], list[float], list[float], list[int], list[float], list[float]]: # Calculate and return the segments_lengths, segments_sigmas, segments_torques, segments_dna_states, segments_writhe_fractions, and segments_plectoneme_thresholds based on the current state_vector; segments are ordered from right to left on the DNA
	RNAP_count = len(RNAP_gene_index)
	x_vector = state_vector[0:RNAP_count]
	theta_vector = state_vector[RNAP_count:2*RNAP_count]
	Lk_vector = state_vector[2*RNAP_count:-1]

	topological_barrier_proteins_positions = []
	for i in range(len(model.binding_proteins)):
		if model.binding_proteins[i].is_topological_barrier:
			topological_barrier_proteins_positions = topological_barrier_proteins_positions + model.binding_proteins_positions[i]
	x_vector = list(x_vector) + topological_barrier_proteins_positions
	x_vector = sorted(x_vector, reverse = True)

	segments_lengths = []
	if len(x_vector) == 0: # No RNAPs or topological barrier proteins on the DNA; single segment between clamps
		segments_lengths.append(model.genomic_setup.clamp_right - model.genomic_setup.clamp_left)
	else:
		for i in range(len(x_vector)):
			if i == 0:
				segments_lengths.append(model.genomic_setup.clamp_right - x_vector[i])
			else:
				segments_lengths.append(x_vector[i - 1] - x_vector[i])
		segments_lengths.append(x_vector[-1] - model.genomic_setup.clamp_left)
	for segment_length in segments_lengths:
		if segment_length < 0.0:
			raise ValueError('Negative segment length calculated, which is invalid.')
	assert len(segments_lengths) == len(x_vector) + 1, 'Number of segments should be one more than the number of boundaries (RNAPs and topological barrier proteins) on the DNA.'
	
	segments_LK0 = []
	for length in segments_lengths:
		segments_LK0.append(length / model.model_setup.h_dna)
	
	segments_sigmas = []
	for Lk, Lk0 in zip(Lk_vector, segments_LK0):
		segments_sigmas.append((Lk - Lk0) / Lk0)
	
	segments_nucleosome_densities = [0.0 for _ in segments_lengths]
	if model.genomic_setup.chromatin_type == 'eukaryotic':
		for i in range(len(segments_lengths)):
			segments_nucleosome_densities[i] = get_nucleosome_occupied_fraction_per_segment(model, segments_lengths, i)
			assert segments_nucleosome_densities[i] >= 0.0 and segments_nucleosome_densities[i] <= 1.0, 'Invalid nucleosome density for segment ' + str(i) + ': ' + str(segments_nucleosome_densities[i])
	
	segments_torques = []
	for segment_length, segment_sigma, segment_psi in zip(segments_lengths, segments_sigmas, segments_nucleosome_densities):
		if model.genomic_setup.chromatin_type == 'prokaryotic':
			segments_torques.append(get_prokaryotic_torque(model.model_setup.w0, model.model_setup.force, model.model_setup.kBT, segment_length, segment_sigma, model.model_setup.finite_size_effect_flag, model.model_setup.finite_size_effect_length))
		else:
			segments_torques.append(get_eukaryotic_torque(model.model_setup.force, segment_length, segment_psi, segment_sigma, model.model_setup.finite_size_effect_flag, model.model_setup.finite_size_effect_length))
	
	segments_plectoneme_thresholds = [val[3] for val in segments_torques]
	segments_writhe_fractions = [val[2] for val in segments_torques]
	segments_dna_states = [val[1] for val in segments_torques]
	segments_torques = [val[0] for val in segments_torques]
	
	return (segments_lengths, segments_sigmas, segments_torques, segments_dna_states, segments_writhe_fractions, segments_plectoneme_thresholds)

def get_RNAP_velocities(model: Model, state_vector: list[float], RNAP_gene_index: list[int], segments_lengths: list[float], segments_torques: list[float]) -> list[float]: # Get the velocities of all RNAPs based on the segments lengths and torques
	RNAP_count = len(RNAP_gene_index)
	if RNAP_count == 0: # No RNAPs on the DNA; return empty list
		return []

	sorted_positions, sorted_ids = get_ordering_of_RNAPs_and_proteins(model, RNAP_gene_index, state_vector, topological_barrier_proteins_only = True)	
	RNAP_velocities = []
	left_segment_length = 0.0
	right_segment_length = 0.0
	left_torque = 0.0
	right_torque = 0.0
	curr_RNAP_index = 0
	for i in range(len(sorted_ids)):
		if 'RNAP' not in sorted_ids[i]:
			continue
		left_segment_length = segments_lengths[i + 1]
		right_segment_length = segments_lengths[i]
		left_torque = segments_torques[i + 1]
		right_torque = segments_torques[i]

		RNAP_velocities.append(get_RNAP_velocity(model, RNAP_gene_index[curr_RNAP_index], left_segment_length, right_segment_length, left_torque, right_torque))
		curr_RNAP_index += 1
	
	assert len(RNAP_velocities) == RNAP_count, 'Number of calculated RNAP velocities does not match number of RNAPs on the DNA.'
	
	sorted_positions, sorted_ids = get_ordering_of_RNAPs_and_proteins(model, RNAP_gene_index, state_vector)
	for i in range(RNAP_count):
		left_obstacle_position, left_obstacle_id, right_obstacle_position, right_obstacle_id = get_nearest_steric_obstacles_for_RNAP(i, state_vector[i], sorted_positions, sorted_ids)
		if model.genomic_setup.gene_directions[RNAP_gene_index[i]] == 1:
			if right_obstacle_position is not None:
				steric_hindrance_distance = model.model_setup.RNAP_diameter / 2.0
				if 'RNAP' in right_obstacle_id:
					steric_hindrance_distance += model.model_setup.RNAP_diameter / 2.0
				elif model.binding_proteins[int(right_obstacle_id)].is_a_nucleosome:
					steric_hindrance_distance += (model.genomic_setup.per_nucleosome_DNA_length + model.genomic_setup.nucleosome_linker_length) / 2.0
				else:
					steric_hindrance_distance += model.model_setup.generic_binding_protein_diameter / 2.0
				RNAP_velocities[i] = RNAP_velocities[i]*get_steric_hindrance_factor(model, right_obstacle_position - state_vector[i], steric_hindrance_distance)
		else:
			if left_obstacle_position is not None:
				steric_hindrance_distance = model.model_setup.RNAP_diameter / 2.0
				if 'RNAP' in left_obstacle_id:
					steric_hindrance_distance += model.model_setup.RNAP_diameter / 2.0
				elif model.binding_proteins[int(left_obstacle_id)].is_a_nucleosome:
					steric_hindrance_distance += (model.genomic_setup.per_nucleosome_DNA_length + model.genomic_setup.nucleosome_linker_length) / 2.0
				else:
					steric_hindrance_distance += model.model_setup.generic_binding_protein_diameter / 2.0
				RNAP_velocities[i] = RNAP_velocities[i]*get_steric_hindrance_factor(model, state_vector[i] - left_obstacle_position, steric_hindrance_distance)
	
	return RNAP_velocities

def get_RNAP_angular_velocities(model: Model, RNAP_gene_index: list[int], state_vector: list[float], segments_lengths: list[float], segments_torques: list[float], RNAP_velocities: list[float]) -> list[float]: # Get the angular velocities of all RNAPs based on the state_vector, segments lengths and torques, and RNAP linear velocities
	RNAP_count = len(RNAP_gene_index)
	if RNAP_count == 0: # No RNAPs on the DNA; return empty list
		return []
	
	x_vector = state_vector[0:RNAP_count]
	sorted_positions, sorted_ids = get_ordering_of_RNAPs_and_proteins(model, RNAP_gene_index, state_vector, topological_barrier_proteins_only = True)	
	
	RNAP_angular_velocities = []
	left_torque = 0.0
	right_torque = 0.0
	curr_RNAP_index = 0
	for i in range(len(sorted_ids)):
		if 'RNAP' not in sorted_ids[i]:
			continue
		left_torque = segments_torques[i + 1]
		right_torque = segments_torques[i]

		RNAP_angular_velocities.append(get_RNAP_angular_velocity(model, RNAP_gene_index[curr_RNAP_index], x_vector[curr_RNAP_index], RNAP_velocities[curr_RNAP_index], left_torque, right_torque))
		curr_RNAP_index += 1
	
	return RNAP_angular_velocities

def get_segments_Lk_dynamics(model: Model, RNAP_gene_index: list[int], state_vector: list[float], dx_dt: list[float], dtheta_dt: list[float], segments_lengths: list[float], segments_sigmas: list[float], segments_torques: list[float], segments_writhe_fractions: list[float]) -> list[float]: # Get the rates of change of linking number for all DNA segments based on the RNAP angular velocities
	sorted_positions, sorted_ids = get_ordering_of_RNAPs_and_proteins(model, RNAP_gene_index, state_vector, topological_barrier_proteins_only = True)
	assert len(segments_lengths) == len(segments_sigmas) == len(segments_torques) == len(segments_writhe_fractions) == len(sorted_positions) + 1, 'Length of segments attributes lists should be one more than the number of boundaries (RNAPs and topological barrier proteins) on the DNA.'

	if len(sorted_positions) == 0: # No RNAPs or topological barrier proteins on the DNA; single segment between clamps
		return [0.0]
	
	dLk_dt = []
	dx_dt_front = 0.0
	dx_dt_back = 0.0
	dtheta_dt_front = 0.0
	dtheta_dt_back = 0.0
	is_rightmost_segment = False
	is_leftmost_segment = False
	curr_RNAP_index = 0
	left_barrier_type = None
	right_barrier_type = None
	for i in range(len(segments_lengths)):
		if i == 0: # Rightmost segment
			if 'RNAP' in sorted_ids[i]:
				dx_dt_front = 0.0
				dx_dt_back = dx_dt[curr_RNAP_index]
				dtheta_dt_front = 0.0
				dtheta_dt_back = dtheta_dt[curr_RNAP_index]
				curr_RNAP_index += 1
				left_barrier_type = 'RNAP'
				right_barrier_type = 'clamp'
			else:
				dx_dt_front = 0.0
				dx_dt_back = 0.0
				dtheta_dt_front = 0.0
				dtheta_dt_back = 0.0
				left_barrier_type = sorted_ids[i]
				right_barrier_type = 'clamp'
			is_rightmost_segment = True
			is_leftmost_segment = False
		elif i == len(segments_lengths) - 1: # Leftmost segment
			if 'RNAP' in sorted_ids[i - 1]:
				dx_dt_front = dx_dt[curr_RNAP_index - 1]
				dx_dt_back = 0.0
				dtheta_dt_front = dtheta_dt[curr_RNAP_index - 1]
				dtheta_dt_back = 0.0
				left_barrier_type = 'clamp'
				right_barrier_type = 'RNAP'
			else:
				dx_dt_front = 0.0
				dx_dt_back = 0.0
				dtheta_dt_front = 0.0
				dtheta_dt_back = 0.0
				left_barrier_type = 'clamp'
				right_barrier_type = sorted_ids[i - 1]
			is_rightmost_segment = False
			is_leftmost_segment = True
		else:
			if 'RNAP' in sorted_ids[i - 1]:
				dx_dt_front = dx_dt[curr_RNAP_index - 1]
				dtheta_dt_front = dtheta_dt[curr_RNAP_index - 1]
				left_barrier_type = 'RNAP'
			else:
				dx_dt_front = 0.0
				dtheta_dt_front = 0.0
				left_barrier_type = sorted_ids[i - 1]
			if 'RNAP' in sorted_ids[i]:
				dx_dt_back = dx_dt[curr_RNAP_index]
				dtheta_dt_back = dtheta_dt[curr_RNAP_index]
				curr_RNAP_index += 1
				right_barrier_type = 'RNAP'
			else:
				dx_dt_back = 0.0
				dtheta_dt_back = 0.0
				right_barrier_type = sorted_ids[i]
			is_rightmost_segment = False
			is_leftmost_segment = False
		dLk_dt.append(get_segment_Lk_dynamics(model, dx_dt_front, dx_dt_back, dtheta_dt_front, dtheta_dt_back, left_barrier_type, right_barrier_type, is_rightmost_segment, is_leftmost_segment))
	
	return dLk_dt

def get_RNAP_recruitment_rates(model: Model, RNAP_gene_index: list[int], state_vector: list[float], segments_lengths: list[float], segments_sigmas: list[float]) -> list[float]: # Get the RNAP recruitment rates at all TSSes based on the current state_vector and segments attributes
	RNAP_count = len(RNAP_gene_index)
	x_vector = state_vector[0:RNAP_count]

	TSS_segments_indices = [get_spot_segment_index(model.genomic_setup.TSSes[i], segments_lengths) for i in range(len(model.genomic_setup.gene_names))] # Get the segment indices for all TSSes

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

def get_mRNA_degradation_rates(model: Model) -> list[float]: # Get the mRNA degradation rates for all genes based on their current mRNA counts
	mRNA_degradation_rates = []
	for i in range(len(model.genomic_setup.gene_names)):
		mRNA_degradation_rates.append(get_mRNA_degradation_rate(model, model.mRNA_counts[i]))
	return mRNA_degradation_rates

def get_binding_proteins_on_rates(model: Model, segments_lengths: list[float], segments_sigmas: list[float]) -> list[list[float]]:
	binding_proteins_on_rates = []
	for i in range(len(model.binding_proteins)):
		protein = model.binding_proteins[i]
		bound_protein_count = len(model.binding_proteins_positions[i])
		unbound_protein_count = protein.total_copy_number - bound_protein_count
		per_segment_on_rates = [protein.on_rate_func(segment_length, segment_sigma) for segment_length, segment_sigma in zip(segments_lengths, segments_sigmas)]
		binding_proteins_on_rates.append([unbound_protein_count*per_segment_on_rate for per_segment_on_rate in per_segment_on_rates])
	
	return binding_proteins_on_rates

def get_binding_proteins_off_rates(model: Model, segments_lengths: list[float], segments_sigmas: list[float]) -> list[list[float]]:
	binding_proteins_off_rates = []
	for i in range(len(model.binding_proteins)):
		protein = model.binding_proteins[i]
		segment_of_bound_proteins = [get_spot_segment_index(model.binding_proteins_positions[i][j], segments_lengths) for j in range(len(model.binding_proteins_positions[i]))]
		per_segment_off_rates = [protein.off_rate_func(segment_length, segment_sigma) for segment_length, segment_sigma in zip(segments_lengths, segments_sigmas)]
		binding_proteins_off_rates.append([per_segment_off_rates[segment_index] for segment_index in segment_of_bound_proteins])
	
	return binding_proteins_off_rates

def update_Lk_vector_after_RNAP_or_protein_recruitment(model: Model, recruitment_position: float, RNAP_gene_index: list[int], state_vector: list[float], segments_lengths: list[float], segments_sigmas: list[float]) -> None: # Update the model's Lk vector after an RNAP or a topological barrier protein is recruited at the given position; Lk for the segments on either side of the recruitment position are calculated such that the supercoiling density in the segments is the same as in the segment spanning the recruitment position before recruitment
	segment_index = get_spot_segment_index(recruitment_position, segments_lengths)
	segment_sigma = segments_sigmas[segment_index]

	x_vector = list(state_vector[0:len(RNAP_gene_index)])
	topological_barrier_proteins_positions = []
	for i in range(len(model.binding_proteins)):
		if model.binding_proteins[i].is_topological_barrier:
			topological_barrier_proteins_positions = topological_barrier_proteins_positions + model.binding_proteins_positions[i]
	x_vector = x_vector + topological_barrier_proteins_positions
	x_vector = sorted(x_vector, reverse = True)

	left_segment_boundary = model.genomic_setup.clamp_left
	right_segment_boundary = model.genomic_setup.clamp_right

	if len(segments_lengths) > 1:
		if segment_index == 0:
			left_segment_boundary = x_vector[0]
		elif segment_index == len(segments_lengths) - 1:
			right_segment_boundary = x_vector[-1]
		else:
			left_segment_boundary = x_vector[segment_index]
			right_segment_boundary = x_vector[segment_index - 1]
	
	left_segment_length = recruitment_position - left_segment_boundary
	right_segment_length = right_segment_boundary - recruitment_position
	assert left_segment_length >= 0.0 and right_segment_length >= 0.0, 'Error in calculating segment lengths to the left and right of the RNAP / protein recruitment position.'

	left_Lk0 = left_segment_length / model.model_setup.h_dna
	right_Lk0 = right_segment_length / model.model_setup.h_dna

	# Based on sigma = (Lk - Lk0) / Lk0 => Lk = Lk0 * (1 + sigma)
	left_Lk = left_Lk0*(1.0 + segment_sigma)
	right_Lk = right_Lk0*(1.0 + segment_sigma)

	model.Lk = model.Lk[:segment_index] + [right_Lk, left_Lk] + model.Lk[segment_index + 1:]

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
			model.mRNA_counts[RNAP_gene_index[i]] += 1
			simulation_setup_and_state.RNAPs_finished_transcription[RNAP_gene_index[i]] += 1
			simulation_setup_and_state.RNAPs_exit_positions[RNAP_gene_index[i]].append(x_vector[i])
			simulation_setup_and_state.RNAP_exit_times[RNAP_gene_index[i]].append(t)

	new_RNAP_gene_index = [RNAP_gene_index[i] for i in range(RNAP_count) if RNAPs_alive_status[i] == 1]
	new_x_vector = [x_vector[i] for i in range(RNAP_count) if RNAPs_alive_status[i] == 1]
	new_theta_vector = [theta_vector[i] for i in range(RNAP_count) if RNAPs_alive_status[i] == 1]

	# When an RNAP finishes transcription, the two segments on either side of it merge into one segment; thus, we need to update the Lk_vector accordingly. Lk of the merged segment is the sum of the Lk of the two segments.
	sorted_positions, sorted_ids = get_ordering_of_RNAPs_and_proteins(model, RNAP_gene_index, state_vector, topological_barrier_proteins_only = True)
	new_Lk_vector = []
	curr_barrier_index = 0
	curr_RNAP_index = 0
	while curr_barrier_index < len(sorted_ids):
		Lk_front_segment = Lk_vector[curr_barrier_index]
		Lk_back_segment = Lk_vector[curr_barrier_index + 1]
		if 'RNAP' not in sorted_ids[curr_barrier_index]: # Current barrier is a topological barrier protein; just add the Lk of the front segment and move to the next barrier
			new_Lk_vector.append(Lk_front_segment)
			curr_barrier_index += 1
			continue
		if RNAPs_alive_status[curr_RNAP_index] == 1: # Current barrier is an alive RNAP; just add the Lk of the front segment and move to the next barrier and RNAP
			new_Lk_vector.append(Lk_front_segment)
			curr_barrier_index += 1
			curr_RNAP_index += 1
		else: # Current barrier is a dead RNAP; need to merge segments until the next alive RNAP or topological barrier protein and add the total Lk of the merged segment; then move to the next barrier and the next alive RNAP
			index_of_next_alive_RNAP_or_barrier = curr_barrier_index
			index_of_next_alive_RNAP = curr_RNAP_index
			next_alive_is_RNAP = False
			total_Lk_merged_segment = 0.0
			while index_of_next_alive_RNAP_or_barrier < len(sorted_ids):
				if 'RNAP' not in sorted_ids[index_of_next_alive_RNAP_or_barrier]: # Next barrier is a topological barrier protein; stop merging segments
					next_alive_is_RNAP = False
					break
				if RNAPs_alive_status[index_of_next_alive_RNAP] == 1: # Next alive RNAP found; stop merging segments
					next_alive_is_RNAP = True
					break
				Lk_front_segment = Lk_vector[index_of_next_alive_RNAP_or_barrier]
				Lk_back_segment = Lk_vector[index_of_next_alive_RNAP_or_barrier + 1]
				total_Lk_merged_segment += Lk_front_segment
				if 'RNAP' in sorted_ids[index_of_next_alive_RNAP_or_barrier]:
					index_of_next_alive_RNAP += 1
				index_of_next_alive_RNAP_or_barrier += 1
			total_Lk_merged_segment += Lk_vector[index_of_next_alive_RNAP_or_barrier] # Add the Lk of the segment after the last merged segment to the total Lk of the merged segment
			new_Lk_vector.append(total_Lk_merged_segment)
			curr_barrier_index = index_of_next_alive_RNAP_or_barrier + 1
			if next_alive_is_RNAP:
				curr_RNAP_index = index_of_next_alive_RNAP + 1
			else:
				curr_RNAP_index = index_of_next_alive_RNAP
	if 'RNAP' in sorted_ids[-1] and RNAPs_alive_status[-1] == 1:
		new_Lk_vector.append(Lk_vector[-1])
	if 'RNAP' not in sorted_ids[-1]:
		new_Lk_vector.append(Lk_vector[-1])
	
	state_vector[:] = new_x_vector + new_theta_vector + new_Lk_vector + [state_vector[-1]]
	RNAP_gene_index[:] = new_RNAP_gene_index

def update_Lk_vector_after_protein_unbinding(model: Model, unbinding_position: float, RNAP_gene_index: list[int], state_vector: list[float], segments_lengths: list[float], segments_sigmas: list[float]) -> None:
	x_vector = list(state_vector[0:len(RNAP_gene_index)])
	topological_barrier_proteins_positions = []
	for i in range(len(model.binding_proteins)):
		if model.binding_proteins[i].is_topological_barrier:
			topological_barrier_proteins_positions = topological_barrier_proteins_positions + model.binding_proteins_positions[i]
	x_vector = x_vector + topological_barrier_proteins_positions
	x_vector = sorted(x_vector, reverse = True)
	assert len(x_vector) + 1 == len(segments_lengths), 'Error: number of segments should be one more than the number of boundaries (RNAPs and topological barrier proteins) on the DNA.'

	protein_index_in_x_vector = -1
	for i in range(len(x_vector)):
		if abs(x_vector[i] - unbinding_position) < 1e-6: # Found the position of the unbinding protein in the x_vector; use this to determine the segment index for updating Lk
			protein_index_in_x_vector = i
			break
	assert protein_index_in_x_vector != -1, 'Error: unbinding position not found in x_vector when updating Lk vector after protein unbinding.'

	right_segment_index = protein_index_in_x_vector
	left_segment_index = protein_index_in_x_vector + 1

	model.Lk = model.Lk[:right_segment_index] + [model.Lk[right_segment_index] + model.Lk[left_segment_index]] + model.Lk[left_segment_index + 1:]

def get_events_rates(model: Model, RNAP_gene_index: list[int], state_vector: list[float]) -> tuple[list[float], list[int]]: # Get the rates of all possible events and the indices that separate different event types in the rates_vector
	segments_lengths, segments_sigmas, segments_torques, segments_dna_states, segments_writhe_fractions, segments_plectoneme_thresholds = calculate_segments_attributes(model, RNAP_gene_index, state_vector)

	RNAP_recruitment_rates = get_RNAP_recruitment_rates(model, RNAP_gene_index, state_vector, segments_lengths, segments_sigmas)
	model_observation_event_rate = [model.model_setup.model_observation_event_rate]
	global_supercoiling_relaxation_rate = [model.model_setup.global_supercoiling_relaxation_rate]
	local_supercoiling_relaxation_rates = model.model_setup.local_supercoiling_relaxation_rates
	TOPO_activity_rates = [model.model_setup.TOP1_effective_relaxation_rate, model.model_setup.TOP2_effective_relaxation_rate]
	mRNA_degradation_rates = get_mRNA_degradation_rates(model)
	binding_proteins_on_rates = [sum(per_protein_on_rates) for per_protein_on_rates in get_binding_proteins_on_rates(model, segments_lengths, segments_sigmas)]
	binding_proteins_off_rates = [sum(per_protein_off_rates) for per_protein_off_rates in get_binding_proteins_off_rates(model, segments_lengths, segments_sigmas)]
	promoter_on_rates = [get_promoter_on_rate(model, gene_index, segments_sigmas[get_spot_segment_index(model.genomic_setup.TSSes[gene_index], segments_lengths)]) for gene_index in range(len(model.genomic_setup.gene_names))]
	promoter_off_rates = [get_promoter_off_rate(model, gene_index, segments_sigmas[get_spot_segment_index(model.genomic_setup.TSSes[gene_index], segments_lengths)]) for gene_index in range(len(model.genomic_setup.gene_names))]

	rates_vector = RNAP_recruitment_rates + model_observation_event_rate + global_supercoiling_relaxation_rate + local_supercoiling_relaxation_rates + TOPO_activity_rates + mRNA_degradation_rates + binding_proteins_on_rates + binding_proteins_off_rates + promoter_on_rates + promoter_off_rates
	assert all(rate >= 0.0 for rate in rates_vector), 'Negative rate encountered in calculating events rates.'

	events_indices = []
	events_indices.append(len(RNAP_recruitment_rates))
	events_indices.append(events_indices[-1] + len(model_observation_event_rate))
	events_indices.append(events_indices[-1] + len(global_supercoiling_relaxation_rate))
	events_indices.append(events_indices[-1] + len(local_supercoiling_relaxation_rates))
	events_indices.append(events_indices[-1] + len(TOPO_activity_rates))
	events_indices.append(events_indices[-1] + len(mRNA_degradation_rates))
	events_indices.append(events_indices[-1] + len(binding_proteins_on_rates))
	events_indices.append(events_indices[-1] + len(binding_proteins_off_rates))
	events_indices.append(events_indices[-1] + len(promoter_on_rates))
	events_indices.append(events_indices[-1] + len(promoter_off_rates))

	return rates_vector, events_indices

def model_dynamics(t: float, state_vector: list[float], RNAP_gene_index: list[int], model: Model, simulation_setup_and_state: SimulationSetupAndState) -> list[float]: # Calculate the time derivatives of the state_vector and the cumulative propensity for event selection
	RNAP_count = len(RNAP_gene_index)

	segments_lengths, segments_sigmas, segments_torques, segments_dna_states, segments_writhe_fractions, segments_plectoneme_thresholds = calculate_segments_attributes(model, RNAP_gene_index, state_vector)

	rates_vector, _ = get_events_rates(model, RNAP_gene_index, state_vector)

	dx_dt = get_RNAP_velocities(model, state_vector, RNAP_gene_index, segments_lengths, segments_torques)
	dtheta_dt = get_RNAP_angular_velocities(model, RNAP_gene_index, state_vector, segments_lengths, segments_torques, dx_dt)
	dLk_dt = get_segments_Lk_dynamics(model, RNAP_gene_index, state_vector, dx_dt, dtheta_dt, segments_lengths, segments_sigmas, segments_torques, segments_writhe_fractions)

	return dx_dt + dtheta_dt + dLk_dt + [sum(rates_vector)]

def integrate(model: Model, simulation_setup_and_state: SimulationSetupAndState, t_start: float, state_vector: list[float], RNAP_gene_index: list[int], p0: float, print_at_each_integration_step: Union[Callable, None]) -> tuple[float, float]: # Integrate the model dynamics until the next event occurs; return the time step dt and cumulative propensity a0 at event time
	dt = simulation_setup_and_state.RNAP_alive_status_check_interval
	t = t_start
	a0_ini = state_vector[-1]
	t_event = -1.0
	t_event_index = -1

	while True:
		if print_at_each_integration_step is not None:
			print_at_each_integration_step(model, simulation_setup_and_state, t, state_vector)
		t_eval = []
		t_first = t
		while t_first < t + dt:
			t_eval.append(t_first)
			t_first += simulation_setup_and_state.integration_time_resolution
		t_eval.append(t + dt)
		sol = solve_ivp(lambda t, y: model_dynamics(t, y, RNAP_gene_index, model, simulation_setup_and_state), (t, t + dt), state_vector, t_eval = t_eval, method = simulation_setup_and_state.integration_method, rtol = simulation_setup_and_state.integration_rtol, atol = simulation_setup_and_state.integration_atol)
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
			# assert len(state_vector) == len(RNAP_gene_index) + len(RNAP_gene_index) + len(RNAP_gene_index) + 1 + 1, 'State vector length mismatch after integration and removing dead RNAPs (when event occurred during integration).'
			return (dt_current, a0)
		# No event has occurred; update state_vector and continue integration
		state_vector[:] = sol.y[:, -1]
		update_state_vector_to_remove_dead_RNAPs(model, RNAP_gene_index, t + dt, state_vector, simulation_setup_and_state)
		# assert len(state_vector) == len(RNAP_gene_index) + len(RNAP_gene_index) + len(RNAP_gene_index) + 1 + 1, 'State vector length mismatch after integration and removing dead RNAPs (when no event occurred during integration).'
		t += dt