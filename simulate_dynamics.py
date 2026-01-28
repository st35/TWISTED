from model_setup import *
from model_dynamics import *
from utilities import *

from typing import Callable
from random import random
from math import log

def simulate_dynamics(model: Model, simulation_setup_and_state: SimulationSetupAndState, print_at_each_integration_step: Union[Callable, None] = None, print_at_each_simulation_step: Union[Callable, None] = None, print_at_end_of_simulation: Union[Callable, None] = None) -> None: # Main function to simulate the dynamics of the model until the end condition is met
	while True:
		if print_at_each_simulation_step is not None:
			print_at_each_simulation_step(model, simulation_setup_and_state) # Print current state if a print function is provided

		RNAP_gene_index, state_vector = get_state_vectors_from_dicts(model) # Get current state vectors from model dictionaries

		p0 = random() # Generate a random number for event time calculation
		dt, _ = integrate(model, simulation_setup_and_state, simulation_setup_and_state.curr_simulation_time, state_vector, RNAP_gene_index, p0, print_at_each_integration_step) # Integrate dynamics until next event
		assert dt is not None, 'Integration failed to return a valid time step.'
		update_dicts_from_state_vector(model, RNAP_gene_index, state_vector) # Update model dictionaries from updated state vectors
		simulation_setup_and_state.curr_simulation_time += dt # Update current simulation time
		segments_lengths, segments_sigmas, segments_torques, segments_dna_states, segments_writhe_fractions = calculate_segments_attributes(model, RNAP_gene_index, state_vector) # Recalculate segment attributes after integration

		rates_vector, events_indices = get_events_rates(model, RNAP_gene_index, state_vector) # Get event rates and event indices
		p1 = random() # Generate a random number for event selection
		event_index = select_event_based_on_propensities(rates_vector, p1) # Select event based on propensities

		if event_index < events_indices[0]: # RNAP recruitment event
			event = event_index - 0
			if get_TSS_steric_hindrance_status(model, model.genomic_setup.TSSes[event], RNAP_gene_index, state_vector) == 1: # Steric hindrance at TSS; recruitment fails
				pass
			else:
				event_gene_index = event
				model.x_dict[event_gene_index].append(model.genomic_setup.TSSes[event_gene_index])
				model.theta_dict[event_gene_index].append(0.0)
				update_Lk_vector_after_RNAP_recruitment(model, event_gene_index, RNAP_gene_index, state_vector, segments_lengths, segments_sigmas) # Update linking number vector after RNAP recruitment
				simulation_setup_and_state.RNAP_recruitment_times[event_gene_index].append(simulation_setup_and_state.curr_simulation_time)
		elif event_index < events_indices[1]: # Model observation event: a dummy event to ensure simulation time progresses even if all biological events rates are zero
			pass
		elif event_index < events_indices[2]: # Global supercoiling relaxation event: supercoiling_relaxation_dynamics_mode in ['global_overall', 'global_per_segment']
			segments_Lk0 = [segments_length / model.model_setup.h_dna for segments_length in segments_lengths]
			if model.model_setup.supercoiling_relaxation_dynamics_mode == 'global_overall':
				model.Lk = [Lk0 for Lk0 in segments_Lk0]
			elif model.model_setup.supercoiling_relaxation_dynamics_mode == 'global_per_segment':
				per_segment_propensity = [segments_length / (model.genomic_setup.clamp_right - model.genomic_setup.clamp_left) for segments_length in segments_lengths]
				chosen_segment_index = select_event_based_on_propensities(per_segment_propensity, random())
				if chosen_segment_index is not None:
					model.Lk[chosen_segment_index] = segments_Lk0[chosen_segment_index]
			else:
				raise ValueError('Invalid type for global supercoiling relaxation.')
		elif event_index < events_indices[3]: # Type-specific supercoiling relaxation event: supercoiling_relaxation_dynamics_mode in ['global_by_type', 'per_segment_by_type']
			event = event_index - events_indices[2]
			segments_Lk0 = [segments_length / model.model_setup.h_dna for segments_length in segments_lengths]
			if model.model_setup.supercoiling_relaxation_dynamics_mode == 'global_by_type':
				if event == 0: # Relax positive supercoiling only
					model.Lk = [segments_Lk0[i] if segments_sigmas[i] > 0.0 else model.Lk[i] for i in range(len(segments_lengths))]
				elif event == 1: # Relax negative supercoiling only
					model.Lk = [segments_Lk0[i] if segments_sigmas[i] < 0.0 else model.Lk[i] for i in range(len(segments_lengths))]
			else:
				chosen_segment_index = None
				if event == 0: # Relax positive supercoiling only
					per_segment_propensity = [segments_lengths[i] / (model.genomic_setup.clamp_right - model.genomic_setup.clamp_left) if segments_sigmas[i] > 0.0 else 0.0 for i in range(len(segments_lengths))]
					chosen_segment_index = select_event_based_on_propensities(per_segment_propensity, random())
				elif event == 1: # Relax negative supercoiling only
					per_segment_propensity = [segments_lengths[i] / (model.genomic_setup.clamp_right - model.genomic_setup.clamp_left) if segments_sigmas[i] < 0.0 else 0.0 for i in range(len(segments_lengths))]
					chosen_segment_index = select_event_based_on_propensities(per_segment_propensity, random())
				if chosen_segment_index is not None:
					model.Lk[chosen_segment_index] = segments_Lk0[chosen_segment_index]
		elif event_index < events_indices[4]: # Topoisomerase-mediated supercoiling relaxation: TOP1 / TOP2 binding event
			event = event_index - events_indices[3]
			if model.topoisomerase_status[event] == 1: # Topoisomerase is already bound; cannot bind again
				raise ValueError('Binding event selected for already bound topoisomerase.')
			else:
				if model.topoisomerase_type[event] == 0:
					TOP1_on_rates_per_segment = get_per_TOP1_binding_rate_for_each_segment(model, segments_lengths, segments_sigmas)
					chosen_segment_index = select_event_based_on_propensities(TOP1_on_rates_per_segment, random())
					binding_position = model.genomic_setup.clamp_left + sum(segments_lengths[:chosen_segment_index]) + (segments_lengths[chosen_segment_index]*uniform_random_in_interval(0.0, 1.0))
					if is_TOPO_binding_blocked(model, state_vector, binding_position) == 0:
						model.topoisomerase_status[event] = 1 # Bind topoisomerase
						model.topoisomerase_segment_indices[event] = chosen_segment_index
						model.topoisomerase_positions[event] = binding_position
				else:
					TOP2_on_rates_per_segment = get_per_TOP2_binding_rate_for_each_segment(model, segments_lengths, segments_sigmas)
					chosen_segment_index = select_event_based_on_propensities(TOP2_on_rates_per_segment, random())
					binding_position = model.genomic_setup.clamp_left + sum(segments_lengths[:chosen_segment_index]) + (segments_lengths[chosen_segment_index]*uniform_random_in_interval(0.0, 1.0))
					if is_TOPO_binding_blocked(model, state_vector, binding_position) == 0:
						model.topoisomerase_status[event] = 1 # Bind topoisomerase
						model.topoisomerase_segment_indices[event] = chosen_segment_index
						model.topoisomerase_positions[event] = binding_position
		elif event_index < events_indices[5]: # Topoisomerase-mediated supercoiling relaxation: TOP1 / TOP2 unbinding event
			event = event_index - events_indices[4]
			if model.topoisomerase_status[event] == 0: # Topoisomerase is already unbound; cannot unbind again
				raise ValueError('Unbinding event selected for already unbound topoisomerase.')
			else:
				model.topoisomerase_status[event] = 0 # Unbind topoisomerase
				model.topoisomerase_segment_indices[event] = -1
				model.topoisomerase_positions[event] = -1.0
		else:
			raise ValueError('Event index out of bounds during simulation.')
		
		if simulation_setup_and_state.simulation_end_mode == 0: # End simulation based on time
			if simulation_setup_and_state.curr_simulation_time >= simulation_setup_and_state.simulation_end_time:
				simulation_setup_and_state.simulation_completed = True
				break
		else: # End simulation based on number of RNAPs that have completed transcription for each gene
			all_genes_completed = True
			for gene_index in range(len(model.genomic_setup.gene_names)):
				if simulation_setup_and_state.RNAPs_finished_transcription[gene_index] < simulation_setup_and_state.simulation_end_event_counts[gene_index]:
					all_genes_completed = False
					break
			if all_genes_completed:
				simulation_setup_and_state.simulation_completed = True
				break
	if print_at_end_of_simulation is not None:
		print_at_end_of_simulation(model, simulation_setup_and_state) # Print final state / outcome if a print function is provided