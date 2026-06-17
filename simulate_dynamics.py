from model_setup import *
from model_dynamics import *
from utilities import *

from typing import Callable
from random import random
from math import log

def simulate_dynamics(model: Model, simulation_setup_and_state: SimulationSetupAndState, print_at_each_integration_step: Union[Callable, None] = None, print_at_each_simulation_step: Union[Callable, None] = None, print_at_end_of_simulation: Union[Callable, None] = None) -> None: # Main function to simulate the dynamics of the model until the end condition is met
	if not simulation_setup_and_state.state_has_been_initialized:
		simulation_setup_and_state.setup_simulation_state(model.genomic_setup) # Set up initial simulation state based on genomic setup
	while True:
		if print_at_each_simulation_step is not None:
			print_at_each_simulation_step(model, simulation_setup_and_state) # Print current state if a print function is provided

		RNAP_gene_index, state_vector = get_state_vectors_from_dicts(model) # Get current state vectors from model dictionaries

		p0 = random() # Generate a random number for event time calculation
		dt, _ = integrate(model, simulation_setup_and_state, simulation_setup_and_state.curr_simulation_time, state_vector, RNAP_gene_index, p0, print_at_each_integration_step) # Integrate dynamics until next event
		assert dt is not None, 'Integration failed to return a valid time step.'
		update_dicts_from_state_vector(model, RNAP_gene_index, state_vector) # Update model dictionaries from updated state vectors
		simulation_setup_and_state.curr_simulation_time += dt # Update current simulation time
		segments_lengths, segments_sigmas, segments_torques, segments_dna_states, segments_writhe_fractions, segments_plectoneme_thresholds = calculate_segments_attributes(model, RNAP_gene_index, state_vector) # Recalculate segment attributes after integration

		rates_vector, events_indices = get_events_rates(model, RNAP_gene_index, state_vector) # Get event rates and event indices
		p1 = random() # Generate a random number for event selection
		event_index = select_event_based_on_propensities(rates_vector, p1) # Select event based on propensities
		simulation_setup_and_state.last_event_index = event_index # Update the last event index

		if event_index < events_indices[0]: # RNAP recruitment event
			simulation_setup_and_state.last_event_type = 'RNAP_recruitment'
			event = event_index - 0
			TSS_steric_hindrance_status, blocking_entity_position, blocking_entity_id = get_TSS_steric_hindrance_status(model, model.genomic_setup.TSSes[event], RNAP_gene_index, state_vector) # Check for steric hindrance at TSS
			if TSS_steric_hindrance_status == 1 and blocking_entity_id == -1: # Steric hindrance from another RNAP; recruitment fails
				pass
			elif TSS_steric_hindrance_status == 1 and blocking_entity_id >= 0 and blocking_entity_id < len(model.binding_proteins) and model.binding_proteins[blocking_entity_id].can_be_displaced_at_TSS_by_RNAP is False: # Steric hindrance from a bound protein that cannot be displaced; recruitment fails
				pass
			elif simulation_setup_and_state.max_RNAPs_to_recruit is not None and len(simulation_setup_and_state.RNAP_recruitment_times[event]) >= simulation_setup_and_state.max_RNAPs_to_recruit[event]: # Maximum number of RNAPs already recruited for this gene; recruitment fails
				pass
			else:
				event_gene_index = event
				model.x_dict[event_gene_index].append(model.genomic_setup.TSSes[event_gene_index])
				model.theta_dict[event_gene_index].append(0.0)
				update_Lk_vector_after_RNAP_or_protein_recruitment(model, model.genomic_setup.TSSes[event_gene_index], RNAP_gene_index, state_vector, segments_lengths, segments_sigmas) # Update linking number vector after RNAP recruitment
				simulation_setup_and_state.RNAP_recruitment_times[event_gene_index].append(simulation_setup_and_state.curr_simulation_time)
				if TSS_steric_hindrance_status == 1:
					assert blocking_entity_id >= 0 and blocking_entity_id < len(model.binding_proteins) and model.binding_proteins[blocking_entity_id].can_be_displaced_at_TSS_by_RNAP is True, 'There must be a bound protein that can be displaced if steric hindrance is present but recruitment is successful.'
					find_and_remove_from_list(model.binding_proteins_positions[blocking_entity_id], blocking_entity_position) # Displace the bound protein causing steric hindrance by removing it from the list of bound positions for that protein
				simulation_setup_and_state.last_event_type = 'RNAP_recruitment_successful'
		elif event_index < events_indices[1]: # Model observation event: a dummy event to ensure simulation time progresses even if all biological events rates are zero
			simulation_setup_and_state.last_event_type = 'model_observation'
		elif event_index < events_indices[2]: # Global supercoiling relaxation event: supercoiling_relaxation_dynamics_mode in ['global_overall', 'global_per_segment']
			segments_Lk0 = [segments_length / model.model_setup.h_dna for segments_length in segments_lengths]
			if model.model_setup.supercoiling_relaxation_dynamics_mode == 'global_overall':
				model.Lk = [Lk0 for Lk0 in segments_Lk0]
				simulation_setup_and_state.last_event_type = 'global_supercoiling_relaxation_overall'
			elif model.model_setup.supercoiling_relaxation_dynamics_mode == 'global_per_segment':
				per_segment_propensity = [segments_length / (model.genomic_setup.clamp_right - model.genomic_setup.clamp_left) for segments_length in segments_lengths]
				chosen_segment_index = select_event_based_on_propensities(per_segment_propensity, random())
				if chosen_segment_index is not None:
					model.Lk[chosen_segment_index] = segments_Lk0[chosen_segment_index]
					simulation_setup_and_state.last_event_type = 'global_supercoiling_relaxation_per_segment'
			else:
				raise ValueError('Invalid type for global supercoiling relaxation.')
		elif event_index < events_indices[3]: # Type-specific supercoiling relaxation event: supercoiling_relaxation_dynamics_mode in ['global_by_type', 'per_segment_by_type']
			event = event_index - events_indices[2]
			segments_Lk0 = [segments_length / model.model_setup.h_dna for segments_length in segments_lengths]
			if model.model_setup.supercoiling_relaxation_dynamics_mode == 'global_by_type':
				if event == 0: # Relax positive supercoiling only
					model.Lk = [segments_Lk0[i] if segments_sigmas[i] > 0.0 else model.Lk[i] for i in range(len(segments_lengths))]
					simulation_setup_and_state.last_event_type = 'global_supercoiling_relaxation_positive_only'
				elif event == 1: # Relax negative supercoiling only
					model.Lk = [segments_Lk0[i] if segments_sigmas[i] < 0.0 else model.Lk[i] for i in range(len(segments_lengths))]
					simulation_setup_and_state.last_event_type = 'global_supercoiling_relaxation_negative_only'
			else:
				chosen_segment_index = None
				if event == 0: # Relax positive supercoiling only
					per_segment_propensity = [segments_lengths[i] / (model.genomic_setup.clamp_right - model.genomic_setup.clamp_left) if segments_sigmas[i] > 0.0 else 0.0 for i in range(len(segments_lengths))]
					chosen_segment_index = select_event_based_on_propensities(per_segment_propensity, random())
					simulation_setup_and_state.last_event_type = 'per_segment_supercoiling_relaxation_positive_only'
				elif event == 1: # Relax negative supercoiling only
					per_segment_propensity = [segments_lengths[i] / (model.genomic_setup.clamp_right - model.genomic_setup.clamp_left) if segments_sigmas[i] < 0.0 else 0.0 for i in range(len(segments_lengths))]
					chosen_segment_index = select_event_based_on_propensities(per_segment_propensity, random())
					simulation_setup_and_state.last_event_type = 'per_segment_supercoiling_relaxation_negative_only'
				if chosen_segment_index is not None:
					model.Lk[chosen_segment_index] = segments_Lk0[chosen_segment_index]
		elif event_index < events_indices[4]: # Topoisomerase-mediated supercoiling relaxation; approximating TOP1 / TOP2 activity
			event = event_index - events_indices[3] # event = 0 for TOP1-mediated relaxation, event = 1 for TOP2-mediated relaxation
			segments_Lk0 = [segments_length / model.model_setup.h_dna for segments_length in segments_lengths]
			per_segment_propensity = [segments_lengths[i] / (model.genomic_setup.clamp_right - model.genomic_setup.clamp_left) for i in range(len(segments_lengths))]
			chosen_segment_index = select_event_based_on_propensities(per_segment_propensity, random())
			if chosen_segment_index is not None:
				if event == 0: # TOP1 activity: relax supercoiling to a relaxed state if no writhe, otherwise TOP1 cannot act
					simulation_setup_and_state.last_event_type = 'TOP1_supercoiling_relaxation_approximation'
					if segments_writhe_fractions[chosen_segment_index] > 0.0: # Non-zero writhe; TOP1 cannot act
						pass
					else:
						model.Lk[chosen_segment_index] = segments_Lk0[chosen_segment_index] # Relax supercoiling to a relaxed state
				elif event == 1: # TOP2 activity: relax supercoiling if writhe is present, otherwise TOP2 cannot act
					simulation_setup_and_state.last_event_type = 'TOP2_supercoiling_relaxation_approximation'
					if segments_writhe_fractions[chosen_segment_index] > 0.0: # Non-zero writhe; TOP2 can act:
						model.Lk[chosen_segment_index] = segments_Lk0[chosen_segment_index]*(1.0 + segments_plectoneme_thresholds[chosen_segment_index]) # Relax supercoiling to the threshold beyond which plectonemes form, since TOP2 relaxes supercoiling only until writhe is removed
		elif event_index < events_indices[5]: # mRNA degradation event
			simulation_setup_and_state.last_event_type = 'mRNA_degradation'
			event = event_index - events_indices[4]
			if model.mRNA_counts[event] <= 0: # No mRNA to degrade; cannot degrade
				raise ValueError('mRNA degradation event selected for gene with zero mRNA count.')
			else:
				model.mRNA_counts[event] -= 1 # Degrade one mRNA molecule for the gene
		elif event_index < events_indices[6]: # Binding protein binding event
			event = event_index - events_indices[5]
			simulation_setup_and_state.last_event_type = model.binding_proteins[event].protein_name + '_binding'
			per_segment_on_rates = get_binding_proteins_on_rates(model, segments_lengths, segments_sigmas)[event]
			chosen_segment_index = select_event_based_on_propensities(per_segment_on_rates, random())
			binding_position = model.genomic_setup.clamp_left + sum(segments_lengths[chosen_segment_index + 1:]) + (segments_lengths[chosen_segment_index]*uniform_random_in_interval(0.0, 1.0))
			if is_protein_binding_blocked(model, RNAP_gene_index, state_vector, event, binding_position) == 0:
				if model.binding_proteins[event].is_topological_barrier:
					update_Lk_vector_after_RNAP_or_protein_recruitment(model, binding_position, RNAP_gene_index, state_vector, segments_lengths, segments_sigmas) # Update linking number vector after binding of a topological barrier protein
				model.binding_proteins_positions[event].append(binding_position) # Bind the binding protein at the chosen position
		elif event_index < events_indices[7]: # Binding protein unbinding event
			event = event_index - events_indices[6]
			simulation_setup_and_state.last_event_type = model.binding_proteins[event].protein_name + '_unbinding'
			binding_proteins_off_rates = get_binding_proteins_off_rates(model, segments_lengths, segments_sigmas)[event]
			assert len(binding_proteins_off_rates) == len(model.binding_proteins_positions[event]), 'Length of binding proteins off rates does not match number of bound proteins for the event.'
			chosen_bound_protein_index = select_event_based_on_propensities(binding_proteins_off_rates, random())
			if chosen_bound_protein_index is not None:
				if model.binding_proteins[event].is_topological_barrier:
					update_Lk_vector_after_protein_unbinding(model, model.binding_proteins_positions[event][chosen_bound_protein_index], RNAP_gene_index, state_vector, segments_lengths, segments_sigmas) # Update linking number vector after unbinding of a topological barrier protein
				model.binding_proteins_positions[event].pop(chosen_bound_protein_index) # Unbind the chosen binding protein
		elif event_index < events_indices[8]: # Promoter ON event
			simulation_setup_and_state.last_event_type = 'promoter_ON'
			event = event_index - events_indices[7]
			if model.promoter_status[event] == 1: # Promoter already ON; cannot turn ON again
				raise ValueError('Promoter ON event selected for a promoter that is already ON.')
			else:
				model.promoter_status[event] = 1 # Turn the promoter ON
		elif event_index < events_indices[9]: # Promoter OFF event
			simulation_setup_and_state.last_event_type = 'promoter_OFF'
			event = event_index - events_indices[8]
			if model.promoter_status[event] == 0: # Promoter already OFF; cannot turn OFF again
				raise ValueError('Promoter OFF event selected for a promoter that is already OFF.')
			else:
				model.promoter_status[event] = 0 # Turn the promoter OFF
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