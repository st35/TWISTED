from model_setup import *
from model_dynamics import *
from utilities import *

from random import random
from math import log

def simulate_dynamics(model: Model, simulation_setup_and_state: SimulationSetupAndState) -> None: # Main function to simulate the dynamics of the model until the end condition is met
	while True:
		RNAP_gene_index, state_vector = get_state_vectors_from_dicts(model) # Get current state vectors from model dictionaries

		p0 = random() # Generate a random number for event time calculation
		dt, _ = integrate(model, simulation_setup_and_state, simulation_setup_and_state.curr_simulation_time, state_vector, RNAP_gene_index, p0) # Integrate dynamics until next event
		assert dt is not None, 'Integration failed to return a valid time step.'
		update_dicts_from_state_vector(model, RNAP_gene_index, state_vector) # Update model dictionaries from updated state vectors
		simulation_setup_and_state.curr_simulation_time += dt # Update current simulation time
		segments_lengths, segments_sigmas, segments_torques, segments_dna_states, segments_writhe_fractions = calculate_segments_attributes(model, RNAP_gene_index, state_vector) # Recalculate segment attributes after integration

		rates_vector, events_indices = get_events_rates(model, RNAP_gene_index, state_vector) # Get event rates and event indices
		p1 = random() # Generate a random number for event selection
		event_index = select_event_based_on_propensities(rates_vector, p1) # Select event based on propensities

		if event_index < events_indices[0]: # RNAP recruitment event
			event = event_index - 0
			if get_TSS_steric_hindrance_status(model.genomic_setup.TSSes[event], RNAP_gene_index, state_vector, model.model_setup.between_RNAPs_steric_effect_cutoff) == 1: # Steric hindrance at TSS; recruitment fails
				pass
			else:
				event_gene_index = event
				model.x_dict[event_gene_index].append(model.genomic_setup.TSSes[event_gene_index])
				model.theta_dict[event_gene_index].append(0.0)
				update_Lk_vector_after_RNAP_recruitment(model, event_gene_index, RNAP_gene_index, state_vector, segments_lengths, segments_sigmas) # Update linking number vector after RNAP recruitment
				simulation_setup_and_state.RNAP_recruitment_times[event_gene_index].append(simulation_setup_and_state.curr_simulation_time)
		elif event_index < events_indices[1]: # Model observation event: a dummy event to ensure simulation time progresses even if all biological events rates are zero
			pass
		elif event_index < events_indices[2]: # Global supercoiling relaxation event; supercoiling_relaxation_dynamics_mode = 0
			segments_Lk0 = [segments_length / model.model_setup.h_dna for segments_length in segments_lengths]
			model.Lk = [Lk0 for Lk0 in segments_Lk0]
		elif event_index < events_indices[3]: # Local supercoiling relaxation event: supercoiling_relaxation_dynamics_mode = 1
			event = event_index - events_indices[2]
			segments_Lk0 = [segments_length / model.model_setup.h_dna for segments_length in segments_lengths]
			if event == 0: # Relax positive supercoiling only
				model.Lk = [model.Lk0[i] if segments_sigmas[i] > 0.0 else model.Lk[i] for i in range(len(segments_lengths))]
			elif event == 1: # Relax negative supercoiling only
				model.Lk = [model.Lk0[i] if segments_sigmas[i] < 0.0 else model.Lk[i] for i in range(len(segments_lengths))]
			else:
				raise ValueError('Invalid event index for local supercoiling relaxation.')
		elif event_index < events_indices[4]: # Topoisomerase-mediated supercoiling relaxation: supercoiling_relaxation_dynamics_mode = 2
			raise NotImplementedError('Topoisomerase-mediated supercoiling relaxation not implemented yet.')
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
	
	transcription_rates = simulation_setup_and_state.calculate_RNAP_transcription_times(model) # Calculate average RNAP velocities for each gene