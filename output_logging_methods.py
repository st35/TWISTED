from typing import TextIO

from utilities import *
from model_setup import *
from simulate_dynamics import *

def print_at_each_integration_step(model: Model, simulation_setup_and_state: SimulationSetupAndState, t: float, state_vector: list[float], RNAPs_pos_file: TextIO = None, RNAPs_dx_dt_file: TextIO = None, RNAPs_dtheta_dt_file: TextIO = None, sigma_file: TextIO = None, torque_file: TextIO = None, dLk_dt_file: TextIO = None, nucleosome_position_file: TextIO = None, promoter_status_file: TextIO = None, mRNA_count_file: TextIO = None, barrier_position_file: TextIO = None, binding_proteins_file: TextIO = None, screen_logging = False) -> None:
	RNAP_gene_index, state_vector = get_state_vectors_from_dicts(model)
	RNAP_count = len(RNAP_gene_index)
	dstate_dt = model_dynamics(t, state_vector, RNAP_gene_index, model, simulation_setup_and_state)
	_, segments_sigmas, segments_torques, _, _, _ = calculate_segments_attributes(model, RNAP_gene_index, state_vector)

	if RNAPs_pos_file is not None:
		RNAPs_pos_file.write(str(t) + '\t' + '\t'.join([str(state_vector[i]) for i in range(RNAP_count)]) + '\n')
	if RNAPs_dx_dt_file is not None:
		RNAPs_dx_dt_file.write(str(t) + '\t' + '\t'.join([str(dstate_dt[i]) for i in range(RNAP_count)]) + '\n')
	if RNAPs_dtheta_dt_file is not None:
		RNAPs_dtheta_dt_file.write(str(t) + '\t' + '\t'.join([str(dstate_dt[i + RNAP_count]) for i in range(RNAP_count)]) + '\n')
	if sigma_file is not None:
		sigma_file.write(str(t) + '\t' + '\t'.join([str(segments_sigmas[i]) for i in range(len(segments_sigmas))]) + '\n')
	if torque_file is not None:
		torque_file.write(str(t) + '\t' + '\t'.join([str(segments_torques[i]) for i in range(len(segments_sigmas))]) + '\n')
	if dLk_dt_file is not None:
		dLk_dt_file.write(str(t) + '\t' + '\t'.join([str(state_vector[2 * RNAP_count + i]) for i in range(len(segments_sigmas))]) + '\n')
	if nucleosome_position_file is not None:
		if model.genomic_setup.chromatin_type == 'eukaryotic':
			nucleosome_position_file.write(str(t) + '\t' + '\t'.join([str(model.binding_proteins_positions[0][i]) for i in range(len(model.binding_proteins_positions[0]))]) + '\n')
	if barrier_position_file is not None and model.genomic_setup.are_multiple_chromosomes_present:
		if len(model.binding_proteins) > 1:
			barrier_position_file.write(str(t) + '\t' + '\t'.join([str(model.binding_proteins_positions[-1][i]) for i in range(len(model.binding_proteins_positions[-1]))]) + '\n')
	if binding_proteins_file is not None:
		if model.genomic_setup.chromatin_type == 'prokaryotic':
			if model.genomic_setup.are_multiple_chromosomes_present:
				binding_proteins_file.write(str(t) + '\t')
				for protein_index in range(len(model.binding_proteins) - 1):
					binding_proteins_file.write('\t'.join([str(model.binding_proteins_positions[protein_index][i]) for i in range(len(model.binding_proteins_positions[protein_index]))]) + '\t')
				binding_proteins_file.write('\n')
			else:
				binding_proteins_file.write(str(t) + '\t')
				for protein_index in range(len(model.binding_proteins)):
					binding_proteins_file.write('\t'.join([str(model.binding_proteins_positions[protein_index][i]) for i in range(len(model.binding_proteins_positions[protein_index]))]) + '\t')
				binding_proteins_file.write('\n')
		else:
			if model.genomic_setup.are_multiple_chromosomes_present:
				binding_proteins_file.write(str(t) + '\t')
				for protein_index in range(1, len(model.binding_proteins) - 1):
					binding_proteins_file.write('\t'.join([str(model.binding_proteins_positions[protein_index][i]) for i in range(len(model.binding_proteins_positions[protein_index]))]) + '\t')
				binding_proteins_file.write('\n')
			else:
				binding_proteins_file.write(str(t) + '\t')
				for protein_index in range(1, len(model.binding_proteins)):
					binding_proteins_file.write('\t'.join([str(model.binding_proteins_positions[protein_index][i]) for i in range(len(model.binding_proteins_positions[protein_index]))]) + '\t')
				binding_proteins_file.write('\n')
	if promoter_status_file is not None:
		promoter_status_file.write(str(t) + '\t' + '\t'.join([str(model.promoter_status[i]) for i in range(len(model.genomic_setup.gene_names))]) + '\n')
	if mRNA_count_file is not None:
		currently_transcribing_RNAP_counts = []
		for gene_index in range(len(model.genomic_setup.gene_names)):
			count = 0
			for RNAP_index in range(RNAP_count):
				if RNAP_gene_index[RNAP_index] == gene_index:
					count += 1
			currently_transcribing_RNAP_counts.append(count)
		mRNA_count_file.write(str(t) + '\t' + '\t'.join([str(count) for count in currently_transcribing_RNAP_counts]) + '\t' + '\t'.join([str(simulation_setup_and_state.RNAPs_finished_transcription[i]) for i in range(len(model.genomic_setup.gene_names))]) + '\t' + '\t'.join([str(model.mRNA_counts[i]) for i in range(len(model.genomic_setup.gene_names))]) + '\n')

	if screen_logging:
		print(str(t) + '\t' + str(RNAP_count) + '\t' + '\t'.join([str(simulation_setup_and_state.RNAPs_finished_transcription[i]) for i in range(len(model.genomic_setup.gene_names))]) + '\t' + '\t'.join([str(model.mRNA_counts[i]) for i in range(len(model.genomic_setup.gene_names))]) + '\t' + '\t'.join([str(status) for status in model.promoter_status]))

def print_at_each_Gillespie_step(model: Model, simulation_setup_and_state: SimulationSetupAndState, events_log_file: TextIO = None, log_nucleosomal_events = False) -> None:
	if events_log_file is not None:
		last_event_index = str(simulation_setup_and_state.last_event_index) if simulation_setup_and_state.last_event_index >= 0 else 'NA'
		last_event_type = simulation_setup_and_state.last_event_type if simulation_setup_and_state.last_event_type is not None else 'NA'
		if 'nucleosome' in last_event_type and not log_nucleosomal_events:
			return
		events_log_file.write(str(simulation_setup_and_state.curr_simulation_time) + '\t' + last_event_index + '\t' + last_event_type + '\n')

def print_at_end_of_simulation(model: Model, simulation_setup_and_state: SimulationSetupAndState, transcription_rates_file: TextIO = None, final_mRNA_counts_file: TextIO = None) -> None:
	transcription_rates = simulation_setup_and_state.calculate_RNAP_transcription_rates(model)
	if transcription_rates_file is not None:
		for gene_index, gene_name in enumerate(model.genomic_setup.gene_names):
			transcription_rates_file.write(gene_name + '\t' + '\t'.join([str(rate) for rate in transcription_rates[gene_index]]) + '\n')
	if final_mRNA_counts_file is not None:
		final_mRNA_counts_file.write('\t'.join([str(count) for count in model.mRNA_counts]) + '\n')