from model_setup import *

from random import random
from bisect import bisect
from math import tanh

def read_genes_information(filename: str) -> tuple[list[str], list[float], list[float], list[int], list[float]]: # Read gene information from a tab-delimited file
	gene_names = []
	TSSes = []
	gene_lengths = []
	gene_directions = []
	RNAP_on_rates = []

	with open(filename, 'r') as f:
		for line in f:
			l = line.strip().split('\t')
			gene_names.append(l[0].strip())
			TSSes.append(float(l[1].strip())*0.34) # Input location in bp; convert to nm
			gene_lengths.append(float(l[2].strip())*0.34) # Input length in bp; convert to nm
			gene_directions.append(int(l[3].strip())) # +1 for + strand, -1 for - strand
			RNAP_on_rates.append(float(l[4].strip())) # in 1 / s
	
	return gene_names, TSSes, gene_lengths, gene_directions, RNAP_on_rates

def construct_genomic_setup(filename: str, chromatin_type: str, promoter_mode: str = 'constitutive', buffer_length: float = 10000.0*0.34, **kwargs) -> GenomicSetup: # Return a GenomicSetup object constructed from a gene information file; buffer_length in nm
	gene_names, TSSes, gene_lengths, gene_directions, RNAP_on_rates = read_genes_information(filename)

	if 'explicit_RNAP_on_rates' in kwargs: # Modify RNAP_on_rates if explicit rates are provided
		explicit_RNAP_on_rates = kwargs.pop('explicit_RNAP_on_rates')
		if len(explicit_RNAP_on_rates) != len(gene_names):
			raise ValueError('Length of explicit_RNAP_on_rates must match number of genes.')
		RNAP_on_rates = [RNAP_on_rates[i]*explicit_RNAP_on_rates[i] for i in range(len(gene_names))]

	if promoter_mode == 'non-constitutive':
		if 'TF_on_off_rates' not in kwargs:
			raise ValueError('For promoter_mode "non-constitutive", "TF_on_off_rates" argument must be provided.')
		TF_on_off_rates = kwargs['TF_on_off_rates']
		if len(TF_on_off_rates) != len(gene_names):
			raise ValueError('Length of TF_on_off_rates must match number of genes.')

	return GenomicSetup(chromatin_type, gene_names, TSSes, gene_lengths, gene_directions, RNAP_on_rates, promoter_mode, buffer_length, **kwargs)
	
def uniform_random_in_interval(start: float, end: float) -> float: # Generate a uniform random number in [start, end)
	return start + (end - start)*random()

def get_spot_segment_index(spot: float, segments_lengths: list[float]) -> int: # Get the index of the segment containing the given spot
	spot_segment_index = -1
	segment_index = len(segments_lengths) - 1
	cumulative_length = 0.0

	while segment_index >= 0:
		cumulative_length += segments_lengths[segment_index]
		if spot < cumulative_length:
			spot_segment_index = segment_index
			break
		segment_index -= 1
	
	return spot_segment_index

def get_nucleosome_occupied_fraction_per_segment(model: Model, segments_lengths: list[float], segment_index: int) -> float:
	nucl_positions = []
	for i in range(len(model.binding_proteins)):
		if model.binding_proteins[i].is_a_nucleosome:
			nucl_positions = nucl_positions + model.binding_proteins_positions[i]
	segment_left_end = model.genomic_setup.clamp_left + sum(segments_lengths[segment_index + 1:])
	segment_right_end = segment_left_end + segments_lengths[segment_index]
	assert segment_left_end < segment_right_end, 'Invalid segment boundaries.'

	occupied_length = 0.0
	for nucl_pos in nucl_positions:
		nucl_left_edge = nucl_pos - (model.genomic_setup.per_nucleosome_DNA_length + model.genomic_setup.nucleosome_linker_length) / 2.0
		nucl_right_edge = nucl_pos + (model.genomic_setup.per_nucleosome_DNA_length + model.genomic_setup.nucleosome_linker_length) / 2.0
		assert nucl_left_edge < nucl_right_edge, 'Invalid nucleosome boundaries.'

		overlap_left = max(segment_left_end, nucl_left_edge)
		overlap_right = min(segment_right_end, nucl_right_edge)
		if overlap_left < overlap_right:
			occupied_length += (overlap_right - overlap_left)

	return occupied_length / segments_lengths[segment_index]

def get_TSS_steric_hindrance_status(model: Model, TSS_position: float, RNAP_gene_index: list[int], state_vector: list[float]) -> tuple[int, float, int]: # Check for steric hindrance at the TSS position; return a tuple of (steric_hindrance_status, blocking_entity_position, blocking_entity_id); steric_hindrance_status: 0 for no hindrance, 1 for hindrance; blocking_entity_id: -1 for RNAP, 0 to len(model.binding_proteins) - 1 for bound protein index; blocking_entity_position: position of the entity causing steric hindrance if present, None if no hindrance
	RNAP_count = len(RNAP_gene_index)
	for x in state_vector[:RNAP_count]:
		if abs(x - TSS_position) < model.model_setup.RNAP_diameter:
			return 1, x, -1
		
	nucl_positions = []
	nucl_ids = []
	for i in range(len(model.binding_proteins)):
		if model.binding_proteins[i].is_a_nucleosome and model.binding_proteins[i].is_steric_barrier_to_RNAPs:
			nucl_positions = nucl_positions + model.binding_proteins_positions[i]
			nucl_ids = nucl_ids + [i for _ in model.binding_proteins_positions[i]]
	for nucl_pos, nucl_id in zip(nucl_positions, nucl_ids):
		if abs(nucl_pos - TSS_position) < (model.model_setup.RNAP_diameter + model.genomic_setup.per_nucleosome_DNA_length + model.genomic_setup.nucleosome_linker_length) / 2.0:
			return 1, nucl_pos, nucl_id
	
	bound_protein_positions = []
	bound_protein_ids = []
	for i in range(len(model.binding_proteins)):
		if model.binding_proteins[i].is_a_nucleosome is False and model.binding_proteins[i].is_steric_barrier_to_RNAPs:
			bound_protein_positions = bound_protein_positions + model.binding_proteins_positions[i]
			bound_protein_ids = bound_protein_ids + [i for _ in model.binding_proteins_positions[i]]
	for protein_pos, protein_id in zip(bound_protein_positions, bound_protein_ids):
		if abs(protein_pos - TSS_position) < (model.model_setup.RNAP_diameter + model.model_setup.generic_binding_protein_diameter) / 2.0:
			return 1, protein_pos, protein_id
	
	return 0, None, None

def is_protein_binding_blocked(model: Model, RNAP_gene_index: list[int], state_vector: list[float], protein_index: int, binding_position: float) -> int:
	RNAP_count = len(RNAP_gene_index)
	for x in state_vector[:RNAP_count]:
		if model.binding_proteins[protein_index].is_a_nucleosome:
			if abs(x - binding_position) < (model.model_setup.RNAP_diameter + model.genomic_setup.per_nucleosome_DNA_length + model.genomic_setup.nucleosome_linker_length) / 2.0:
				return 1
		else:
			if abs(x - binding_position) < (model.model_setup.RNAP_diameter + model.model_setup.generic_binding_protein_diameter) / 2.0:
				return 1
	
	if model.binding_proteins[protein_index].is_a_nucleosome:
		bound_nucl_positions = []
		for i in range(len(model.binding_proteins)):
			if model.binding_proteins[i].is_a_nucleosome:
				bound_nucl_positions = bound_nucl_positions + model.binding_proteins_positions[i]
		for nucl_pos in bound_nucl_positions:
			if abs(nucl_pos - binding_position) < (model.genomic_setup.per_nucleosome_DNA_length + model.genomic_setup.nucleosome_linker_length):
				return 1
	
	if model.binding_proteins[protein_index].is_a_nucleosome is False:
		bound_protein_positions = []
		for i in range(len(model.binding_proteins)):
			bound_protein_positions = bound_protein_positions + model.binding_proteins_positions[i]
		for protein_pos in bound_protein_positions:
			if abs(protein_pos - binding_position) < model.model_setup.generic_binding_protein_diameter:
				return 1

	return 0

def get_ordering_of_RNAPs_and_proteins(model: Model, RNAP_gene_index: list[int], state_vector: list[float]) -> tuple[list[float], list[str]]:
	RNAP_positions = list(state_vector[:len(RNAP_gene_index)])
	RNAP_ids = ['RNAP_' + str(i) for i in range(len(RNAP_gene_index))]

	bound_protein_positions = []
	bound_protein_ids = []
	for i in range(len(model.binding_proteins)):
		if model.binding_proteins[i].is_steric_barrier_to_RNAPs:
			bound_protein_positions = bound_protein_positions + model.binding_proteins_positions[i]
			bound_protein_ids = bound_protein_ids + [str(i) for _ in range(len(model.binding_proteins_positions[i]))]
	
	all_positions = RNAP_positions + bound_protein_positions
	all_ids = RNAP_ids + bound_protein_ids

	zipped_pairs = sorted(zip(all_positions, all_ids), key = lambda pair: pair[0], reverse = False)
	sorted_positions, sorted_ids = zip(*zipped_pairs)

	return list(sorted_positions), list(sorted_ids)

def get_nearest_steric_obstacles_for_RNAP(RNAP_index: int, RNAP_pos:float, sorted_positions: list[float], sorted_ids: list[str]) -> tuple[Union[float, None], Union[str, None], Union[float, None], Union[str, None]]:
	RNAP_id = 'RNAP_' + str(RNAP_index)
	index_in_ordering = bisect(sorted_positions, RNAP_pos)
	index_in_ordering = index_in_ordering - 1 # bisect returns the insertion point, so subtract 1 to get the index of the RNAP itself
	assert index_in_ordering >= 0 and index_in_ordering < len(sorted_ids) and sorted_ids[index_in_ordering] == RNAP_id, 'RNAP ID not found in sorted IDs list.'

	left_obstacle_index = index_in_ordering - 1 if index_in_ordering - 1 >= 0 else None
	right_obstacle_index = index_in_ordering + 1 if index_in_ordering + 1 < len(sorted_positions) else None

	left_obstacle_position = sorted_positions[left_obstacle_index] if left_obstacle_index is not None else None
	left_obstacle_id = sorted_ids[left_obstacle_index] if left_obstacle_index is not None else None

	right_obstacle_position = sorted_positions[right_obstacle_index] if right_obstacle_index is not None else None
	right_obstacle_id = sorted_ids[right_obstacle_index] if right_obstacle_index is not None else None

	return left_obstacle_position, left_obstacle_id, right_obstacle_position, right_obstacle_id

def check_separation_between_nucleosomes(model: Model) -> None:
	nucl_positions = []
	for i in range(len(model.binding_proteins)):
		if model.binding_proteins[i].is_a_nucleosome:
			nucl_positions = nucl_positions + model.binding_proteins_positions[i]
	nucl_positions.sort()
	for i in range(len(nucl_positions) - 1):
		if abs(nucl_positions[i + 1] - nucl_positions[i]) < model.genomic_setup.per_nucleosome_DNA_length + model.genomic_setup.nucleosome_linker_length:
			raise ValueError('Nucleosomes are too close to each other.')

def check_separation_between_nucleosomes_and_RNAPs(model: Model, RNAP_gene_index: list[int], state_vector: list[float]) -> None:
	RNAP_count = len(RNAP_gene_index)
	nucl_positions = []
	for i in range(len(model.binding_proteins)):
		if model.binding_proteins[i].is_a_nucleosome:
			nucl_positions = nucl_positions + model.binding_proteins_positions[i]
	for x in state_vector[:RNAP_count]:
		for nucl_pos in nucl_positions:
			if abs(x - nucl_pos) < (model.model_setup.RNAP_diameter + model.genomic_setup.per_nucleosome_DNA_length) / 2.0:
				raise ValueError('RNAP and nucleosome are too close to each other.')

def get_steric_hindrance_factor(model:Model, separation: float, steric_hindrance_distance: float) -> float:
	return (1.0 / 2.0)*(1 + tanh((separation - steric_hindrance_distance) / model.model_setup.steric_hindrance_constraint_parameter))

def find_and_remove_from_list(lst: list[float], value: float, tolerance: float = 1e-6) -> None: # Find the first occurrence of a value in a list within a given tolerance and remove it; raise an error if not found
	for i in range(len(lst)):
		if abs(lst[i] - value) < tolerance:
			lst.pop(i)
			return

	raise ValueError('Value not found in list within tolerance.')

def select_event_based_on_propensities(rates_vector: list[float], p: float) -> Union[int, None]: # Select an event index based on the given rates vector and random number p in [0, 1)
	a0 = sum(rates_vector)
	if a0 > 0.0:
		event_index = -1

		sum_prev = 0.0
		sum_new = 0.0
		for i in range(len(rates_vector)):
			sum_new = sum_prev + (rates_vector[i] / a0)
			if p >= sum_prev and p < sum_new:
				event_index = i
				break
			sum_prev = sum_new
		
		return event_index
	return None

def print_list(name: str, lst: list[float]) -> None: # Utility function to print a list with a name
	print(f'{name}:', end = ' ')
	for val in lst:
		print(f'{val}', end = ' ')
	print()