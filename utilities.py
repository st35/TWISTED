from model_setup import *

from random import random

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
		explicit_RNAP_on_rates = kwargs['explicit_RNAP_on_rates']
		if len(explicit_RNAP_on_rates) != len(gene_names):
			raise ValueError('Length of explicit_RNAP_on_rates must match number of genes.')
		RNAP_on_rates = [RNAP_on_rates[i]*explicit_RNAP_on_rates[i] for i in range(len(gene_names))]

	if promoter_mode == 'non-constitutive':
		if 'TF_on_off_rates' not in kwargs:
			raise ValueError('For promoter_mode "non-constitutive", "TF_on_off_rates" argument must be provided.')
		TF_on_off_rates = kwargs['TF_on_off_rates']
		if len(TF_on_off_rates) != len(gene_names):
			raise ValueError('Length of TF_on_off_rates must match number of genes.')
		
		return GenomicSetup(chromatin_type, gene_names, TSSes, gene_lengths, gene_directions, RNAP_on_rates, promoter_mode, buffer_length, TF_on_off_rates = TF_on_off_rates)
	
	return GenomicSetup(chromatin_type, gene_names, TSSes, gene_lengths, gene_directions, RNAP_on_rates, promoter_mode, buffer_length)
	
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

def get_TSS_steric_hindrance_status(TSS_position: float, RNAP_gene_index: list[int], state_vector: list[float], between_RNAPs_steric_effect_cutoff: float) -> int: # Check if there is steric hindrance at the TSS position due to existing RNAPs; return 1 if hindered, 0 if not
	RNAP_count = len(RNAP_gene_index)
	for x in state_vector[:RNAP_count]:
		if abs(x - TSS_position) < between_RNAPs_steric_effect_cutoff:
			return 1
	return 0

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