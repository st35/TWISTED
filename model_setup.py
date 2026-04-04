from utilities import *

from typing import Union
from math import floor

class GenomicSetup: # Class to hold genomic setup information
	def __init__(self, chromatin_type: str, gene_names: list[str], TSSes: list[float], gene_lengths: list[float], gene_directions: list[int], RNAP_on_rates: list[float], promoter_mode: str, buffer_length: float, **kwargs) -> None:
		self.chromatin_type = chromatin_type
		assert chromatin_type in ['prokaryotic', 'eukaryotic'], 'chromatin_type must be either "prokaryotic" or "eukaryotic".'

		self.gene_names = list(gene_names)
		self.TSSes = list(TSSes)
		self.gene_lengths = list(gene_lengths)
		self.gene_directions = list(gene_directions)
		self.RNAP_on_rates = list(RNAP_on_rates)
		assert all(direction in [1, -1] for direction in self.gene_directions), 'gene_directions must be either +1 (positive strand) or -1 (reverse strand).'
		assert len(self.gene_names) == len(self.TSSes) == len(self.gene_lengths) == len(self.gene_directions) == len(self.RNAP_on_rates), 'All gene-related lists must have the same length.'

		self.promoter_mode = promoter_mode
		assert promoter_mode in ['constitutive', 'non-constitutive'], 'promoter_mode must be either "constitutive" or "non-constitutive".'
		self.TF_on_off_rates = [(0.0, 0.0) for _ in self.gene_names]
		if self.promoter_mode == 'non-constitutive': # Non-constitutive promoters
			if 'TF_on_off_rates' not in kwargs:
				raise ValueError("For promoter_mode 'non-constitutive', 'TF_on_off_rates' argument must be provided.")
			self.TF_on_off_rates = [(rates[0], rates[1]) for rates in kwargs['TF_on_off_rates']]
		
		if self.chromatin_type == 'eukaryotic':
			if 'per_nucleosome_DNA_length' not in kwargs:
				self.per_nucleosome_DNA_length = 147.0*0.34 # Default: 147 bp
			else:
				self.per_nucleosome_DNA_length = float(kwargs['per_nucleosome_DNA_length'])*0.34
			if 'nucleosome_linker_length' not in kwargs:
				self.nucleosome_linker_length = 30.0*0.34 # Default: 30 bp
			else:
				self.nucleosome_linker_length = float(kwargs['nucleosome_linker_length'])*0.34
			if 'nucleosomes_are_steric_barriers_to_RNAPs' not in kwargs:
				self.nucleosomes_are_steric_barriers_to_RNAPs = True
			else:
				self.nucleosomes_are_steric_barriers_to_RNAPs = bool(kwargs['nucleosomes_are_steric_barriers_to_RNAPs'])
			if 'nucleosome_count' in kwargs:
				self.explicit_nucleosome_count = int(kwargs['nucleosome_count'])
			else:
				self.explicit_nucleosome_count = None
		
		self.clamp_left = 0.0 # Left end of DNA is at position 0 nm
		self.clamp_right = TSSes[0] + gene_lengths[0] + buffer_length if gene_directions[0] == 1 else TSSes[0] + buffer_length # Right end of DNA is at position beyond the last gene plus buffer length
	
	def get_total_nucleosome_count(self) -> int:
		if self.chromatin_type == 'prokaryotic':
			return 0
		if self.explicit_nucleosome_count is not None:
			return self.explicit_nucleosome_count
		nucl_count = 0
		start_pos = self.clamp_left + (self.nucleosome_linker_length / 2.0)
		while start_pos < self.clamp_right:
			if start_pos + self.per_nucleosome_DNA_length > self.clamp_right:
				break
			nucl_count += 1
			start_pos += self.per_nucleosome_DNA_length + self.nucleosome_linker_length
		return nucl_count
	
	def print_genomic_setup(self) -> None: # Utility function to print genomic setup information
		print('Chromatin type:', self.chromatin_type.capitalize())
		print('Promoter mode:', 'Constitutive' if self.promoter_mode == 'constitutive' else 'Non-constitutive')
		print('=' * 40)
		print('Name\tTSS (nm)\tLength (nm)\tDirection\tRNAP on-rate (1 / s)')
		print('-' * 40)
		for i in range(len(self.gene_names)):
			print(f'{self.gene_names[i]}\t{self.TSSes[i]}\t{self.gene_lengths[i]}\t{self.gene_directions[i]}\t{self.RNAP_on_rates[i]}')
			if self.promoter_mode == 'non-constitutive':
				print(f'\tTF on-rate: {self.TF_on_off_rates[i][0]}, TF off-rate: {self.TF_on_off_rates[i][1]}')
		print('=' * 40)

class ModelSetup: # Class to hold model setup parameters
	def __init__(self, w0: float = 1.85, chi: float = 0.05, eta: float = 0.0005, alpha: float = 1.5, v0: float = 20.0, tau_c: float = 12.0, force: float = 1.0, kBT: float = 4.1, TOP1_k0: float = 11.0, TOP1_theta: float = 0.25, TOP2_V0: float = 2.6, TOP2_k12: float = 2.0, RNAP_diameter: float = 15.0, TOPO_diameter: float = 15.0, generic_binding_protein_diameter: float = 15.0, clamps_status: tuple[str, str] = ('clamped', 'clamped'), finite_size_effect_flag: int = 1, supercoiling_relaxation_dynamics_mode: str = 'global_overall', mRNA_dynamics_mode: int = 0, model_observation_event_rate: float = 1.0 / 2.0, **kwargs) -> None:
		self.w0 = w0 # Default: 1.85 1 / nm
		self.h_dna = (2.0*3.14) / w0 # From w0*h_dna = 2*pi
		self.chi = chi # Default: 0.05 pN*nm*s
		self.eta = eta # Default: 0.0005 pN*nm^(1 - alpha)*s
		self.alpha = alpha # Default: 1.5
		self.v0 = v0 # Default: 20.0 nm / s
		self.tau_c = tau_c # Default: 12.0 pN*nm
		self.force = force # Default: 1.0 pN
		self.kBT = kBT # Default: 4.1 pN*nm (room temperature)
		self.TOP1_k0 = TOP1_k0 # Default: 11.0 1 / s
		self.TOP1_theta = TOP1_theta # Default: 0.25
		self.TOP2_V0 = TOP2_V0 # Default: 2.6 1 / s
		self.TOP2_k12 = TOP2_k12 # Default: 2.0
		self.RNAP_diameter = RNAP_diameter # Default: 15.0 nm; used for calculating steric hindrance effects involving RNAPs
		self.TOPO_diameter = TOPO_diameter # Default: 15.0 nm; used for calculating steric hindrance effects involving topoisomerases
		self.generic_binding_protein_diameter = generic_binding_protein_diameter # Default: 15.0 nm; used for calculating steric hindrance effects involving generic DNA-binding proteins
		assert len(clamps_status) == 2, 'clamps_status must be a tuple of two strings representing the status of the left and right clamps, respectively.'
		assert all(status in ['clamped', 'free'] for status in clamps_status), 'Each clamp status in clamps_status must be either "clamped" or "free".'
		self.left_clamp_status = 0 if clamps_status[0] == 'free' else 1
		self.right_clamp_status = 0 if clamps_status[1] == 'free' else 1

		self.finite_size_effect_flag = finite_size_effect_flag # Default: 1 (enabled); indiacates whether to consider finite size effect in torque calculations
		assert finite_size_effect_flag in [0, 1], 'finite_size_effect_flag must be either 0 (disabled) or 1 (enabled).'
		self.finite_size_effect_length = -1.0 # Length scale for finite size effect (in nm); only relevant if finite_size_effect_flag = 1
		if self.finite_size_effect_flag == 1:
			if 'finite_size_effect_length' not in kwargs:
				self.finite_size_effect_length = 1000.0*0.34 # Default: 1000 bp = 340 nm
			else:
				self.finite_size_effect_length = float(kwargs['finite_size_effect_length'])
		
		self.supercoiling_relaxation_dynamics_mode = supercoiling_relaxation_dynamics_mode
		assert supercoiling_relaxation_dynamics_mode in ['global_overall', 'global_per_segment', 'global_by_type', 'per_segment_by_type', 'topoisomerase_approximated', 'topoisomerase_based'], 'supercoiling_relaxation_dynamics_mode must be one of "global_overall", "global_per_segment", "global_by_type", "per_segment_by_type", "topoisomerase_approximated", or "topoisomerase_based".'
		self.global_supercoiling_relaxation_rate = 0.0 # Rate for global supercoiling relaxation (in 1 / s)
		if self.supercoiling_relaxation_dynamics_mode in ['global_overall', 'global_per_segment']: # Global supercoiling relaxation: supercoiling is relaxed at once throughout the genomic segment
			if 'global_supercoiling_relaxation_rate' not in kwargs:
				raise ValueError('For supercoiling_relaxation_dynamics_mode "global_overall" or "global_per_segment", "global_supercoiling_relaxation_rate" argument must be provided.')
			self.global_supercoiling_relaxation_rate = float(kwargs['global_supercoiling_relaxation_rate'])
		self.local_supercoiling_relaxation_rates = [0.0, 0.0] # Rates for local supercoiling relaxation of positive and negative supercoiling (in 1 / s)
		if self.supercoiling_relaxation_dynamics_mode in ['global_by_type', 'per_segment_by_type']: # Local supercoiling relaxation with specified rates; positive and negative supercoiling are relaxed independently at specified rates
			if 'local_supercoiling_relaxation_rates' not in kwargs:
				raise ValueError('For supercoiling_relaxation_dynamics_mode "global_by_type" or "per_segment_by_type", "local_supercoiling_relaxation_rates" argument must be provided.')
			if len(kwargs['local_supercoiling_relaxation_rates']) != 2:
				raise ValueError('"local_supercoiling_relaxation_rates" must be a list or tuple of two floats: [rate_positive, rate_negative].')
			self.local_supercoiling_relaxation_rates = [float(rate) for rate in kwargs['local_supercoiling_relaxation_rates']]
		self.TOP1_effective_relaxation_rate = 0.0 # Effective relaxation rate for TOP1 (in 1 / s)
		self.TOP2_effective_relaxation_rate = 0.0 # Effective relaxation rate for TOP2 (in 1 / s)
		if self.supercoiling_relaxation_dynamics_mode == 'topoisomerase_approximated': # Approximate topoisomerase-based supercoiling relaxation: approximate TOP1 / TOP2 activity by effective relaxation rates instead of explicitly modeling topoisomerase binding and unbinding dynamics
			if 'TOP1_effective_relaxation_rate' not in kwargs:
				raise ValueError('For supercoiling_relaxation_dynamics_mode "topoisomerase_approximated", "TOP1_effective_relaxation_rate" argument must be provided.')
			if 'TOP2_effective_relaxation_rate' not in kwargs:
				raise ValueError('For supercoiling_relaxation_dynamics_mode "topoisomerase_approximated", "TOP2_effective_relaxation_rate" argument must be provided.')
			self.TOP1_effective_relaxation_rate = float(kwargs['TOP1_effective_relaxation_rate'])
			self.TOP2_effective_relaxation_rate = float(kwargs['TOP2_effective_relaxation_rate'])
		self.topoisomerase_copy_numbers = [0, 0] # Copy numbers for TOP1 and TOP2
		self.topoisomerase_on_off_rates = [(0.0, 0.0), (0.0, 0.0)] # On and off rates for TOP1 and TOP2
		if self.supercoiling_relaxation_dynamics_mode == 'topoisomerase_based': # Supercoiling relaxation per segment based on topoisomerase binding and unbinding; requires specifying topoisomerase copy numbers and on/off rates
			if 'topoisomerase_copy_numbers' not in kwargs:
				raise ValueError('For supercoiling_relaxation_dynamics_mode "topoisomerase_based", "topoisomerase_copy_numbers" argument must be provided.')
			if len(kwargs['topoisomerase_copy_numbers']) != 2:
				raise ValueError('"topoisomerase_copy_numbers" must be a list or tuple of two integers: [num_TOP1, num_TOP2].')
			self.topoisomerase_copy_numbers = [int(num) for num in kwargs['topoisomerase_copy_numbers']]
			if 'topoisomerase_on_off_rates' not in kwargs:
				raise ValueError('For supercoiling_relaxation_dynamics_mode "topoisomerase_based", "topoisomerase_on_off_rates" argument must be provided.')
			if len(kwargs['topoisomerase_on_off_rates']) != 2:
				raise ValueError('"topoisomerase_on_off_rates" must be a list or tuple of two tuples: [(TOP1_on_rate, TOP1_off_rate), (TOP2_on_rate, TOP2_off_rate)].')
			self.topoisomerase_on_off_rates = [(float(rates[0]), float(rates[1])) for rates in kwargs['topoisomerase_on_off_rates']]

		self.mRNA_dynamics_mode = mRNA_dynamics_mode # 0: no mRNA degradation; 1: with mRNA degradation
		assert mRNA_dynamics_mode in [0, 1], 'mRNA_dynamics_mode must be either 0 (no mRNA degradation) or 1 (with mRNA degradation).'
		self.mRNA_degradation_rate = 0.0 # Rate for mRNA degradation (in 1 / s); only relevant if mRNA_dynamics_mode = 1
		if self.mRNA_dynamics_mode == 1:
			if 'mRNA_degradation_rate' not in kwargs:
				raise ValueError('For mRNA_dynamics_mode 1 (i.e., with mRNA degradation), "mRNA_degradation_rate" argument must be provided.')
			self.mRNA_degradation_rate = float(kwargs['mRNA_degradation_rate'])
		
		self.model_observation_event_rate = model_observation_event_rate # Rate for model observation events (in 1 / s)
		assert model_observation_event_rate > 0.0, 'model_observation_event_rate must be a positive float.'

		self.supercoiling_relaxation_dynamics_modes_with_no_steric_hindrance = ['global_overall', 'global_per_segment', 'global_by_type', 'per_segment_by_type', 'topoisomerase_approximated'] # List of supercoiling relaxation dynamics modes that do explicitly model topoisomerase binding and unbinding dynamics and therefore do not exert steric hindrance effects on RNAPs

class BindingProtein:
	def __init__(self, protein_name: str, total_copy_number: int, is_steric_barrier_to_RNAPs: bool, is_topological_barrier: bool, basal_on_rate: float, basal_off_rate: float, on_rate_func: callable = None, off_rate_func: callable = None, is_a_nucleosome: bool = False) -> None:
		self.protein_name = protein_name
		self.total_copy_number = total_copy_number
		self.is_steric_barrier_to_RNAPs = is_steric_barrier_to_RNAPs
		self.is_topological_barrier = is_topological_barrier
		self.basal_on_rate = basal_on_rate
		self.basal_off_rate = basal_off_rate
		if on_rate_func is None:
			self.on_rate_func = lambda segment_length, segment_sigma, *args: basal_on_rate*segment_length
		else:
			if not callable(on_rate_func):
				raise ValueError('on_rate_func must be a callable function if provided.')
			self.on_rate_func = lambda segment_length, segment_sigma, *args: on_rate_func(segment_length, segment_sigma, *args)*basal_on_rate*segment_length
		if off_rate_func is None:
			self.off_rate_func = lambda segment_length, segment_sigma, *args: basal_off_rate
		else:
			if not callable(off_rate_func):
				raise ValueError('off_rate_func must be a callable function if provided.')
			self.off_rate_func = lambda segment_length, segment_sigma, *args: off_rate_func(segment_length, segment_sigma, *args)*basal_off_rate
		self.is_a_nucleosome = is_a_nucleosome

class Model: # Class to hold the model, including genomic setup, model setup, and dynamic state variables
	def __init__(self, genomic_setup: GenomicSetup, model_setup: ModelSetup, binding_proteins: list[BindingProtein] = None) -> None:
		self.genomic_setup = genomic_setup # GenomicSetup object
		self.model_setup = model_setup # ModelSetup object
		self.x_dict = [[] for _ in genomic_setup.gene_names] # List of lists to hold positions of RNAPs for each gene
		self.theta_dict = [[] for _ in genomic_setup.gene_names] # List of lists to hold angular positions of RNAPs for each gene
		self.Lk = [(genomic_setup.clamp_right - genomic_setup.clamp_left) / model_setup.h_dna] # List to hold linking number of each DNA segment; initially, only one segment spanning the entire DNA
		self.promoter_status = [1 for _ in genomic_setup.gene_names] if genomic_setup.promoter_mode == 'constitutive' else [floor(uniform_random_in_interval(0.0, 2.0)) for _ in genomic_setup.gene_names] # List to hold promoter status (1: ON, 0: OFF) for each gene
		self.mRNA_counts = [0 for _ in genomic_setup.gene_names] # List to hold mRNA counts for each gene

		if model_setup.supercoiling_relaxation_dynamics_mode == 'topoisomerase_based': # State variables for topoisomerases
			self.topoisomerase_type = [0 for _ in range(model_setup.topoisomerase_copy_numbers[0])] + [1 for _ in range(model_setup.topoisomerase_copy_numbers[1])] # 0: TOP1, 1: TOP2
			self.topoisomerase_positions = [-1.0 for _ in range(model_setup.topoisomerase_copy_numbers[0] + model_setup.topoisomerase_copy_numbers[1])] # Positions of topoisomerases; -1.0 indicates unbound; initially all unbound
			self.topoisomerase_segment_indices = [-1 for _ in range(model_setup.topoisomerase_copy_numbers[0] + model_setup.topoisomerase_copy_numbers[1])] # Segment indices of bound topoisomerases; -1 indicates unbound; initially all unbound
			self.topoisomerase_status = [0 for _ in range(model_setup.topoisomerase_copy_numbers[0] + model_setup.topoisomerase_copy_numbers[1])] # Topoisomerase binding status; 0: unbound, 1: bound; initially all unbound
		
		if binding_proteins is None:
			binding_proteins = []
		self.binding_proteins = binding_proteins # List of BindingProtein objects representing other DNA-binding proteins in the system
		if genomic_setup.chromatin_type == 'eukaryotic':
			nucl_count = genomic_setup.get_total_nucleosome_count()
			nucleosomes = BindingProtein(protein_name = 'nucleosome', total_copy_number = nucl_count, is_steric_barrier_to_RNAPs = genomic_setup.nucleosomes_are_steric_barriers_to_RNAPs, is_topological_barrier = False, basal_on_rate = 1.2 / (genomic_setup.clamp_right - genomic_setup.clamp_left), basal_off_rate = 0.4, is_a_nucleosome = True)
			self.binding_proteins = [nucleosomes] + self.binding_proteins
		self.binding_proteins_positions = [[] for _ in self.binding_proteins] # List of lists to hold positions of each bound protein; each sublist corresponds to a binding protein type and contains the positions of all bound proteins of that type

class SimulationSetupAndState: # Class to hold simulation setup parameters
	def __init__(self, genomic_setup: GenomicSetup, simulation_end_mode: int, simulation_end_criterion: Union[float, list[int]], integration_time_resolution: float = 1.0e-1, RNAP_alive_status_check_interval: float = 1.0, max_RNAPs_to_recruit: list[int] = None) -> None:
		self.simulation_end_mode = simulation_end_mode # 0: time-based, 1: event-based
		assert simulation_end_mode in [0, 1], 'simulation_end_mode must be either 0 (time-based) or 1 (event-based).'

		self.simulation_end_criterion = simulation_end_criterion # Float for time-based (in s) or list of integers for event-based (number of RNAPs that must finish transcription for each gene)
		if simulation_end_mode == 0:
			self.simulation_end_time = float(simulation_end_criterion)
		else:
			if not isinstance(simulation_end_criterion, list):
				raise ValueError('For simulation_end_mode 1 (i.e., event-based), simulation_end_criterion must be a list of integers representing event counts for each gene.')
			if len(simulation_end_criterion) != len(genomic_setup.gene_names):
				raise ValueError('Length of simulation_end_criterion list must match the number of genes in genomic_setup.')
			self.simulation_end_event_counts = list(simulation_end_criterion)

		self.integration_time_resolution = integration_time_resolution # Time resolution for integration (in s)
		assert integration_time_resolution > 0.0, 'integration_time_resolution must be a positive float.'

		self.RNAP_alive_status_check_interval = RNAP_alive_status_check_interval # Interval for checking RNAP alive status (in s)
		assert RNAP_alive_status_check_interval > 0.0, 'RNAP_alive_status_check_interval must be a positive float.'

		if max_RNAPs_to_recruit is not None:
			assert len(max_RNAPs_to_recruit) == len(genomic_setup.gene_names), 'Length of max_RNAPs_to_recruit list must match the number of genes in genomic_setup.'
		self.max_RNAPs_to_recruit = max_RNAPs_to_recruit # Maximum number of RNAPs to recruit for each gene
		if self.max_RNAPs_to_recruit is not None and self.simulation_end_mode == 1:
			for i in range(len(self.simulation_end_event_counts)):
				if self.simulation_end_event_counts[i] > self.max_RNAPs_to_recruit[i]:
					raise ValueError(f'For gene index {i}, simulation_end_event_counts ({self.simulation_end_event_counts[i]}) cannot be greater than max_RNAPs_to_recruit ({self.max_RNAPs_to_recruit[i]}).')

		self.RNAPs_finished_transcription = [0 for _ in genomic_setup.gene_names] # List to hold counts of RNAPs that have finished transcription for each gene
		self.RNAPs_exit_positions = [[] for _ in genomic_setup.gene_names] # List of lists to hold exit positions of RNAPs for each gene
		self.RNAP_recruitment_times = [[] for _ in genomic_setup.gene_names] # List of lists to hold recruitment times of RNAPs for each gene
		self.RNAP_exit_times = [[] for _ in genomic_setup.gene_names] # List of lists to hold exit times of RNAPs for each gene

		self.curr_simulation_time = 0.0 # Current simulation time (in s)
		self.simulation_completed = False # Flag indicating whether the simulation has completed
	
	def calculate_RNAP_transcription_rates(self, model: Model) -> list[list[float]]: # Calculate and return the transcription rates (in bp / s) for each RNAP that has finished transcription for each gene
		transcription_rates = [[] for _ in model.genomic_setup.gene_names]
		for gene_index in range(len(model.genomic_setup.gene_names)):
			for i in range(len(self.RNAP_exit_times[gene_index])):
				recruitment_time = self.RNAP_recruitment_times[gene_index][i]
				exit_time = self.RNAP_exit_times[gene_index][i]
				time_interval = exit_time - recruitment_time
				assert time_interval > 0.0, 'Exit time must be greater than recruitment time for transcription rate calculation.'

				distance_covered = self.RNAPs_exit_positions[gene_index][i] - model.genomic_setup.TSSes[gene_index] if model.genomic_setup.gene_directions[gene_index] == 1 else model.genomic_setup.TSSes[gene_index] - self.RNAPs_exit_positions[gene_index][i]
				assert distance_covered > 0.0 and distance_covered >= model.genomic_setup.gene_lengths[gene_index], 'Distance covered must be positive and at least equal to gene length for transcription rate calculation.'

				distance_covered = distance_covered / 0.34 # Convert nm to bp
				transcription_rates[gene_index].append(distance_covered / time_interval)
		
		return transcription_rates