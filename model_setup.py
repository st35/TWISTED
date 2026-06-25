from utilities import *

from typing import Union
from math import floor
import warnings
import random
import dill

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
		if promoter_mode == 'non-constitutive':
			if 'TF_on_off_rates' not in kwargs:
				raise ValueError('For promoter_mode "non-constitutive", "TF_on_off_rates" argument must be provided.')
			if len(kwargs['TF_on_off_rates']) != len(self.gene_names):
				raise ValueError('Length of TF_on_off_rates must match number of genes.')
			for rates in kwargs['TF_on_off_rates']:
				if len(rates) != 2:
					raise ValueError('Each element in TF_on_off_rates must be a tuple or list of two floats: (TF_on_rate, TF_off_rate).')
			self.TF_on_off_rates = [tuple(rates) for rates in kwargs['TF_on_off_rates']]
		else:
			self.TF_on_off_rates = [(0.0, 0.0) for _ in self.gene_names]
		
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
			if 'nucleosome_on_rate_func' not in kwargs:
				self.nucleosome_on_rate_func = None
			else:
				if not callable(kwargs['nucleosome_on_rate_func']):
					raise ValueError('nucleosome_on_rate_func must be a callable function if provided.')
				self.nucleosome_on_rate_func = kwargs['nucleosome_on_rate_func']
			if 'nucleosome_off_rate_func' not in kwargs:
				self.nucleosome_off_rate_func = None
			else:
				if not callable(kwargs['nucleosome_off_rate_func']):
					raise ValueError('nucleosome_off_rate_func must be a callable function if provided.')
				self.nucleosome_off_rate_func = kwargs['nucleosome_off_rate_func']
			if 'nucleosomes_can_be_displaced_at_TSS_by_RNAP' not in kwargs:
				self.nucleosomes_can_be_displaced_at_TSS_by_RNAP = False
			else:
				self.nucleosomes_can_be_displaced_at_TSS_by_RNAP = bool(kwargs['nucleosomes_can_be_displaced_at_TSS_by_RNAP'])
		
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
	def __init__(self, w0: float = 1.85, chi: float = 0.05, eta: float = 0.0005, alpha: float = 1.5, v0: float = 20.0, tau_c: float = 12.0, force: float = 1.0, kBT: float = 4.1, TOP1_k0: float = 11.0, TOP1_theta: float = 0.25, TOP2_V0: float = 2.6, TOP2_k12: float = 2.0, RNAP_diameter: float = 15.0, generic_binding_protein_diameter: float = 15.0, steric_hindrance_constraint_parameter: float = 2.0, clamps_status: tuple[str, str] = ('clamped', 'clamped'), finite_size_effect_flag: int = 1, supercoiling_relaxation_dynamics_mode: str = 'global_overall', mRNA_dynamics_mode: int = 0, model_observation_event_rate: float = 1.0 / 2.0, **kwargs) -> None:
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
		self.generic_binding_protein_diameter = generic_binding_protein_diameter # Default: 15.0 nm; used for calculating steric hindrance effects involving generic DNA-binding proteins
		self.steric_hindrance_constraint_parameter = steric_hindrance_constraint_parameter # Default: 2.0; used for calculating steric hindrance effects
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
		if self.supercoiling_relaxation_dynamics_mode == 'topoisomerase_based': # Supercoiling relaxation per segment based on topoisomerase binding and unbinding; requires specifying topoisomerase copy numbers and on/off rates
			raise NotImplementedError('supercoiling_relaxation_dynamics_mode "topoisomerase_based" is not yet implemented.')

		self.mRNA_dynamics_mode = mRNA_dynamics_mode # 0: no mRNA degradation; 1: with mRNA degradation
		assert mRNA_dynamics_mode in [0, 1], 'mRNA_dynamics_mode must be either 0 (no mRNA degradation) or 1 (with mRNA degradation).'
		self.mRNA_degradation_rate = 0.0 # Rate for mRNA degradation (in 1 / s); only relevant if mRNA_dynamics_mode = 1
		if self.mRNA_dynamics_mode == 1:
			if 'mRNA_degradation_rate' not in kwargs:
				raise ValueError('For mRNA_dynamics_mode 1 (i.e., with mRNA degradation), "mRNA_degradation_rate" argument must be provided.')
			self.mRNA_degradation_rate = float(kwargs['mRNA_degradation_rate'])
		
		self.model_observation_event_rate = model_observation_event_rate # Rate for model observation events (in 1 / s)
		assert model_observation_event_rate > 0.0, 'model_observation_event_rate must be a positive float.'

class BindingProtein:
	def __init__(self, protein_name: str, total_copy_number: int, is_steric_barrier_to_RNAPs: bool, is_topological_barrier: bool, basal_on_rate: float, basal_off_rate: float, on_rate_func: callable = None, off_rate_func: callable = None, is_a_nucleosome: bool = False, can_be_displaced_at_TSS_by_RNAP: bool = False) -> None:
		self.protein_name = protein_name
		self.total_copy_number = total_copy_number

		self.is_steric_barrier_to_RNAPs = is_steric_barrier_to_RNAPs
		self.is_topological_barrier = is_topological_barrier
		if self.is_topological_barrier and not self.is_steric_barrier_to_RNAPs:
			raise ValueError('A protein that is a topological barrier must also be a steric barrier to RNAPs.')

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
		self.can_be_displaced_at_TSS_by_RNAP = can_be_displaced_at_TSS_by_RNAP

class Model: # Class to hold the model, including genomic setup, model setup, and dynamic state variables
	def __init__(self, genomic_setup: GenomicSetup, model_setup: ModelSetup, binding_proteins: list[BindingProtein] = None) -> None:
		self.genomic_setup = genomic_setup # GenomicSetup object
		self.model_setup = model_setup # ModelSetup object
		self.x_dict = [[] for _ in genomic_setup.gene_names] # List of lists to hold positions of RNAPs for each gene
		self.theta_dict = [[] for _ in genomic_setup.gene_names] # List of lists to hold angular positions of RNAPs for each gene
		self.Lk = [(genomic_setup.clamp_right - genomic_setup.clamp_left) / model_setup.h_dna] # List to hold linking number of each DNA segment; initially, only one segment spanning the entire DNA
		self.promoter_status = [1 for _ in genomic_setup.gene_names] if genomic_setup.promoter_mode == 'constitutive' else [0 for _ in genomic_setup.gene_names] # List to hold promoter status (1: ON, 0: OFF) for each gene
		self.mRNA_counts = [0 for _ in genomic_setup.gene_names] # List to hold mRNA counts for each gene
		
		if binding_proteins is None:
			binding_proteins = []
		self.binding_proteins = binding_proteins # List of BindingProtein objects representing other DNA-binding proteins in the system
		if genomic_setup.chromatin_type == 'eukaryotic':
			nucl_count = genomic_setup.get_total_nucleosome_count()
			nucleosomes = BindingProtein(protein_name = 'nucleosome', total_copy_number = nucl_count, is_steric_barrier_to_RNAPs = genomic_setup.nucleosomes_are_steric_barriers_to_RNAPs, is_topological_barrier = False, basal_on_rate = 1.2 / (genomic_setup.clamp_right - genomic_setup.clamp_left), basal_off_rate = 0.4, is_a_nucleosome = True, can_be_displaced_at_TSS_by_RNAP = genomic_setup.nucleosomes_can_be_displaced_at_TSS_by_RNAP, on_rate_func = genomic_setup.nucleosome_on_rate_func, off_rate_func = genomic_setup.nucleosome_off_rate_func)
			self.binding_proteins = [nucleosomes] + self.binding_proteins
		self.binding_proteins_positions = [[] for _ in self.binding_proteins] # List of lists to hold positions of each bound protein; each sublist corresponds to a binding protein type and contains the positions of all bound proteins of that type
	
	def print_model_setup(self) -> None: # Utility function to print model setup information
		self.genomic_setup.print_genomic_setup()
		print('Genome binding proteins:')
		print('Name\tTotal copy number\tSteric barrier to RNAPs\tTopological barrier\tIs a nucleosome\tCan be displaced at TSS by RNAP')
		print('-' * 80)
		for i in range(len(self.binding_proteins)):
			protein = self.binding_proteins[i]
			print(f'{protein.protein_name}\t{protein.total_copy_number}\t{protein.is_steric_barrier_to_RNAPs}\t{protein.is_topological_barrier}\t{protein.is_a_nucleosome}\t{protein.can_be_displaced_at_TSS_by_RNAP}')
		print('=' * 80)

class SimulationSetupAndState: # Class to hold simulation setup parameters
	def __init__(self, simulation_end_mode: int, simulation_end_criterion: Union[float, list[int]], integration_method: str = 'RK23', integration_time_resolution: float = 1.0e-1, integration_rtol: float = 1.0e-6, integration_atol: float = 1.0e-8, RNAP_alive_status_check_interval: float = 1.0, max_RNAPs_to_recruit: list[int] = None, Gillespie_random_seed: int = 42, everything_else_random_seed: int = 42) -> None:
		self.simulation_end_mode = simulation_end_mode # 0: time-based, 1: event-based
		assert simulation_end_mode in [0, 1], 'simulation_end_mode must be either 0 (time-based) or 1 (event-based).'

		self.simulation_end_criterion = simulation_end_criterion # Float for time-based (in s) or list of integers for event-based (number of RNAPs that must finish transcription for each gene)
		if simulation_end_mode == 0:
			self.simulation_end_time = float(simulation_end_criterion)
		else:
			if not isinstance(simulation_end_criterion, list):
				raise ValueError('For simulation_end_mode 1 (i.e., event-based), simulation_end_criterion must be a list of integers representing event counts for each gene.')
			self.simulation_end_event_counts = list(simulation_end_criterion)
		
		self.integration_method = integration_method
		assert integration_method in ['RK23', 'RK45', 'DOP853', 'Radau', 'BDF', 'LSODA'], 'integration_method must be one of: RK23, RK45, DOP853, Radau, BDF, LSODA.'
		if integration_method != 'RK23':
			warnings.warn('Using integration method "' + integration_method + '". Solvers other than RK23 may produce non-physical intermediate states that violate steric constraints, causing the simulation to crash.')

		self.integration_time_resolution = integration_time_resolution # Time resolution for integration (in s)
		assert integration_time_resolution > 0.0, 'integration_time_resolution must be a positive float.'

		self.integration_rtol = integration_rtol # Relative tolerance for ODE integration
		assert integration_rtol > 0.0, 'integration_rtol must be a positive float.'
		self.integration_atol = integration_atol # Absolute tolerance for ODE integration
		assert integration_atol > 0.0, 'integration_atol must be a positive float.'

		self.RNAP_alive_status_check_interval = RNAP_alive_status_check_interval # Interval for checking RNAP alive status (in s)
		assert RNAP_alive_status_check_interval > 0.0, 'RNAP_alive_status_check_interval must be a positive float.'

		if max_RNAPs_to_recruit is not None:
			assert isinstance(max_RNAPs_to_recruit, list), 'max_RNAPs_to_recruit must be a list of integers representing the maximum number of RNAPs to recruit for each gene.'
		self.max_RNAPs_to_recruit = max_RNAPs_to_recruit # Maximum number of RNAPs to recruit for each gene

		self.curr_simulation_time = 0.0 # Current simulation time (in s)
		self.last_event_index = -1 # Index of the last event that occurred in the simulation
		self.last_event_type = None # Type of the last event that occurred in the simulation
		self.simulation_completed = False # Flag indicating whether the simulation has completed

		self.state_has_been_initialized = False # Flag indicating whether the simulation state has been initialized; used to ensure that setup_simulation_state is called before running the simulation

		self.Gillespie_random_seed = Gillespie_random_seed # Random seed for Gillespie events
		self.everything_else_random_seed = everything_else_random_seed # Random seed for all other stochastic processes in the simulation (e.g., choosing segments for supercoiling relaxation events, choosing which bound protein unbinds in a binding protein unbinding event, etc.)

		self.rng_Gillespie = random.Random(self.Gillespie_random_seed)
		self.rng_everything_else = random.Random(self.everything_else_random_seed)
	
	def setup_simulation_state(self, genomic_setup: GenomicSetup) -> None:
		if self.simulation_end_mode == 1:
			if len(self.simulation_end_event_counts) != len(genomic_setup.gene_names):
				raise ValueError('Length of simulation_end_criterion list must match the number of genes in genomic_setup.')
		if self.max_RNAPs_to_recruit is not None:
			if len(self.max_RNAPs_to_recruit) != len(genomic_setup.gene_names):
				raise ValueError('Length of max_RNAPs_to_recruit list must match the number of genes in genomic_setup.')
			if self.simulation_end_mode == 1:
				for i in range(len(self.simulation_end_event_counts)):
					if self.simulation_end_event_counts[i] > self.max_RNAPs_to_recruit[i]:
						raise ValueError(f'For gene index {i}, number of RNAPs that must finish transcription (simulation_end_event_counts) cannot be greater than max_RNAPs_to_recruit ({self.max_RNAPs_to_recruit[i]}).')
		
		self.RNAPs_finished_transcription = [0 for _ in genomic_setup.gene_names]
		self.RNAPs_exit_positions = [[] for _ in genomic_setup.gene_names]
		self.RNAP_recruitment_times = [[] for _ in genomic_setup.gene_names]
		self.RNAP_exit_times = [[] for _ in genomic_setup.gene_names]

		self.state_has_been_initialized = True
	
	def update_simulation_end_criterion(self, new_criterion: Union[float, list[int]]) -> None:
		if self.simulation_end_mode == 0:
			new_simulation_end_time = float(new_criterion)
			if new_simulation_end_time < self.curr_simulation_time:
				raise ValueError('New simulation end time cannot be earlier than the current simulation time.')
			if new_simulation_end_time < self.simulation_end_time:
				warnings.warn('New simulation end time is earlier than the previous end time. The simulation will now end at the earlier time.')
			self.simulation_end_time = new_simulation_end_time
			self.simulation_end_criterion = new_simulation_end_time
		elif self.simulation_end_mode == 1:
			if not isinstance(new_criterion, list):
				raise ValueError('For simulation_end_mode 1 (i.e., event-based), new_criterion must be a list of integers representing event counts for each gene.')
			if len(new_criterion) != len(self.simulation_end_event_counts):
				raise ValueError('Length of new_criterion list must match the number of genes in the simulation.')
			if self.max_RNAPs_to_recruit is not None:
				for i in range(len(new_criterion)):
					if new_criterion[i] > self.max_RNAPs_to_recruit[i]:
						raise ValueError('For any gene, number of RNAPs that must finish transcription (per the new simulation end criterion) cannot be greater than max_RNAPs_to_recruit (specified earlier).')
			invalid_criterion = True
			for i in range(len(new_criterion)):
				if new_criterion[i] > self.RNAPs_finished_transcription[i]:
					invalid_criterion = False
					break
			if invalid_criterion:
				raise ValueError('New simulation end criterion must require at least one more RNAP to finish transcription for at least one gene.')
			self.simulation_end_event_counts = [new_criterion[i] for i in range(len(new_criterion))]
			self.simulation_end_criterion = list(new_criterion)
		self.simulation_completed = False # Allow the simulation to resume under the updated criterion
	
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

def save_simulation_state_to_file(model: Model, simulation_setup_and_state: SimulationSetupAndState, filename: str) -> None: # Save the current state of the simulation to a file
	simulation_state = {
		'model': model,
		'simulation_setup_and_state': simulation_setup_and_state
	}
	with open(filename, 'wb') as f:
		dill.dump(simulation_state, f)

def load_simulation_state_from_file(filename: str) -> tuple[Model, SimulationSetupAndState]: # Load the simulation state from a file and return the model and simulation setup/state
	with open(filename, 'rb') as f:
		simulation_state = dill.load(f)
	if not simulation_state['simulation_setup_and_state'].state_has_been_initialized:
		raise ValueError('The loaded simulation state has not been initialized. Cannot resume simulation from this state.')
	return simulation_state['model'], simulation_state['simulation_setup_and_state']