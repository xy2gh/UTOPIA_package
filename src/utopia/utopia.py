import copy
import json
import os
from datetime import datetime
import math
import numpy as np
from objects_generation import generate_objects


class utopia():
    """The class that controls usage of the UTOPIA model
    Parameters
    ----------
    config : dict
        Configuration dictionary with model configuration options regarding to parameterization of the unit world compartments
    data : dict
        Data dictionary with input data required to run the model. This contains the selected options from the imput menues of the user interface including microplastics properties and emission scenario selected (in the future to also include compartment properties or  the selected compartment parameterization scenario).
    validate : bool. deafult=True
        Validate input data and configuration options
    """

    def __init__(self, config: None, data: None, validate=True):
        """Initialize the UTOPIA model"""
        
        # Load default config and data if not provided
        if config is None:
            config = self.load_json_file("data/default_config.json")
        if data is None:
            data = self.load_json_file("data/default_data.json")
        
        # Validate the config and data inputs
        if validate:
            config, data = self.validate_inputs(config, data)
        
        # Assign attributes
        self.config = config
        self.data = data

        # Load parameters
        self.load_parameters()

    @staticmethod
    def load_json_file(filename):
        """Load a JSON file from the package data directory."""
        filepath = os.path.join(os.path.dirname(__file__), filename)
        with open(filepath, "r") as file:
            return json.load(file)
        
    @staticmethod
    def check_required_keys(dictionary, required_keys, dict_name):
        """Checks if all required keys are present in a dictionary."""
        missing_keys = [key for key in required_keys if key not in dictionary]
        if missing_keys:
            raise KeyError(f"Missing keys in {dict_name}: {', '.join(missing_keys)}")

    def validate_inputs(self, config, data):
        """Validates input dictionaries to ensure all required keys exist."""
        required_data_keys = [
            "MPdensity_kg_m3", "MP_composition", "shape", "FI", "t_half_deg_free",
            "heter_deg_factor", "biof_deg_factor", "factor_deepWater_soilSurface",
            "factor_sediment", "t_frag_gen_FreeSurfaceWater", "biof_frag_factor",
            "heter_frag_factor", "emiss_dict_g_s"
        ]

        required_config_keys = [
            "big_bin_diameter_um", "N_sizeBins", "vol_algal_cell_m3",
            "spm_density_kg_m3", "comp_input_file_name",
            "comp_interactions_file_name", "boxName", "compartment_names",
            "compartment_types", "MPforms_list"
        ]

        self.check_required_keys(data, required_data_keys, "data")
        self.check_required_keys(config, required_config_keys, "config")

    def load_parameters(self):
        """Loads required parameters from config and data dictionaries."""
        # Microplastics physical properties
        self.MPdensity_kg_m3 = self.data["MPdensity_kg_m3"]
        self.MP_composition = self.data["MP_composition"]
        self.shape = self.data["shape"]
        self.big_bin_diameter_um = self.config["big_bin_diameter_um"]
        self.N_sizeBins = self.config["N_sizeBins"]
        self.FI = self.data["FI"]
        self.t_half_deg_free = self.data["t_half_deg_free"]
        self.t_frag_gen_FreeSurfaceWater = self.data["t_frag_gen_FreeSurfaceWater"]
        
        self.heter_deg_factor = self.data["heter_deg_factor"]
        self.biof_deg_factor = self.data["biof_deg_factor"]
        self.factor_deepWater_soilSurface = self.data["factor_deepWater_soilSurface"]
        self.factor_sediment = self.data["factor_sediment"]
        
        self.biof_frag_factor = self.data["biof_frag_factor"]
        self.heter_frag_factor = self.data["heter_frag_factor"]
        
        self.emiss_dict_g_s = self.data["emiss_dict_g_s"]

        # Environmental characteristics
        
        self.vol_algal_cell_m3 = self.config["vol_algal_cell_m3"]
        self.spm_density_kg_m3 = self.config["spm_density_kg_m3"]
        self.comp_input_file_name = self.config["comp_input_file_name"]
        self.comp_interactions_file_name = self.config["comp_interactions_file_name"]
        self.boxName = self.config["boxName"]
        self.compartment_names = self.config["compartment_names"]
        self.compartment_types = self.config["compartment_types"]
        self.MPforms_list = self.config["MPforms_list"]

        # Derived environmental parameters
        self.radius_algae_m = ((3.0 / 4.0) * (self.vol_algal_cell_m3 / math.pi)) ** (1.0 / 3.0)
        self.spm_radius_um = self.vol_algal_cell_m3 * 1e6

    def run(self):
        """Runs the UTOPIA model. This is a placeholder for model execution logic."""
        print("Running UTOPIA model with configured parameters...")

    def summarize(self):
        """Prints a summary of the model's key parameters."""
        print(f"Model: UTOPIA")
        print(f"Box Name: {self.boxName}")
        print(f"Microplastic Density (kg/m3): {self.MPdensity_kg_m3}")
        print(f"Number of Compartments: {len(self.compartment_names)}")


        

#         # Define model inputs path 
#         inputs_path = os.path.join(os.path.dirname(__file__), "data")

#         # Generate model objects        
#         (
#             system_particle_object_list,
#             SpeciesList,
#             spm,
#             dict_comp,
#             model_lists,
#             particles_df,
#         ) = generate_objects(
#             inputs_path,
#             boxName=self.boxName,
#             MPforms_list=self.MPforms_list,
#             comp_input_file_name=self.comp_input_file_name,
#             comp_interactFile_name=self.comp_interactFile_name,
#             mp_imputFile_name=self.mp_imputFile_name,
#             spm_radius_um=self.spm_radius_um,
#             spm_density_kg_m3=self.spm_density_kg_m3,
#         )
# def run(self)