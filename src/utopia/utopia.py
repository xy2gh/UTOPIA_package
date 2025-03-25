import copy
import string
import json
import pandas as pd
import os
from datetime import datetime
import math
import numpy as np
from objects_generation import generate_objects


class utopia:
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
            self.validate_inputs(config, data)

        # Assign attributes
        self.config = config
        self.data = data

        # Load parameters
        self.load_parameters()
        self.particles_df = self.generate_particles_dataframe()
        self.generate_size_dictionaries()

    @staticmethod
    def load_json_file(filename):
        """Load a JSON file from the package data directory."""
        filepath = os.path.join(os.path.dirname(__file__), filename)
        with open(filepath, "r") as file:
            return json.load(file)

    @staticmethod
    def load_csv_column(file_path, column_name):
        """Load a column from input CSV file: Reads a single column from a CSV file and returns it as a list"""
        df = pd.read_csv(file_path, usecols=[column_name])
        return df[column_name].tolist()

    @staticmethod
    def check_required_keys(dictionary, required_keys, dict_name):
        """Checks if all required keys are present in a dictionary."""
        missing_keys = [key for key in required_keys if key not in dictionary]
        if missing_keys:
            raise KeyError(f"Missing keys in {dict_name}: {', '.join(missing_keys)}")

    def validate_inputs(self, config, data):
        """Validates input dictionaries to ensure all required keys exist."""
        required_data_keys = [
            "MPdensity_kg_m3",
            "MP_composition",
            "shape",
            "FI",
            "t_half_deg_free",
            "heter_deg_factor",
            "biof_deg_factor",
            "factor_deepWater_soilSurface",
            "factor_sediment",
            "t_frag_gen_FreeSurfaceWater",
            "biof_frag_factor",
            "heter_frag_factor",
            "emiss_dict_g_s",
        ]

        required_config_keys = [
            "big_bin_diameter_um",
            "N_sizeBins",
            "vol_algal_cell_m3",
            "spm_density_kg_m3",
            "comp_input_file_name",
            "comp_interactions_file_name",
            "boxName",
            "MPforms_list",
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

        # Environmental characteristics

        self.vol_algal_cell_m3 = self.config["vol_algal_cell_m3"]
        self.spm_density_kg_m3 = self.config["spm_density_kg_m3"]
        self.comp_input_file_name = self.config["comp_input_file_name"]
        self.comp_interactions_file_name = self.config["comp_interactions_file_name"]
        self.boxName = self.config["boxName"]

        self.MPforms_list = self.config["MPforms_list"]

        # Load parameters from config and data dictionaries
        self.MPforms_list = self.config["MPforms_list"]
        self.compartments_list = self.load_csv_column(
            f"data/{self.comp_input_file_name}", "Cname"
        )

        # Derived environmental parameters
        self.radius_algae_m = ((3.0 / 4.0) * (self.vol_algal_cell_m3 / math.pi)) ** (
            1.0 / 3.0
        )
        self.spm_radius_um = self.vol_algal_cell_m3 * 1e6

        # Emission scenario
        self.emiss_dict_g_s = self.data["emiss_dict_g_s"]

    def generate_particles_dataframe(self):
        """Generates the microplastics input DataFrame from Utopia model attributes."""
        MPdensity_kg_m3 = self.MPdensity_kg_m3
        shape = self.shape
        N_sizeBins = self.N_sizeBins
        big_bin_diameter_um = self.big_bin_diameter_um

        # Generate size distribution
        size_distribution = [big_bin_diameter_um]
        for _ in range(N_sizeBins - 1):
            size_distribution.append(size_distribution[-1] / 10)
        size_distribution.reverse()

        # Only supports spherical particles for now
        if shape == "sphere":
            data = {
                "Name": [f"mp{i+1}" for i in range(N_sizeBins)],
                "form": ["freeMP"] * N_sizeBins,
                "shape": [shape] * N_sizeBins,
                "composition": [self.MP_composition] * N_sizeBins,
                "density_kg_m3": [MPdensity_kg_m3] * N_sizeBins,
                "dimensionX_um": [d / 2 for d in size_distribution],
                "dimensionY_um": [d / 2 for d in size_distribution],
                "dimensionZ_um": [d / 2 for d in size_distribution],
            }
            return pd.DataFrame(data)
        else:
            raise ValueError("Shape not supported yet")

    def generate_size_dictionaries(self):
        """Generates size-related dictionaries as attributes."""
        if self.particles_df is None:
            raise ValueError("Particles DataFrame has not been generated.")

        # Dictionary mapping particle names to sizes
        self.dict_size_coding = dict(
            zip(self.particles_df["Name"], self.particles_df["dimensionX_um"] * 2)
        )

        # Generate size codes (a-z based on number of bins)
        self.size_codes = list(string.ascii_lowercase[: self.N_sizeBins])

        # Dictionary mapping size codes to sizes
        self.size_dict = dict(zip(self.size_codes, self.dict_size_coding.values()))

    def run(self):
        """Runs the UTOPIA model. This is a placeholder for model execution logic."""
        print("Running UTOPIA model with configured parameters...")

    def summarize(self):
        """Prints a summary of the model's key parameters."""
        print(f"Model: UTOPIA")
        print(f"Box Name: {self.boxName}")
        print(f"Microplastic Density (kg/m3): {self.MPdensity_kg_m3}")
        print("MP shape: ", self.shape)

        def print_fragmentation_style(F):
            if F == 0:
                print("Fragmentation style: Erosive")
            elif F == 1:
                print("Fragmentation style: Sequential")
            else:
                print(f"Fragmentation style: Mixed (F = {F})")

        print_fragmentation_style(self.FI)
        print("Fragmetation timescale (days): ", self.t_frag_gen_FreeSurfaceWater)
        print("Discorporation timescale (days): ", self.t_half_deg_free)
        for compartment, size_fractions in self.emiss_dict_g_s.items():
            for fraction, value in size_fractions.items():
                if value > 0:
                    print(
                        f"Emissions to {compartment} for size fraction {self.size_dict[fraction]} {chr(181)}m: {value} g/s"
                    )


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
