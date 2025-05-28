import copy
import os
from datetime import datetime
import math


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

    def __init__(self, config: dict, data: dict, validate: bool = True) -> None:
        """Initialize the UTOPIA model"""
        # Validate the config and data inputs
        if validate:
            config, data = self.validate_inputs(config, data)
        # If validation is sucessful, set the config and data attributes
        self.config = config
        self.data = data

        # Define inputs parameters
        # Microplastics physical properties
        self.MPdensity_kg_m3 = self.data["MPdensity_kg_m3"]
        self.MP_composition = self.data["MP_composition"]
        self.shape = self.data["shape"]

        self.big_bin_diameter_um = self.config["big_bin_diameter_um"]
        self.N_sizeBins = self.config["N_sizeBins"]

        # Enviromental characteristics
        # suspended particulate matter properties
        self.radius_algae_m = (
            (3.0 / 4.0) * (self.config["vol_algal_cell_m3"] / math.pi)
        ) ** (1.0 / 3.0)
        self.spm_radius_um = self.config["vol_algal_cell_m3"] * 1e6
        self.spm_density_kg_m3 = self.config["spm_density_kg_m3"]
        # compartments input files
        self.comp_input_file_name = self.config["comp_input_file_name"]
        self.comp_interactions_file_name = self.config["comp_interactions_file_name"]
        # model default configuration
        self.boxName = self.config["boxName"]
        self.compartment_names = self.config["compartment_names"]
        self.compartment_types = self.config["compartment_types"]
        # Microplastics weathering properties and compartment and aggregation state factors
        self.FI = self.data["FI"]
        self.t_half_deg_free = self.data["t_half_deg_free"]
        self.heter_deg_factor = self.data["heter_deg_factor"]
        self.biof_deg_factor = self.data["biof_deg_factor"]
        self.factor_deepWater_soilSurface = self.data["factor_deepWater_soilSurface"]
        self.factor_sediment = self.data["factor_sediment"]
        self.t_frag_gen_FreeSurfaceWater = self.data["t_frag_gen_FreeSurfaceWater"]
        self.biof_frag_factor = self.data["biof_frag_factor"]
        self.heter_frag_factor = self.data["heter_frag_factor"]

        # Load emissions input data
        self.emiss_dict_g_s = self.data["emiss_dict_g_s"]

        # Generate model objects
        (compartment_types)
