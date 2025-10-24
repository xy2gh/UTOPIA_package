# Import the model and helper functions
from utopia.utopia import utopiaModel
from utopia.results_processing.mass_balance_check import *
from utopia.results_processing.process_results import *

# Load example configuration and input data
config_data = utopiaModel.load_json_file("data/default_config.json")
data_data = utopiaModel.load_json_file("data/default_data.json")

print("Loaded Configuration Data:", config_data)
print("Loaded Input Data:", data_data)

# Initialize and run the model
model = utopiaModel(config=config_data, data=data_data)
model.summarize()
model.run()

# Process and check results
massBalance(model)
processor = ResultsProcessor(model)
processor.process_all()

# Display processed outputs
print(processor.processed_results["Overall_exposure_indicators"])
print(processor.processed_results["size_fraction_indicators"])

pd.DataFrame(processor.processed_results["emission_fractions_mass_data"])
