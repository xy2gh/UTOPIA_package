# utopia

Package to run the UTOPIA unit world model to estimate mass and particle distribution of microplastics in a unit world model composed of 17 compartments covering air, water, soil and sediments. Results are given for 5 size fractions ranging from hundreds of nanometers to milimiters and 4 aggregation states. Exposure indicators such as persistence and characteristic travel distance are also given as outputs.

## Usage


Before installing and testing the package, it is highly recommended to set up a virtual environment, this can be made with conda as follows:

```bash
$ conda create --name utopia_env python=3.9 -y
```

Activate the environment:

```bash
$ conda activate utopia_env
```
Use poetry to install our package using the command poetry install at the command line from the root package directory:

```bash
$ poetry install
```

Install the model with `pip` (not ready? use poetry above):

```bash
$ pip install utopia
```

Then run the `utopia` model with example data and plot the results:

```python
from utopia.utopia import utopiaModel
from utopia.results_processing.process_results import*
from utopia.results_processing.mass_balance_check import*



# Create the model, pass it example config and data, then run it
config_data = utopiaModel.load_json_file("data/default_config.json")
data_data = utopiaModel.load_json_file("data/default_data.json")
model = utopiaModel(config=config_data, data=data_data)
model.summarize()
model.run()

# Process results and print and plot them
massBalance(model)
processor = ResultsProcessor(model)
processor.process_all()  # Process all results

# Print exposure indicators
processor.processed_results["Overall_exposure_indicators"]
processor.processed_results["size_fraction_indicators"]
pd.DataFrame(processor.processed_results["emission_fractions_mass_data"])
```

[Access the user step by step guide here.](https://github.com/PradoDomercq/utopia_package/blob/main/docs/model_tutorial.ipynb) 

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`utopia` was created by Prado Domercq. It is licensed under the terms of the MIT license.

## Acknowledgements

Thanks to the European Chemical Industry Council Long-Range Research Initiative ([Cefic-LRI](https://cefic-lri.org/)) for providing funding for this work, under project number ECO56.
