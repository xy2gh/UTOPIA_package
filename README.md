# utopia

Package to run the UTOPIA unit world model to estimate mass and particle distribution of microplastics in a unit world model composed of 17 compartments covering air, water, soil and sediments. Results are given for 5 size fractions ranging from hundreds of nanometers to milimiters and 4 aggregation states. Exposure indicators such as persistence and characteristic travel distance are also given as outputs.

## Usage

Install the model with `pip`:

```bash
$ pip install utopia
```

Then run the `utopia` model with example data and plot the results:

```python
from utopia import utopiaModel 
from results_processing.mass_balance_check import*
from results_processing.process_results import*

# Create the model, pass it example config and data, then run it
config_data = utopiaModel.load_json_file("data/default_config.json")
data_data = utopiaModel.load_json_file("data/default_data.json")
model = utopiaModel(config=config_data, data=data_data)
model.summarize()
model.run()

# Process results and print and plot them
massBalance(model)
processor = ResultsProcessor(model)
processor.estimate_flows()
processor.generate_flows_dict()
processor.process_results()
for fraction in ["mass_fraction", "number_fraction"]:
    processor.plot_fractionDistribution_heatmaps(fraction)

processor.extract_results_by_compartment()
for fraction in ["%_mass", "%_number"]:
    processor.plot_compartment_distribution(fraction)
processor.results_by_comp

for i in range(len(processor.results_by_comp)):
    emissions=sum(processor.model.emiss_dict_g_s[processor.results_by_comp['Compartments'].iloc[i]].values())
    print(f"Mass balance for {processor.results_by_comp['Compartments'].iloc[i]}: {processor.results_by_comp['Total_inflows_g_s'].iloc[i]+emissions-processor.results_by_comp['Total_outflows_g_s'].iloc[i]}")

# Calculate exposure indicators
processor.estimate_exposure_indicators()
processor.overall_exposure_indicators
processor.size_fraction_indicators
processor.estimate_emission_fractions()
```

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

`utopia` was created by Prado Domercq. It is licensed under the terms of the MIT license.

## Credits

`utopia` was created with [`cookiecutter`](https://cookiecutter.readthedocs.io/en/latest/) and the `py-pkgs-cookiecutter` [template](https://github.com/py-pkgs/py-pkgs-cookiecutter).

## Acknowledgements

Thanks to the European Chemical Industry Council Long-Range Research Initiative ([Cefic-LRI](https://cefic-lri.org/)) for providing funding for this work, under project number ECO56.
