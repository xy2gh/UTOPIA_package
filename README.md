# UTOPIA

Package to run the UTOPIA unit world model, a process-based mass balance model that estimates mass and particle number distributions of microplastics in a unit world model composed of 17 compartments covering air, water, soil and sediments though steady state solution of a 340 coupled first-order differential equations. Results are given for 5 size fractions ranging from hundreds of nanometers to milimiters and 4 aggregation states. Exposure indicators such as overall persistence and residence time are also given as outputs. Sensitivity and uncertainty analysis are possible using the [Monaco Monte Carlo framework](https://github.com/microplastics-cluster/utopia_package/blob/main/docs/montecarlo_tutorial.ipynb)

This package version is based on the original UTOPIA model (see [original repo link]). It has been refactored for package use and includes updated process formulations and parameterizations. The package version also includes shape functionalities where users can run the model for [spheres](https://github.com/microplastics-cluster/utopia_package/blob/main/docs/model_tutorial_sphere.ipynb) and [fibers](https://github.com/microplastics-cluster/utopia_package/blob/main/docs/model_tutorial_fiber.ipynb).

## Usage

Before installing and testing the package, it is highly recommended to set up a virtual environment, this can be made with conda as follows:

```bash
$ conda create --name utopia_env python=3.9 -y
```

Activate the environment:

```bash
$ conda activate utopia_env
```

To install Poetry on Windows, run the following command in PowerShell:

```powershell
(Invoke-WebRequest -Uri https://install.python-poetry.org -UseBasicParsing).Content | python -
```

After installation, add the Poetry into your PATH:

```bash
[Environment]::SetEnvironmentVariable("Path", "$env:Path;%USERPROFILE%\.poetry\bin", "User")
```

Then, restart your terminal and verify the installation:

```bash
poetry --version
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

[Access the user step by step guide here.](https://github.com/microplastics-cluster/utopia_package/blob/main/docs/montecarlo_tutorial.ipynb) 

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License

Developed by [Prado Domercq](https://github.com/PradoDomercq).  
Transferred to [microplastics-cluster](https://github.com/microplastics-cluster) for continued collaborative development and maintenance
`utopia` is licensed under the terms of the MIT license.


## Acknowledgements

Thanks to the European Chemical Industry Council Long-Range Research Initiative ([Cefic-LRI](https://cefic-lri.org/)) for providing funding for this work, under project number ECO56.
