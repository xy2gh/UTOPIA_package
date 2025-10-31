# utopia_package

Package to run the [UTOPIA unit world model](https://github.com/microplastics-cluster/UTOPIA_model), a process-based mass balance model that estimates mass and particle number distributions of microplastics in a unit world model composed of 17 compartments covering air, water, soil and sediments though steady state solution of a 340 coupled first-order differential equations. Results are given for 5 size fractions ranging from hundreds of nanometers to milimiters and 4 aggregation states. Exposure indicators such as overall persistence and residence time are also given as outputs. Sensitivity and uncertainty analysis are possible using the [Monaco Monte Carlo framework](https://github.com/microplastics-cluster/utopia_package/blob/main/docs/montecarlo_tutorial.ipynb)

This package version is based on the original UTOPIA model (see [https://github.com/microplastics-cluster/UTOPIA_model]). It has been refactored for package use and includes updated process formulations and parameterizations. The package version also includes shape functionalities where users can run the model for [spheres](https://github.com/microplastics-cluster/utopia_package/blob/main/docs/model_tutorial_sphere.ipynb) and [fibers](https://github.com/microplastics-cluster/utopia_package/blob/main/docs/model_tutorial_fiber.ipynb).

## Installation

Before installing and testing the package, it is highly recommended to set up a virtual environment, this can be made as follows:

1. Create and activate a python environment
   
```bash
python -m venv .venv
.venv\Scripts\activate
```

2. Install the package

You can install directly once the package is published to PyPI (not ready? use poetry):

```bash
pip install utopia-pkg
```
3. To verify the installation:
   
```python
import utopia_pkg
print(utopia_pkg.__version__)
```

## Quick Start

You can run the UTOPIA model either:

- From Jupyter notebooks (recommended for new users — see docs/), or
- Using the provided example script (docs/run_utopia_example.py).
  
### Run the example python script

After installation, you can execute the model using the provided example script:

```bash
poetry run python docs/run_utopia_example.py
```

This script:

- Loads example configuration and data (data/default_config.json and data/default_data.json),
- Initializes and runs the UTOPIA model,
- Checks the mass balance (Difference inflow-outflow),
- Processes the results,
- Prints exposure indicators and mass and particle number distribution figures.

If successful, you should see printed summaries of model results in your terminal and output files printed such as the mentioned plots.

### Run the model interactively using the step by step guide jupyter notebook

The notebooks in docs/ provide guided examples that explain each step, configuration parameter, and output indicator in detail.

Explore the following tutorials:

-[Generic utopia step by step guide](https://github.com/microplastics-cluster/utopia_package/blob/main/docs/model_tutorial.ipynb) 
-[Fiber specific step by step guide](https://github.com/microplastics-cluster/utopia_package/blob/main/docs/model_tutorial_fiber.ipynb)
-[Sensitivity and Uncertainty analysis tutorial](https://github.com/microplastics-cluster/utopia_package/blob/main/docs/montecarlo_tutorial.ipynb)

## Contributing

Interested in contributing? Check out the contributing guidelines. Please note that this project is released with a Code of Conduct. By contributing to this project, you agree to abide by its terms.

## License and authorship

Developed by [Prado Domercq](https://github.com/PradoDomercq) with contributions from [Xiaoyu Zhang](https://github.com/xy2gh)
Transferred to [microplastics-cluster](https://github.com/microplastics-cluster) for continued collaborative development and maintenance
`utopia` is licensed under the terms of the MIT license.


## Acknowledgements

Thanks to the European Chemical Industry Council Long-Range Research Initiative ([Cefic-LRI](https://cefic-lri.org/)) for providing funding for this work, under project number ECO56.

