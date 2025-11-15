# PlasticSim Python Analysis

Python project for iterative analysis using Jupyter notebooks. Shared code modules are organized in the `src/` directory.

## Setup

1. **Create and activate the conda environment:**
   ```bash
   conda env create -f environment.yml
   conda activate plasticsim-python
   ```

2. **Install the package in editable mode (optional, for importing from src/):**
   ```bash
   pip install -e .
   ```

3. **Start Jupyter:**
   ```bash
   jupyter notebook
   ```

   Or use JupyterLab:
   ```bash
   jupyter lab
   ```

## Project Structure

```
.
├── notebooks/          # Jupyter notebooks for analysis
├── src/               # Shared Python modules
│   └── __init__.py
├── environment.yml    # Conda environment definition
├── requirements.txt   # Pip dependencies (backup/alternative)
├── pyproject.toml     # Package configuration
└── README.md          # This file
```

## Using Shared Modules in Notebooks

To import modules from `src/` in your notebooks, add this at the top:

```python
import sys
sys.path.insert(0, '../src')  # Adjust path if needed

from your_module import your_function
```

Or install the package in editable mode (recommended):
```bash
pip install -e .
```

Then you can import directly:
```python
from src.your_module import your_function
```

## Dependency Management

- All dependencies are tracked in `environment.yml` (conda) and `requirements.txt` (pip backup)
- When adding new conda packages, update `environment.yml`:
  ```bash
  conda env export > environment.yml
  ```
- When adding pip-only packages, add them to the `pip:` section in `environment.yml`
- To update the environment after changes:
  ```bash
  conda env update -f environment.yml --prune
  ```
- For maximum reproducibility, consider pinning exact versions in `environment.yml`

