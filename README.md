# PyDynSM

PyDynSM is a Python implementation of the **Dynamic Stiffness Matrix (DSM) method** for structural dynamics.

---

## Installation (Users)

This section is for users who want to **install and use PyDynSM**.

### 1. Create the Conda environment


```bash
conda env create -f environment.yml
conda activate pydynsm_env
```

### 2. Install PyDynSM

Install the package using `pip`:

```bash
pip install pydynsm
```

## Development Setup (Developers)

This section is for contributors and developers working on PyDynSM.

### 1. Create the development environment

The development environment includes additional tooling for building, testing, and maintaining the package.

```bash
conda env create -f dev_environment.yml
conda activate pydynsm_dev_env
```

### 2. Install PyDynSM in editable mode


```bash
cd path/to/PyDynSM

pip install -e .
```

This ensures local code changes are immediately reflected without reinstalling.

---

## Building the Package (Developers)

To build source and wheel distributions:

```bash
python -m build
```

The built artifacts will be placed in the `dist/` directory.

---



