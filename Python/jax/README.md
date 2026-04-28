# modular_Rothe (JAX)

JAX + SciPy implementation of Gaussian-basis Rothe time stepping for the time-dependent Schrodinger equation (TDSE).

The solver represents the wavefunction in a generally non-orthogonal Gaussian basis and advances one implicit Rothe step at a time:

- nonlinear update of Gaussian parameters,
- linear coefficient solve in the updated basis,
- optional HDF5 checkpointing for resume and analysis.

## Current capabilities

- Core solver in `rothe/solver.py` with optional frozen basis functions.
- JIT objective and gradient construction in `rothe/objective.py`.
- Shared block-matrix assembly helpers in `rothe/block_assembly.py`.
- Shared regularized linear solve/scoring helpers in `rothe/linear_solve.py`.
- HDF5 append-only per-step persistence in `rothe/io.py`.
- System helpers for Hydrogen, Henon-Heiles, and Double Well examples under `rothe/systems/` and `scripts/`.

## Installation

From this directory:

```bash
pip install -e .
```

For tests and plotting extras:

```bash
pip install -e ".[test,plot]"
```

Note: `pyproject.toml` currently lists `jax[cuda12]`. If you run on CPU-only machines, install the JAX variant that matches your platform before running examples.

## Quick start

All propagators require an explicit regularization strength via `--regularization_lambda`.

### Hydrogen

Installed console entry point:

```bash
propagate-hydrogen \
  --num_gauss_wavefunction 21 \
  --num_gauss_potential 21 \
  --epsilon 1e-2 \
  --regularization_lambda 1e-8
```

Useful options:

- `--splitting_type {none,kinetic}`
- `--dt 0.2` (real-time) or `--dt 0.2j` (imaginary-time)
- `--frozen` / `--no-frozen`
- `--num_dynamic 10`
- `--t_start <time>` to resume from an existing file

### Henon-Heiles

Installed console entry point:

```bash
propagate-henon-heiles \
  --dim 2 \
  --dt 0.01 \
  --num_gaussians 5 \
  --epsilon 1e-3 \
  --regularization_lambda 1e-8
```

### 1D Double Well

Module runner (no console entry point in `pyproject.toml`):

```bash
python -m scripts.propagate_double_well \
  --num_gaussians 3 \
  --epsilon 1e-3 \
  --regularization_lambda 1e-8
```

## HDF5 output format

Saved by `rothe.io.save_rothe_state` to:

- `wave_function_data/{sim_name}__{splitting_type}.h5` by default.

File-level layout:

- root attrs: `sim_name`, `splitting_type`, `dt`, `epsilon`, `regularization_lambda`, `created`
- scalar dataset: `polynomial_string`
- group: `steps/`
  - `steps/0000/`, `steps/0001/`, ...
  - each step stores exact-shape arrays for that time point:
    - `params` with shape `(n_total, D, 4)`
    - `coeffs` with shape `(n_total,)`
  - step attrs include `t`, `rothe_error`, `n_dynamic`
  - optional step data: `rng_state` (attr), `dipole` (dataset), `double_autocorrelation` (attr)

This schema allows variable Gaussian counts across time steps without padding.

## Repository map

- `rothe/solver.py`: high-level `RotheSolver`, optimization loop, candidate augmentation.
- `rothe/objective.py`: Rothe objective assembly, value/gradient factories.
- `rothe/block_assembly.py`: shared formulas for A, rho, B and block layouts.
- `rothe/linear_solve.py`: shared regularized linear solves (JAX and NumPy paths).
- `rothe/optimization.py`: objective wrappers, parameter transforms, optimizer helpers.
- `rothe/io.py`: checkpoint save/load/resume-trim helpers.
- `rothe/wavefunction.py`: Gaussian overlap/H/H2 construction and propagation helpers.
- `rothe/systems/`: problem-specific setup utilities.
- `scripts/`: runnable propagation entry points.
- `tests/`: pytest suite.

## Testing

Run from repository root:

```bash
pytest
```

## Notes and limitations

- First call overhead can be noticeable due to JAX compilation.
- SciPy optimization is host-side; performance depends on objective/gradient call cost and basis size.
- Physical units/conventions are controlled by the chosen potential strings and system setup.
