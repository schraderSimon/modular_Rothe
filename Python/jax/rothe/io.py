"""
HDF5 I/O for the Rothe solver.

Provides append-only simulation state persistence and file loading/trimming
for checkpoint-restart workflows.
"""

import json
import os
import time

import h5py
import numpy as np
from attr import dataclass

from rothe._config import PROJECT_ROOT


@dataclass
class OutputConfig:
    """Configuration for saving solver output."""

    name: None | str = None
    polynomial_string: None | str = None
    out_dir: str = str(PROJECT_ROOT / "wave_function_data")
    compression: str = "gzip"
    compression_opts: int = 4
    epsilon: float = 0.0


def _resolve_output_filename(sim_name, splitting_type, path=None, create_dir=True):
    if path is None:
        path = str(PROJECT_ROOT / "wave_function_data")
    if create_dir:
        os.makedirs(path, exist_ok=True)
    return os.path.join(path, f"{sim_name}__{splitting_type}.h5")


def _compression_kwargs(compression, compression_opts):
    return dict(
        compression=compression if compression else None,
        compression_opts=compression_opts if compression else None,
        shuffle=bool(compression),
    )


def _write_metadata(f, sim_name, splitting_type, dt_scalar, params_np, coeffs_np, epsilon):
    f.attrs.setdefault("sim_name", sim_name)
    f.attrs.setdefault("splitting_type", str(splitting_type))
    f.attrs.setdefault("dt", np.complex128(dt_scalar))
    f.attrs.setdefault("epsilon", float(epsilon))
    f.attrs.setdefault("created", time.time())
    f.attrs["p_dtype"] = str(params_np.dtype)
    f.attrs["c_dtype"] = str(coeffs_np.dtype)


def _ensure_polynomial_string_dataset(f, polynomial_string):
    if "polynomial_string" in f:
        return
    dset = f.create_dataset("polynomial_string", shape=(), dtype=h5py.string_dtype(encoding="utf-8"))
    dset[()] = polynomial_string


def _ensure_step_array_dataset(g, name, tail_shape, dtype, kwargs):
    if name in g:
        return g[name]
    return g.create_dataset(
        name, shape=(0,) + tail_shape, maxshape=(None,) + tail_shape, chunks=(1,) + tail_shape, dtype=dtype, **kwargs
    )


def _maybe_grow_step_array_dataset(g, name, new_tail_shape, dtype, kwargs):
    """Return the dataset for ``name``, growing its tail shape if necessary.

    When the number of Gaussians increases (Gaussian addition), the existing
    dataset is read, deleted, recreated with the larger tail shape, and all
    previous steps are written back with zero-padding in the new dimension.
    """
    if name not in g:
        return g.create_dataset(
            name,
            shape=(0,) + new_tail_shape,
            maxshape=(None,) + new_tail_shape,
            chunks=(1,) + new_tail_shape,
            dtype=dtype,
            **kwargs,
        )

    ds = g[name]
    existing_tail = ds.shape[1:]

    if existing_tail == new_tail_shape:
        return ds

    if len(existing_tail) != len(new_tail_shape):
        raise ValueError(f"Dataset '{name}' ndim changed: {existing_tail} vs {new_tail_shape}")
    if existing_tail[1:] != new_tail_shape[1:]:
        raise ValueError(f"Dataset '{name}' non-Gaussian dimensions changed: {existing_tail} vs {new_tail_shape}")

    old_n = existing_tail[0]
    new_n = new_tail_shape[0]
    if new_n < old_n:
        raise ValueError(f"Dataset '{name}' Gaussian count shrank: {old_n} -> {new_n}")

    old_data = ds[...]  # shape (n_steps, old_n, ...)
    n_steps = old_data.shape[0]

    new_data = np.zeros((n_steps,) + new_tail_shape, dtype=old_data.dtype)
    # Copy existing Gaussians into the new array; new ones default to zero
    idx_slice = (slice(None), slice(None, old_n)) + (slice(None),) * (len(new_tail_shape) - 1)
    new_data[idx_slice] = old_data

    del g[name]
    ds_new = g.create_dataset(
        name,
        shape=(n_steps,) + new_tail_shape,
        maxshape=(None,) + new_tail_shape,
        chunks=(1,) + new_tail_shape,
        dtype=dtype,
        **kwargs,
    )
    if n_steps > 0:
        ds_new[...] = new_data
    return ds_new


def _ensure_time_dataset(g, complex_time, kwargs):
    target_dtype = np.complex128 if complex_time else np.float64
    if "t" not in g:
        return g.create_dataset("t", shape=(0,), maxshape=(None,), chunks=(1024,), dtype=target_dtype, **kwargs)

    ds_existing = g["t"]
    if ds_existing.dtype == target_dtype:
        return ds_existing

    data_existing = ds_existing[...]
    del g["t"]
    return g.create_dataset(
        "t", data=data_existing.astype(target_dtype), maxshape=(None,), chunks=(1024,), dtype=target_dtype, **kwargs
    )


def _ensure_scalar_step_dataset(g, name, dtype, kwargs):
    if name in g:
        return g[name]
    return g.create_dataset(name, shape=(0,), maxshape=(None,), chunks=(1024,), dtype=dtype, **kwargs)


def _ensure_rng_state_dataset(g):
    if "rng_state" in g:
        return g["rng_state"]
    return g.create_dataset(
        "rng_state", shape=(0,), maxshape=(None,), chunks=(1024,), dtype=h5py.string_dtype(encoding="utf-8")
    )


def _json_default(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def _encode_rng_state(rng_state):
    return json.dumps(rng_state, default=_json_default, separators=(",", ":"))


def _decode_rng_state(encoded):
    if isinstance(encoded, bytes):
        encoded = encoded.decode("utf-8")
    if encoded is None or encoded == "":
        return None
    return json.loads(encoded)


def _append_optional_dipole(g, idx, dipole_moment, kwargs):
    if dipole_moment is None:
        return

    dipole_np = np.asarray(dipole_moment)
    dipole_shape = tuple(dipole_np.shape)

    if "dipole" not in g:
        g.create_dataset(
            "dipole",
            shape=(idx + 1,) + dipole_shape,
            maxshape=(None,) + dipole_shape,
            chunks=(1,) + dipole_shape,
            dtype=dipole_np.dtype,
            **kwargs,
        )
        g["dipole"][...] = np.nan

    ds_dipole = g["dipole"]
    if ds_dipole.shape[0] < idx + 1:
        ds_dipole.resize(idx + 1, axis=0)
    ds_dipole[idx, ...] = dipole_np


def _append_optional_double_autocorrelation(g, idx, double_autocorrelation, kwargs):
    if double_autocorrelation is None:
        return

    c2 = np.asarray(double_autocorrelation, dtype=np.complex128).item()
    if "double_autocorrelation" not in g:
        g.create_dataset(
            "double_autocorrelation", shape=(idx + 1,), maxshape=(None,), chunks=(1024,), dtype=np.complex128, **kwargs
        )
        g["double_autocorrelation"][...] = np.nan + 0.0j

    ds_c2 = g["double_autocorrelation"]
    if ds_c2.shape[0] < idx + 1:
        ds_c2.resize(idx + 1, axis=0)
    ds_c2[idx] = c2


def save_rothe_state(
    sim_name,
    splitting_type,
    polynomial_string,
    epsilon,
    dt,
    t,
    rothe_error,
    params,
    coeffs,
    rng_state=None,
    dipole_moment=None,
    double_autocorrelation=None,
    compression="gzip",
    compression_opts=4,
    path=None,
):
    """
    Append one step of the simulation to a single HDF5 file.

    Supports both real and imaginary time propagation. Complex `dt`/`t` values
    are stored losslessly (dataset ``t`` is upgraded to complex if needed).

    Data layout (dictionary-style):
      - attrs: sim_name, splitting_type, dt (complex), created, dtypes
      - dataset "polynomial_string": full input string
    - group "data": datasets "p", "c", "t", "rothe_error" (first dim = step index)
    - optional datasets: "dipole", "double_autocorrelation"
    """
    filename = _resolve_output_filename(sim_name, splitting_type, path=path)

    params_np = np.asarray(params)
    coeffs_np = np.asarray(coeffs)

    dt_scalar = np.asarray(dt, dtype=np.complex128).item()
    t_scalar = np.asarray(t, dtype=np.complex128).item()
    complex_time = bool(np.abs(np.imag(dt_scalar)) > 0 or np.abs(np.imag(t_scalar)) > 0)

    kwargs = _compression_kwargs(compression, compression_opts)

    with h5py.File(filename, "a") as f:
        g = f.require_group("data")

        _write_metadata(f, sim_name, splitting_type, dt_scalar, params_np, coeffs_np, epsilon)
        _ensure_polynomial_string_dataset(f, polynomial_string)

        # --- ensure datasets exist, growing tail shape if Gaussians were added ---
        params_shape = tuple(params_np.shape)
        coeffs_shape = tuple(coeffs_np.shape)

        ds_p = _maybe_grow_step_array_dataset(g, "p", params_shape, params_np.dtype, kwargs)
        ds_c = _maybe_grow_step_array_dataset(g, "c", coeffs_shape, coeffs_np.dtype, kwargs)
        ds_t = _ensure_time_dataset(g, complex_time, kwargs)
        ds_RE = _ensure_scalar_step_dataset(g, "rothe_error", np.float64, kwargs)
        ds_rng = _ensure_rng_state_dataset(g) if (rng_state is not None or "rng_state" in g) else None

        idx = ds_t.shape[0]
        ds_p.resize(idx + 1, axis=0)
        ds_c.resize(idx + 1, axis=0)
        ds_t.resize(idx + 1, axis=0)
        ds_RE.resize(idx + 1, axis=0)
        if ds_rng is not None:
            ds_rng.resize(idx + 1, axis=0)

        ds_p[idx, ...] = params_np
        ds_c[idx, ...] = coeffs_np
        ds_t[idx] = t_scalar if complex_time else float(np.real(t_scalar))
        ds_RE[idx] = float(rothe_error)
        if ds_rng is not None:
            ds_rng[idx] = _encode_rng_state(rng_state) if rng_state is not None else ""

        _append_optional_dipole(g, idx, dipole_moment, kwargs)
        _append_optional_double_autocorrelation(g, idx, double_autocorrelation, kwargs)

    return filename


def load_rothe_file(sim_name, splitting_type, path=None):
    """
    Read the entire HDF5 into a dict (numpy arrays).
    """
    filename = _resolve_output_filename(sim_name, splitting_type, path=path, create_dir=False)
    if not os.path.exists(filename):
        raise FileNotFoundError(filename)
    with h5py.File(filename, "r") as f:
        g = f["data"]
        poly = f["polynomial_string"][()]
        if isinstance(poly, bytes):
            poly = poly.decode("utf-8")
        params = g["p"][...]
        coeffs = g["c"][...]
        out = {
            "filename": filename,
            "attrs": dict(f.attrs),
            "polynomial_string": poly,
            "t": g["t"][...],
            "rothe_error": g["rothe_error"][...],
            "params": params,
            "coeffs": coeffs,
        }
        # Load dipole moment if it exists
        if "dipole" in g:
            out["dipole"] = g["dipole"][...]
        if "double_autocorrelation" in g:
            out["double_autocorrelation"] = g["double_autocorrelation"][...]
        if "rng_state" in g:
            out["rng_state"] = [_decode_rng_state(v) for v in g["rng_state"][...]]
    return out


def _select_and_trim_to_time(filename, t_target):
    """
    Open HDF5 in append mode, trim datasets in-place according to the rules:
      - If t_target exactly exists (within atol), keep up to and including it;
        delete everything after.
      - If t_target is beyond the latest saved time, raise an error.
      - If it doesn't exist but is within the saved range:
        - For real time: keep up to the last t < t_target.
        - For complex/imaginary time: keep up to the last step whose |t| < |t_target|.
    Returns dict with last-kept state (params, coeffs, t) to resume from.
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(filename)

    with h5py.File(filename, "a") as f:
        dt_attr = f.attrs.get("dt", 0.0)
        dt_mag = float(np.abs(dt_attr)) if dt_attr is not None else 0.0
        g = f["data"]
        ds_t, ds_p, ds_c, ds_RE = g["t"], g["p"], g["c"], g["rothe_error"]
        ds_dipole = g["dipole"] if "dipole" in g else None
        ds_c2 = g["double_autocorrelation"] if "double_autocorrelation" in g else None
        ds_rng = g["rng_state"] if "rng_state" in g else None
        t_arr = ds_t[...]
        if t_arr.size == 0:
            raise ValueError("File has no time steps.")

        if t_target is None:
            idx_keep = t_arr.size - 1  # resume at last available
            params_prev = ds_p[idx_keep - 1, ...] if idx_keep >= 1 else None
            rng_state = (
                _decode_rng_state(ds_rng[idx_keep]) if ds_rng is not None and ds_rng.shape[0] > idx_keep else None
            )
            return dict(
                t=t_arr[idx_keep],
                idx=int(idx_keep),
                params=ds_p[idx_keep, ...],
                coeffs=ds_c[idx_keep, ...],
                params_prev=params_prev,
                rng_state=rng_state,
                trimmed=False,
            )

        complex_time = np.iscomplexobj(t_arr) or np.iscomplexobj(t_target)
        atol = max(1e-12, 1e-6 * dt_mag if dt_mag > 0 else 1e-9)
        exact = np.where(np.isclose(t_arr, t_target, atol=atol, rtol=0.0))[0]
        if exact.size > 0:
            keep_len = int(exact[-1] + 1)  # include the exact match
        else:
            t_max = t_arr[-1]
            if complex_time:
                target_beyond = np.abs(t_target) > np.abs(t_max) + atol
            else:
                target_beyond = np.real(t_target) > np.real(t_max) + atol

            if target_beyond:
                raise ValueError(
                    f"Requested t={t_target} is beyond the latest saved time t={t_max} "
                    f"in {os.path.basename(filename)}. "
                    f"Use t_start=None to resume from the latest."
                )

            if complex_time:
                lt = np.where(np.abs(t_arr) < np.abs(t_target) - atol)[0]
            else:
                lt = np.where(t_arr < t_target - atol)[0]
            keep_len = int(lt[-1] + 1) if lt.size > 0 else 0

        trimmed = keep_len < t_arr.size
        if trimmed:
            ds_p.resize((keep_len,) + ds_p.shape[1:])
            ds_c.resize((keep_len,) + ds_c.shape[1:])
            ds_t.resize((keep_len,))
            ds_RE.resize((keep_len,))
            if ds_dipole is not None:
                ds_dipole.resize((keep_len,) + ds_dipole.shape[1:])
            if ds_c2 is not None:
                ds_c2.resize((keep_len,))
            if ds_rng is not None:
                ds_rng.resize((keep_len,))

        if keep_len == 0:
            raise ValueError(f"No saved state at or before t={t_target} in {os.path.basename(filename)}.")

        idx = keep_len - 1
        params_prev = ds_p[idx - 1, ...] if idx >= 1 else None
        rng_state = _decode_rng_state(ds_rng[idx]) if ds_rng is not None and ds_rng.shape[0] > idx else None
        return dict(
            t=ds_t[idx],
            idx=int(idx),
            params=ds_p[idx, ...],
            coeffs=ds_c[idx, ...],
            params_prev=params_prev,
            rng_state=rng_state,
            trimmed=trimmed,
        )
