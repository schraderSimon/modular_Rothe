import os
import time

import h5py
import numpy as np
from attr import dataclass


@dataclass
class OutputConfig:
    """Configuration for saving solver output."""

    name: None | str = None
    polynomial_string: None | str = None
    out_dir: str = "./wave_function_data"
    compression: str = "gzip"
    compression_opts: int = 4


def save_rothe_state(
    sim_name,
    splitting_type,
    polynomial_string,
    dt,
    t,
    rothe_error,
    params,
    coeffs,
    dipole_moment=None,
    compression="gzip",
    compression_opts=4,
    path="./wave_function_data",
):
    """
    Append one step of the simulation to a single HDF5 file.

    Supports both real and imaginary time propagation. Complex `dt`/`t` values
    are stored losslessly (dataset ``t`` is upgraded to complex if needed).

    Data layout (dictionary-style):
      - attrs: sim_name, splitting_type, dt (complex), created, dtypes
      - dataset "polynomial_string": full input string
      - group "data": datasets "p", "c", "t", "rothe_error" (first dim = step index)
    """
    filename = os.path.join(path, f"{sim_name}__{splitting_type}.h5")
    os.makedirs(path, exist_ok=True)

    params_np = np.asarray(params)
    coeffs_np = np.asarray(coeffs)

    dt_scalar = np.asarray(dt, dtype=np.complex128).item()
    t_scalar = np.asarray(t, dtype=np.complex128).item()
    complex_time = bool(np.abs(np.imag(dt_scalar)) > 0 or np.abs(np.imag(t_scalar)) > 0)

    def _kwargs():
        return dict(
            compression=compression if compression else None,
            compression_opts=compression_opts if compression else None,
            shuffle=bool(compression),
        )

    with h5py.File(filename, "a") as f:
        g = f.require_group("data")

        # --- metadata (idempotent) ---
        f.attrs.setdefault("sim_name", sim_name)
        f.attrs.setdefault("splitting_type", str(splitting_type))
        # store dt as complex to keep imaginary time information
        f.attrs.setdefault("dt", np.complex128(dt_scalar))
        f.attrs.setdefault("created", time.time())
        # Keep attribute names stable for older files.
        f.attrs["p_dtype"] = str(params_np.dtype)
        f.attrs["c_dtype"] = str(coeffs_np.dtype)

        # store full polynomial string as scalar vlen UTF-8
        dt_str = h5py.string_dtype(encoding="utf-8")
        if "polynomial_string" not in f:
            dset = f.create_dataset("polynomial_string", shape=(), dtype=dt_str)
            dset[()] = polynomial_string

        # --- ensure datasets exist (robust to half-written files) ---
        params_shape = tuple(params_np.shape)
        coeffs_shape = tuple(coeffs_np.shape)

        if "p" not in g:
            g.create_dataset(
                "p",
                shape=(0,) + params_shape,
                maxshape=(None,) + params_shape,
                chunks=(1,) + params_shape,
                dtype=params_np.dtype,
                **_kwargs(),
            )
        if "c" not in g:
            g.create_dataset(
                "c",
                shape=(0,) + coeffs_shape,
                maxshape=(None,) + coeffs_shape,
                chunks=(1,) + coeffs_shape,
                dtype=coeffs_np.dtype,
                **_kwargs(),
            )

        # Time dataset: upgrade to complex if imaginary time is used.
        def _ensure_time_dataset():
            target_dtype = np.complex128 if complex_time else np.float64
            if "t" not in g:
                return g.create_dataset(
                    "t",
                    shape=(0,),
                    maxshape=(None,),
                    chunks=(1024,),
                    dtype=target_dtype,
                    **_kwargs(),
                )

            ds_existing = g["t"]
            if ds_existing.dtype == target_dtype:
                return ds_existing

            # Upgrade from real -> complex while preserving existing data
            data_existing = ds_existing[...]
            del g["t"]
            return g.create_dataset(
                "t",
                data=data_existing.astype(target_dtype),
                maxshape=(None,),
                chunks=(1024,),
                dtype=target_dtype,
                **_kwargs(),
            )

        if "rothe_error" not in g:
            g.create_dataset(
                "rothe_error",
                shape=(0,),
                maxshape=(None,),
                chunks=(1024,),
                dtype=np.float64,
                **_kwargs(),
            )

        ds_p, ds_c, ds_t, ds_RE = g["p"], g["c"], _ensure_time_dataset(), g["rothe_error"]

        # sanity-check tail shapes
        if ds_p.shape[1:] != params_shape:
            raise ValueError(f"'p' tail shape {ds_p.shape[1:]} != current {params_shape}")
        if ds_c.shape[1:] != coeffs_shape:
            raise ValueError(f"'c' tail shape {ds_c.shape[1:]} != current {coeffs_shape}")

        idx = ds_t.shape[0]
        ds_p.resize(idx + 1, axis=0)
        ds_c.resize(idx + 1, axis=0)
        ds_t.resize(idx + 1, axis=0)
        ds_RE.resize(idx + 1, axis=0)

        ds_p[idx, ...] = params_np
        ds_c[idx, ...] = coeffs_np
        ds_t[idx] = t_scalar if complex_time else float(np.real(t_scalar))
        ds_RE[idx] = float(rothe_error)

        # Store dipole moment if provided
        if dipole_moment is not None:
            dipole_np = np.asarray(dipole_moment)
            dipole_shape = tuple(dipole_np.shape)

            if "dipole" not in g:
                g.create_dataset(
                    "dipole",
                    shape=(0,) + dipole_shape,
                    maxshape=(None,) + dipole_shape,
                    chunks=(1,) + dipole_shape,
                    dtype=dipole_np.dtype,
                    **_kwargs(),
                )

            ds_dipole = g["dipole"]
            ds_dipole.resize(idx + 1, axis=0)
            ds_dipole[idx, ...] = dipole_np

    return filename


def load_rothe_file(sim_name, splitting_type, path="./wave_function_data"):
    """
    Read the entire HDF5 into a dict (numpy arrays).
    """
    filename = os.path.join(path, f"{sim_name}__{splitting_type}.h5")
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
    return out


def _select_and_trim_to_time(filename, t_target):
    """
    Open HDF5 in append mode, trim datasets in-place according to the rules:
      - If t_target exactly exists (within atol), keep up to and including it; delete everything after.
      - If t_target is beyond the latest saved time, raise an error (no data available that far).
      - If it doesn't exist but is within the saved range:
        - For real time: keep up to the last t < t_target; delete everything from t >= t_target.
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
        t_arr = ds_t[...]
        if t_arr.size == 0:
            raise ValueError("File has no time steps.")

        if t_target is None:
            idx_keep = t_arr.size - 1  # resume at last available
            params_prev = ds_p[idx_keep - 1, ...] if idx_keep >= 1 else None
            return dict(
                t=t_arr[idx_keep],
                idx=int(idx_keep),
                params=ds_p[idx_keep, ...],
                coeffs=ds_c[idx_keep, ...],
                params_prev=params_prev,
                trimmed=False,
            )

        complex_time = np.iscomplexobj(t_arr) or np.iscomplexobj(t_target)
        atol = max(1e-12, 1e-6 * dt_mag if dt_mag > 0 else 1e-9)
        exact = np.where(np.isclose(t_arr, t_target, atol=atol, rtol=0.0))[0]
        if exact.size > 0:
            keep_len = int(exact[-1] + 1)  # include the exact match
        else:
            # Check if t_target is beyond the available data
            t_max = t_arr[-1]  # latest saved time
            if complex_time:
                # For imaginary time propagation, "beyond" means larger magnitude
                target_beyond = np.abs(t_target) > np.abs(t_max) + atol
            else:
                # For real time, "beyond" means larger value
                target_beyond = np.real(t_target) > np.real(t_max) + atol

            if target_beyond:
                raise ValueError(
                    f"Requested t={t_target} is beyond the latest saved time t={t_max} "
                    f"in {os.path.basename(filename)}. Use t_start=None to resume from the latest."
                )

            # t_target is within the saved range but not an exact match - trim to before it
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

        if keep_len == 0:
            raise ValueError(
                f"No saved state at or before t={t_target} in {os.path.basename(filename)}."
            )

        idx = keep_len - 1
        params_prev = ds_p[idx - 1, ...] if idx >= 1 else None
        return dict(
            t=ds_t[idx],
            idx=int(idx),
            params=ds_p[idx, ...],
            coeffs=ds_c[idx, ...],
            params_prev=params_prev,
            trimmed=trimmed,
        )
