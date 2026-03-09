"""
HDF5 I/O for the Rothe solver.

Each simulation is stored in a single HDF5 file with the layout::

    {sim_name}__{splitting_type}.h5
    ├── attrs: sim_name, splitting_type, dt (complex128), epsilon, created
    ├── polynomial_string   (scalar string dataset)
    └── steps/
        ├── 0000/
        │   ├── params              (n_total, D, 4)   float64
        │   ├── coeffs              (n_total,)        complex128
        │   ├── attrs: t, rothe_error, n_dynamic
        │   ├── rng_state           (optional, attr, JSON string)
        │   ├── dipole              (optional, dataset)
        │   └── double_autocorrelation  (optional, attr, complex128)
        ├── 0001/
        │   └── ...
        └── ...

Every time step is its own HDF5 group with its own arrays.  No padding,
no reshaping, no growing — each step stores exactly its own data at its
own shape.  Variable Gaussian counts across steps are natural.
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


# ── filename helpers ───────────────────────────────────────────────────────


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


# ── JSON encoding for RNG state ───────────────────────────────────────────


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


# ── step group name ────────────────────────────────────────────────────────

_STEP_FMT = "{:04d}"


def _step_group_name(idx):
    return _STEP_FMT.format(idx)


def _step_count(steps_group):
    """Return the number of step sub-groups already present."""
    return len(steps_group)


# ── save ───────────────────────────────────────────────────────────────────


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
    n_dynamic,
    rng_state=None,
    dipole_moment=None,
    double_autocorrelation=None,
    compression="gzip",
    compression_opts=4,
    path=None,
):
    """
    Append one time-step to the simulation HDF5 file.

    Each step is stored in its own group ``steps/NNNN`` with exact-shape
    datasets — no padding or resizing of earlier steps.
    """
    filename = _resolve_output_filename(sim_name, splitting_type, path=path)
    params_np = np.asarray(params)
    coeffs_np = np.asarray(coeffs)

    dt_scalar = np.asarray(dt, dtype=np.complex128).item()
    t_scalar = np.asarray(t, dtype=np.complex128).item()
    complex_time = bool(np.abs(np.imag(dt_scalar)) > 0 or np.abs(np.imag(t_scalar)) > 0)

    ckw = _compression_kwargs(compression, compression_opts)

    with h5py.File(filename, "a") as f:
        # ── file-level metadata (written once) ──
        f.attrs.setdefault("sim_name", sim_name)
        f.attrs.setdefault("splitting_type", str(splitting_type))
        f.attrs.setdefault("dt", np.complex128(dt_scalar))
        f.attrs.setdefault("epsilon", float(epsilon))
        f.attrs.setdefault("created", time.time())

        if "polynomial_string" not in f:
            ds = f.create_dataset("polynomial_string", shape=(), dtype=h5py.string_dtype(encoding="utf-8"))
            ds[()] = polynomial_string

        # ── step group ──
        steps = f.require_group("steps")
        idx = _step_count(steps)
        sg = steps.create_group(_step_group_name(idx))

        # Scalar attrs
        sg.attrs["t"] = t_scalar if complex_time else float(np.real(t_scalar))
        sg.attrs["rothe_error"] = float(rothe_error)
        sg.attrs["n_dynamic"] = int(n_dynamic)

        # Arrays
        sg.create_dataset("params", data=params_np, **ckw)
        sg.create_dataset("coeffs", data=coeffs_np, **ckw)

        # Optional fields
        if rng_state is not None:
            sg.attrs["rng_state"] = _encode_rng_state(rng_state)

        if dipole_moment is not None:
            sg.create_dataset("dipole", data=np.asarray(dipole_moment), **ckw)

        if double_autocorrelation is not None:
            sg.attrs["double_autocorrelation"] = np.complex128(double_autocorrelation)

    return filename


# ── load ───────────────────────────────────────────────────────────────────


def load_rothe_file(sim_name, splitting_type, path=None):
    """
    Read the entire HDF5 into a dict.

    Returns a dict with keys:
        filename, attrs, polynomial_string,
        t, rothe_error, n_dynamic            (1-D arrays, one entry per step)
        params, coeffs                       (lists of arrays — each element
                                              has its own shape)
        dipole, double_autocorrelation       (lists, entries may be None)
        rng_state                            (list, entries may be None)
    """
    filename = _resolve_output_filename(sim_name, splitting_type, path=path, create_dir=False)
    if not os.path.exists(filename):
        raise FileNotFoundError(filename)

    with h5py.File(filename, "r") as f:
        poly = f["polynomial_string"][()]
        if isinstance(poly, bytes):
            poly = poly.decode("utf-8")

        steps_grp = f["steps"]
        keys_sorted = sorted(steps_grp.keys())

        t_list = []
        re_list = []
        nd_list = []
        params_list = []
        coeffs_list = []
        dipole_list = []
        autocorr_list = []
        rng_list = []

        for k in keys_sorted:
            sg = steps_grp[k]
            t_list.append(sg.attrs["t"])
            re_list.append(sg.attrs["rothe_error"])
            nd_list.append(sg.attrs["n_dynamic"])
            params_list.append(sg["params"][...])
            coeffs_list.append(sg["coeffs"][...])
            dipole_list.append(sg["dipole"][...] if "dipole" in sg else None)
            autocorr_list.append(sg.attrs.get("double_autocorrelation", None))
            rng_enc = sg.attrs.get("rng_state", None)
            rng_list.append(_decode_rng_state(rng_enc) if rng_enc is not None else None)

        file_attrs = dict(f.attrs)

    out = {
        "filename": filename,
        "attrs": file_attrs,
        "polynomial_string": poly,
        "t": np.array(t_list),
        "rothe_error": np.array(re_list, dtype=np.float64),
        "n_dynamic": np.array(nd_list, dtype=np.int32),
        "params": params_list,
        "coeffs": coeffs_list,
        "dipole": dipole_list,
        "double_autocorrelation": autocorr_list,
        "rng_state": rng_list,
    }
    return out


# ── trim / resume ──────────────────────────────────────────────────────────


def _select_and_trim_to_time(filename, t_target):
    """
    Open HDF5 in append mode, locate the resume point and delete later
    step groups.

    Rules:
      - *t_target is None*: resume from the last step (no trimming).
      - *t_target exactly matches* a saved time (within atol): keep up to
        and including that step; delete everything after.
      - *t_target is beyond* the latest saved time: raise ``ValueError``.
      - *t_target is within range* but not exact: keep up to the last step
        whose time is strictly less than t_target.

    Returns a dict::

        t, idx, params, coeffs, n_dynamic,
        params_prev, n_dynamic_prev,
        rng_state, trimmed
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(filename)

    with h5py.File(filename, "a") as f:
        dt_attr = f.attrs.get("dt", 0.0)
        dt_mag = float(np.abs(dt_attr)) if dt_attr is not None else 0.0
        steps = f["steps"]
        keys_sorted = sorted(steps.keys())
        n_steps = len(keys_sorted)
        if n_steps == 0:
            raise ValueError("File has no time steps.")

        # Read all times
        t_arr = np.array([steps[k].attrs["t"] for k in keys_sorted])

        # ── determine keep_len ──
        if t_target is None:
            keep_len = n_steps
        else:
            complex_time = np.iscomplexobj(t_arr) or np.iscomplexobj(t_target)
            atol = max(1e-12, 1e-6 * dt_mag if dt_mag > 0 else 1e-9)

            exact = np.where(np.isclose(t_arr, t_target, atol=atol, rtol=0.0))[0]
            if exact.size > 0:
                keep_len = int(exact[-1] + 1)
            else:
                t_max = t_arr[-1]
                if complex_time:
                    beyond = np.abs(t_target) > np.abs(t_max) + atol
                else:
                    beyond = np.real(t_target) > np.real(t_max) + atol
                if beyond:
                    raise ValueError(
                        f"Requested t={t_target} is beyond the latest saved "
                        f"time t={t_max} in {os.path.basename(filename)}. "
                        f"Use t_start=None to resume from the latest."
                    )
                if complex_time:
                    lt = np.where(np.abs(t_arr) < np.abs(t_target) - atol)[0]
                else:
                    lt = np.where(t_arr < t_target - atol)[0]
                keep_len = int(lt[-1] + 1) if lt.size > 0 else 0

        if keep_len == 0:
            raise ValueError(f"No saved state at or before t={t_target} " f"in {os.path.basename(filename)}.")

        # ── trim ──
        trimmed = keep_len < n_steps
        if trimmed:
            for k in keys_sorted[keep_len:]:
                del steps[k]

        # ── read resume state ──
        idx = keep_len - 1
        sg = steps[keys_sorted[idx]]

        result = dict(
            t=sg.attrs["t"],
            idx=idx,
            params=sg["params"][...],
            coeffs=sg["coeffs"][...],
            n_dynamic=int(sg.attrs["n_dynamic"]),
            rng_state=_decode_rng_state(sg.attrs.get("rng_state", None)),
            trimmed=trimmed,
        )

        # Previous step (for extrapolation)
        if idx >= 1:
            sg_prev = steps[keys_sorted[idx - 1]]
            result["params_prev"] = sg_prev["params"][...]
            result["n_dynamic_prev"] = int(sg_prev.attrs["n_dynamic"])
        else:
            result["params_prev"] = None
            result["n_dynamic_prev"] = None

    return result
