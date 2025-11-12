import os
import time

import h5py
import numpy as np


def save_rothe_state(
    sim_name,
    splitting_type,
    polynomial_string,
    dt,
    t,
    rothe_error,
    p,
    c,
    compression="gzip",
    compression_opts=4,
    path=".",
):
    """
    Append one step of the simulation to a single HDF5 file.
    Data layout (dictionary-style):
      - attrs: sim_name, splitting_type, dt, created, dtypes
      - dataset "polynomial_string": full input string
      - group "data": datasets "p", "c", "t", "rothe_error" (first dim = step index)
    """
    filename = os.path.join(path, f"{sim_name}__{splitting_type}.h5")
    p_np = np.asarray(p)
    c_np = np.asarray(c)

    with h5py.File(filename, "a") as f:
        g = f.require_group("data")

        # --- metadata (idempotent) ---
        f.attrs.setdefault("sim_name", sim_name)
        f.attrs.setdefault("splitting_type", str(splitting_type))
        f.attrs.setdefault("dt", float(dt))
        f.attrs.setdefault("created", time.time())
        f.attrs["p_dtype"] = str(p_np.dtype)
        f.attrs["c_dtype"] = str(c_np.dtype)

        # store full polynomial string as scalar vlen UTF-8
        dt_str = h5py.string_dtype(encoding="utf-8")
        if "polynomial_string" not in f:
            dset = f.create_dataset("polynomial_string", shape=(), dtype=dt_str)
            dset[()] = polynomial_string

        # --- ensure datasets exist (robust to half-written files) ---
        p_shape = tuple(p_np.shape)
        c_shape = tuple(c_np.shape)

        def _kwargs():
            return dict(
                compression=compression if compression else None,
                compression_opts=compression_opts if compression else None,
                shuffle=bool(compression),
            )

        if "p" not in g:
            g.create_dataset(
                "p",
                shape=(0,) + p_shape,
                maxshape=(None,) + p_shape,
                chunks=(1,) + p_shape,
                dtype=p_np.dtype,
                **_kwargs(),
            )
        if "c" not in g:
            g.create_dataset(
                "c",
                shape=(0,) + c_shape,
                maxshape=(None,) + c_shape,
                chunks=(1,) + c_shape,
                dtype=c_np.dtype,
                **_kwargs(),
            )
        if "t" not in g:
            g.create_dataset(
                "t",
                shape=(0,),
                maxshape=(None,),
                chunks=(1024,),
                dtype=np.float64,
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

        ds_p, ds_c, ds_t, ds_RE = g["p"], g["c"], g["t"], g["rothe_error"]

        # sanity-check tail shapes
        if ds_p.shape[1:] != p_shape:
            raise ValueError(f"'p' tail shape {ds_p.shape[1:]} != current {p_shape}")
        if ds_c.shape[1:] != c_shape:
            raise ValueError(f"'c' tail shape {ds_c.shape[1:]} != current {c_shape}")

        idx = ds_t.shape[0]
        ds_p.resize(idx + 1, axis=0)
        ds_c.resize(idx + 1, axis=0)
        ds_t.resize(idx + 1, axis=0)
        ds_RE.resize(idx + 1, axis=0)

        ds_p[idx, ...] = p_np
        ds_c[idx, ...] = c_np
        ds_t[idx] = float(t)
        ds_RE[idx] = float(rothe_error)
    return filename


def load_rothe_file(sim_name, splitting_type, path="."):
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
        out = {
            "filename": filename,
            "attrs": dict(f.attrs),
            "polynomial_string": poly,
            "t": g["t"][...],
            "rothe_error": g["rothe_error"][...],
            "p": g["p"][...],
            "c": g["c"][...],
        }
    return out


def _select_and_trim_to_time(filename, t_target):
    """
    Open HDF5 in append mode, trim datasets in-place according to the rules:
      - If t_target exactly exists (within atol), keep up to and including it; delete everything after.
      - If it doesn't exist, keep up to the last t < t_target; delete everything from t >= t_target.
    Returns dict with last-kept state (p,c,t) to resume from.
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(filename)

    with h5py.File(filename, "a") as f:
        dt = float(f.attrs.get("dt", 0.0))
        g = f["data"]
        ds_t, ds_p, ds_c, ds_RE = g["t"], g["p"], g["c"], g["rothe_error"]
        t_arr = ds_t[...]
        if t_arr.size == 0:
            raise ValueError("File has no time steps.")

        if t_target is None:
            idx_keep = t_arr.size - 1  # resume at last available
            return dict(
                t=float(t_arr[idx_keep]),
                idx=int(idx_keep),
                p=ds_p[idx_keep, ...],
                c=ds_c[idx_keep, ...],
                trimmed=False,
            )

        atol = max(1e-12, 1e-6 * dt if dt > 0 else 1e-9)
        exact = np.where(np.isclose(t_arr, t_target, atol=atol, rtol=0.0))[0]
        if exact.size > 0:
            keep_len = int(exact[-1] + 1)  # include the exact match
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
        return dict(
            t=float(ds_t[idx]),
            idx=int(idx),
            p=ds_p[idx, ...],
            c=ds_c[idx, ...],
            trimmed=trimmed,
        )
