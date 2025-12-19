import os
import time

os.environ["XLA_FLAGS"] = "--xla_gpu_enable_triton_gemm=true --xla_gpu_autotune_level=4"

import numpy as np
from scipy.optimize import minimize

import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)
from calculate_Hessian_coefficients import calculate_Hessian_matrix
from file_handling import _select_and_trim_to_time, save_rothe_state
from general_wf import ND_potentials, propagate_kinetic_analytical
from jax.scipy.linalg import solve


@jax.jit
def evaluate_function_1000_times(func, x):
    start = time.time()
    for _ in range(1000):
        func(x)
    end = time.time()
    time_taken = end - start
    print(f"Time taken for 1000 evaluations: {time_taken:.2f} seconds")
    return time_taken


def build_diagonal_scaling_from_grad(grad_p, epsilon=1e-12):
    """
    Build diagonal scaling S from p-space gradient.

    grad_p: np.ndarray, shape (n, D, 4)
        Gradient in p-space reshaped like p.
    Returns:
        S_flat:     np.ndarray, shape (n*D*4,)   (diagonal of S)
        S_flat_jax: jnp.ndarray, same data, for JAX ops
    """
    abs_grad = np.abs(grad_p)
    # RMS per parameter type (last axis: a, b, mu, p)
    group_rms = np.sqrt(np.mean(abs_grad**2, axis=(0, 1)))  # shape (4,)

    # Target RMS: average of the four
    desired_RMS = np.mean(group_rms)

    # Scale factors: s_k = avg / RMS_k, avoid division by zero
    scale_types = np.where(group_rms > epsilon, desired_RMS / group_rms, 1.0)  # shape (4,)

    # Broadcast to full shape (n, D, 4)
    scale_arr = np.einsum(
        "...k,k->...k", np.ones_like(grad_p), scale_types
    )  # last axis broadcasts
    S_flat = scale_arr.ravel()
    S_flat_jax = jnp.asarray(S_flat)
    return S_flat, S_flat_jax


def params_flat_to_y_flat(params_flat, S_flat):
    """
    Change of basis params = S y  ⇒  y = S^{-1} params.

    params_flat: np.ndarray, flattened params
    S_flat: np.ndarray, same shape, diagonal scaling
    """
    return params_flat / S_flat


def y_flat_to_params_flat(y_flat, S_flat):
    """
    Inverse change of basis y → params = S y.

    y_flat: np.ndarray, flattened y
    S_flat: np.ndarray, same shape, diagonal scaling
    """
    return y_flat * S_flat


def make_objective_params(
    rothe_error_and_gradient, params_shape, params_old, coeffs_old, SHH2, t, dt
):
    """
    Wrap rothe_error_and_gradient into a flat params-space objective.

    Returns:
        f_and_g_params(theta_flat) -> (float value, np.ndarray grad_flat)
        g_only_params(theta_flat)  -> grad_flat
    """

    def f_and_g_params(theta_flat):
        params_new = jnp.asarray(theta_flat).reshape(params_shape)
        val, g = rothe_error_and_gradient(params_new, params_old, coeffs_old, SHH2, t, dt)
        return float(val), np.asarray(jnp.ravel(g))

    def g_only_params(theta_flat):
        _, g = f_and_g_params(theta_flat)
        return g

    return f_and_g_params, g_only_params


def make_objective_y(
    rothe_error_and_gradient,
    params_shape,
    params_old,
    coeffs_old,
    SHH2,
    t,
    dt,
    S_flat_jax,
):
    """
    Wrap rothe_error_and_gradient into y-space objective, with params = S y.

    Returns:
        f_and_g_y(theta_y_flat) -> (float value, np.ndarray grad_y_flat)
    """

    def f_and_g_y(theta_y_flat):
        theta_y = jnp.asarray(theta_y_flat)
        params_flat = theta_y * S_flat_jax  # params = S * y
        params_new = params_flat.reshape(params_shape)
        val, grad_params = rothe_error_and_gradient(
            params_new, params_old, coeffs_old, SHH2, t, dt
        )
        grad_params_flat = jnp.ravel(grad_params)
        g_y_flat = grad_params_flat * S_flat_jax  # g_y = S * grad_params
        return float(val), np.asarray(g_y_flat)

    return f_and_g_y


def compute_gtol(initial_err):
    """
    Stopping criterion for BFGS.
    """
    return min(initial_err / 10.0, 1e-6)


def make_initial_hessian_y(dim_y, grad_y0_norm):
    """
    Simple scaled identity initial inverse Hessian in y-space.
    """
    return np.eye(dim_y) / grad_y0_norm


def maybe_apply_bb1_scaling(
    params_old,
    params_oldold,
    T,
    T_inv,
    grad_theta0_flat,
    g_only_params,
    t,
    hess_inv_default,
):

    if params_oldold is None:
        return hess_inv_default

    # Flatten p's
    params_old_flat = np.asarray(params_old).ravel()
    params_oldold_flat = np.asarray(params_oldold).ravel()

    # Step in p-space
    delta_params = params_old_flat - params_oldold_flat  # s_params

    # Step in theta-space: s_theta = T^{-1} delta_p
    s_theta = T_inv @ delta_params

    # Gradient at old-old point in p-space
    grad_oldold_params_flat = np.asarray(g_only_params(params_oldold_flat))

    # Gradient in theta at old-old: g_theta_oldold = T^T g_p_oldold
    g_theta_oldold = T.T @ grad_oldold_params_flat

    # Difference in theta-gradients
    grad_theta0_flat = np.asarray(grad_theta0_flat)
    y_vec = grad_theta0_flat - g_theta_oldold

    denom = float(s_theta @ y_vec)
    num = float(s_theta @ s_theta)

    if not np.isfinite(denom) or denom <= 0.0:
        print(f"{t} does not fulfill curvature condition: denom={denom}, using default.")
        return hess_inv_default

    BB1 = num / denom
    if not np.isfinite(BB1) or BB1 <= 0.0:
        print(f"{t} invalid BB1={BB1}, using default.")
        return hess_inv_default

    print(f"{t} fulfills curvature condition (theta-space): BB1={BB1}")
    dim = hess_inv_default.shape[0]
    return BB1 * np.eye(dim)


def get_approx_Hessian(params_init, coeffs):
    wf_obj = ND_potentials(params_init, 4)
    S = wf_obj.calculate_S()
    moments = wf_obj.moments
    hess_new = calculate_Hessian_matrix(params_init, S, moments, coeffs)
    n, D, _ = params_init.shape
    h_p_6d = jnp.transpose(hess_new, (0, 2, 4, 1, 3, 5))  # (n, D, 4, n, D, 4)
    i = np.arange(n)
    mask = i[:, None] == i[None, :]  # shape (n, n)

    # reshape to broadcast over (n, D, 4, n, D, 4)
    mask = mask[:, None, None, :, None, None]

    h_p_6d *= mask
    H_p = h_p_6d.reshape(n * D * 4, n * D * 4)
    # H_p = H_p + H_p.T  # Ensure Hermiticity

    return H_p


def hessian_cholesky(hessian):
    """Return a factor L such that H ≈ L^T L.

    Note:
        This currently clips eigenvalues to enforce positive definiteness.
        The function is conservative about numerical stability.

    Caveat:
        The dense-metric path in `RotheSolver.find_next_timestep_solution` is
        experimental; if you rely on this for preconditioning, confirm the
        scaling is actually enabled.
    """
    hessian_np = np.asarray(hessian, dtype=np.float64)
    w, Q = np.linalg.eigh(hessian_np)
    # Clip eigenvalues to avoid insane conditioning / negatives
    lam_med = np.median(w[w > 0])  # median of positive ones
    lam_min = lam_med * 1e-3
    w_clipped = np.clip(w, lam_min, None)
    # w_clipped = np.ones_like(w)  # TEMPORARY: use identity
    # Build L such that H_p ≈ L^T L
    sqrt_w = np.sqrt(w_clipped)
    L = sqrt_w[None, :] * Q.T
    return np.eye(L.shape[0])  # TEMPORARY: use identity


def make_objective_z(
    rothe_error_and_gradient,
    params_shape,
    params_old,
    coeffs_old,
    SHH2,
    t,
    dt,
    L_inv,
):
    """
    z-space objective: params = L^{-1} z, g_z = L^{-T} grad_params
    """
    L_inv_T = L_inv.T

    def f_and_g_z(z_flat):
        z = jnp.asarray(z_flat)
        # p = L^{-1} z
        params_flat = L_inv @ np.asarray(z)  # use np matmul, that's fine here
        params_new = params_flat.reshape(params_shape)
        val, grad_params = rothe_error_and_gradient(
            params_new, params_old, coeffs_old, SHH2, t, dt
        )
        grad_params_flat = np.asarray(jnp.ravel(grad_params))
        # g_z = L^{-T} g_p
        g_z_flat = L_inv_T @ grad_params_flat
        return float(val), g_z_flat

    return f_and_g_z


class RotheSolver:
    """
    Implicit Rothe time stepper for the TDSE in a non-orthogonal Gaussian basis.

    You provide:
        - `SHH2(t, params, splitting_type) -> (S, H, H2)` where
            S, H, H2 are (N,N) matrices in the Gaussian basis.
        - `rothe_grad_jax`: a JAX function returning (value, grad) w.r.t. p_new.
        - `rothe_nograd`: a non-grad version (used to compute c_new).

    State:
        - `params_old`: Gaussian parameters (typically shape (n, D, 4) with (a,b,mu,p)).
        - `coeffs_old`: coefficients (shape (n,), complex).

    Output:
        - Optionally appends to an HDF5 log via `file_handling.save_rothe_state`.
    """

    def __init__(
        self,
        SHH2,
        dt,
        t,
        params_old,
        coeffs_old,
        rothe_grad_jax,
        rothe_nograd=None,
        splitting_type="none",
        name=None,
        polynomial_string=None,
        out_dir="./wave_function_data",
        compression="gzip",
        compression_opts=4,
    ):
        self.SHH2 = SHH2
        self.dt = dt
        self.t = t
        self.params_old = jnp.asarray(params_old)
        self.coeffs_old = jnp.asarray(coeffs_old)
        self.rothe_error_and_gradient = rothe_grad_jax
        self.rothe_nograd = rothe_nograd
        self.splitting_type = splitting_type
        self.total_time = 0
        self.params_oldold = None
        self.params_oldoldold = None
        self.name = name
        self.polynomial_string = polynomial_string
        self.out_dir = out_dir
        self.compression = compression
        self.compression_opts = compression_opts

    def find_next_timestep_solution(self, params_init=None, maxiter=100):
        # 1) Choose initial p and store shape
        if params_init is None:
            params_init = np.asarray(self.params_old)
        self.params_init = np.asarray(params_init)

        params_shape = self.params_init.shape
        x0_params_flat = self.params_init.ravel()  # initial point in params-space
        N = x0_params_flat.size

        # 2) p-space objective + gradient (needed in any case)
        f_and_g_params, g_only_params = make_objective_params(
            self.rothe_error_and_gradient,
            params_shape,
            self.params_old,
            self.coeffs_old,
            self.SHH2,
            self.t,
            self.dt,
        )

        # 3) Evaluate at initial guess in p-space
        initial_err, initial_grad_params_flat = f_and_g_params(x0_params_flat)
        print("Initial squared Rothe error: ", initial_err)
        gtol = compute_gtol(initial_err)

        # ------------------------------------------------------------------
        # CHOOSE COORDINATE SYSTEM:
        #   use_dense_metric = False  -> diagonal scaling
        #   use_dense_metric = True   -> pseudo_Hessian-scaling
        # ------------------------------------------------------------------
        # The dense-metric path is currently a placeholder (see `hessian_cholesky`),
        # and computing the approximate Hessian is expensive. Default to the
        # diagonal/identity path unless you explicitly want to experiment.
        use_dense_metric = False

        if not use_dense_metric:
            # === Diagonal scaling S: p = S y =================================
            grad_params = initial_grad_params_flat.reshape(params_shape)
            S_flat, S_flat_jax = build_diagonal_scaling_from_grad(grad_params)
            # If you actually want scaling, comment out these two lines:
            S_flat_jax = jnp.ones_like(S_flat_jax)
            S_flat = np.ones_like(S_flat)

            # y0 and gradient in y-space
            theta0 = params_flat_to_y_flat(x0_params_flat, S_flat)
            grad_theta0_flat = initial_grad_params_flat * S_flat  # g_y = S g_params
            grad_theta0_norm = np.linalg.norm(grad_theta0_flat)

            # linear map p = T theta with theta=y
            T = np.diag(S_flat)  # p = S y
            T_inv = np.diag(1.0 / S_flat)  # y = S^{-1} p

            # y-space objective
            f_and_g_theta = make_objective_y(
                self.rothe_error_and_gradient,
                params_shape,
                self.params_old,
                self.coeffs_old,
                self.SHH2,
                self.t,
                self.dt,
                S_flat_jax,
            )

        else:
            # === Dense metric from approximate Hessian: p = L^{-1} z =========
            hessian = get_approx_Hessian(self.params_init, self.coeffs_old)
            L = hessian_cholesky(hessian)  # N x N
            L_inv = np.linalg.inv(L)

            # z0 = L p0
            theta0 = L @ x0_params_flat  # this is z0

            # z-space objective: p = L^{-1} z, g_z = L^{-T} g_p
            f_and_g_theta = make_objective_z(
                self.rothe_error_and_gradient,
                params_shape,
                self.params_old,
                self.coeffs_old,
                self.SHH2,
                self.t,
                self.dt,
                L_inv,
            )
            # gradient in z at z0
            _, grad_theta0_flat = f_and_g_theta(theta0)
            grad_theta0_norm = np.linalg.norm(grad_theta0_flat)

            # linear map p = T theta with theta=z
            T = L_inv  # p = L^{-1} z
            T_inv = L  # z = L p

        # 4) Initial inverse Hessian (possibly BB1-scaled) in theta-space
        """
        hess_inv_default = make_initial_hessian_y(N, grad_theta0_norm)

        hess_inv0 = maybe_apply_bb1_scaling(
            params_old=self.params_old,
            params_oldold=getattr(self, "params_oldold", None),
            T=T,
            T_inv=T_inv,
            grad_theta0_flat=grad_theta0_flat,
            g_only_params=g_only_params,
            t=self.t,
            hess_inv_default=hess_inv_default,
        )
        """
        hess_inv0 = np.zeros((N, N), dtype=np.float64)
        for i in range(N):
            hess_inv0[i, i] = 1.0
        hess_inv0 = 1e-8 * hess_inv0
        # 5) Run BFGS in theta-space
        solution = minimize(
            f_and_g_theta,
            theta0,
            jac=True,
            method="BFGS",
            options={
                "gtol": gtol,
                "maxiter": maxiter,
                "disp": False,
                "hess_inv0": hess_inv0,
            },
        )
        niter = solution.nfev
        theta_solved = solution.x

        # 6) Transform solution back to p-space
        if not use_dense_metric:
            # theta = y, p = S y
            params_solved_flat = y_flat_to_params_flat(theta_solved, S_flat)
        else:
            # theta = z, p = L^{-1} z
            params_solved_flat = L_inv @ theta_solved

        params_solved = params_solved_flat.reshape(params_shape)

        # 7) Final Rothe error and new coefficients
        final_RE, coeffs_new = self.rothe_nograd(
            jnp.asarray(params_solved),
            self.params_old,
            self.coeffs_old,
            self.SHH2,
            self.t,
            self.dt,
            return_cnew=True,
        )
        return params_solved, coeffs_new, float(final_RE), niter

    def evaluate_gradient_n_times(self, params_init, coeffs_init, n=100):
        if params_init is not None:
            self.params_old = jnp.asarray(params_init)
        if coeffs_init is not None:
            self.coeffs_old = jnp.asarray(coeffs_init)
        if params_init is None:
            params_init = np.asarray(self.params_old)
        self.params_init = np.asarray(params_init)
        self._p_shape = self.params_old.shape
        self._p_size = int(np.prod(self._p_shape))

        # Jitted “flat” objective+grad: theta_flat -> (val, grad_flat)
        def rothe_fg_flat(theta_flat, params_old, coeffs_old, t, dt):
            params_new = theta_flat.reshape(self._p_shape)
            val, g = self.rothe_error_and_gradient(
                params_new, params_old, coeffs_old, self.SHH2, t, dt
            )
            return val, g.ravel()

        # JIT once; SHH2 is closed over, shapes are static
        self._rothe_fg_flat = jax.jit(rothe_fg_flat)

        def f_and_g(theta_flat_np):
            theta_flat = jnp.asarray(theta_flat_np)
            val, g_flat = self._rothe_fg_flat(
                theta_flat, self.params_old, self.coeffs_old, self.t, self.dt
            )
            # SciPy wants host-side types
            return float(val), np.asarray(g_flat, dtype=np.float64)

        x0 = self.params_init.ravel()
        # hess_inv=np.eye(len(x0)) #Initial Hessian
        val, g = f_and_g(x0)  # warmup for jax compilation
        start = time.time()
        for i in range(n):
            val, g = f_and_g(x0)
        end = time.time()
        time_taken = end - start
        print(f"Time taken for {n} gradient evaluations: {time_taken:.2f} seconds")
        return time_taken

    def propagate(
        self, num_iterations, params_init=None, coeffs_init=None, maxiter=100, renormalize=True
    ):
        if params_init is not None:
            self.params_old = jnp.asarray(params_init)
        if coeffs_init is not None:
            self.coeffs_old = jnp.asarray(coeffs_init)
        startguess = self.params_old
        D = int(self.params_old.shape[1])
        onesD = jnp.ones(D)
        for i in range(num_iterations):
            self.t += self.dt
            dt = self.dt
            start = time.time()
            if self.splitting_type in (None, "none"):
                params_solved, coeffs_new, final_RE, nit = self.find_next_timestep_solution(
                    startguess, maxiter=maxiter
                )
            elif self.splitting_type == "kinetic":
                # Half step with T

                params_old_temp = self.params_old.copy()
                coeffs_old_temp = self.coeffs_old.copy()
                self.params_old, self.coeffs_old = propagate_kinetic_analytical(
                    dt / 2, self.params_old, self.coeffs_old, onesD
                )
                # Need to update coeffs/params_old as find_next_timestep_solution uses it as its base point
                startguess, _ = propagate_kinetic_analytical(
                    -dt / 2, startguess, self.coeffs_old, onesD
                )
                params_solved_temp, coeffs_new_temp, final_RE, nit = (
                    self.find_next_timestep_solution(startguess, maxiter=maxiter)
                )
                params_solved, coeffs_new = propagate_kinetic_analytical(
                    dt / 2, params_solved_temp, coeffs_new_temp, onesD
                )
                self.params_old = params_old_temp  # Restore to original values
                self.coeffs_old = coeffs_old_temp
            else:
                raise ValueError(
                    f"Unknown splitting_type {self.splitting_type!r}; expected 'none' or 'kinetic'."
                )
            if np.isrealobj(self.t) and not np.iscomplexobj(self.t):
                time_str = f"{float(np.real(self.t)):.2f}"
            else:
                time_str = str(self.t)
            print(
                f"Time {time_str}:Squared Rothe error: {final_RE:.3e}, Number of iterations: {nit}"
            )
            if renormalize:
                S, H, _ = self.SHH2(
                    t=self.t, params=params_solved, splitting_type=self.splitting_type
                )
                norm = jnp.real(jnp.conj(coeffs_new) @ S @ coeffs_new)
                print("Norm before renormalization: ", norm)
                coeffs_new = coeffs_new / norm
            if self.name is not None:
                save_rothe_state(
                    self.name,
                    self.splitting_type,
                    self.polynomial_string or "",
                    self.dt,
                    self.t,
                    final_RE,
                    params_solved,
                    coeffs_new,
                    compression=self.compression,
                    compression_opts=self.compression_opts,
                    path=self.out_dir,
                )
            if self.params_old is None:
                startguess = params_solved
            startguess = 2 * params_solved - self.params_old

            # startguess = p_solved

            self.params_oldold = self.params_old
            self.params_old = params_solved

            self.coeffs_old = coeffs_new
            end = time.time()
            time_taken = end - start
            self.total_time += time_taken

            print(f"Time taken: {time_taken:.2f}")
            # print(p_solved)
            # print(c_new)
        print(f"Total time: {self.total_time}")

    def resume_from_file(self, t_start=None):
        """
        Load state from HDF5 and set solver to resume at requested time.
        If t_start is None -> resume from last saved step.
        If t_start exists -> resume at that step and delete later steps.
        If t_start doesn't exist -> resume from latest step < t_start and delete all >= t_start.
        """
        if self.name is None:
            raise ValueError("Cannot resume: 'name' is None (no filename to load).")
        filename = os.path.join(self.out_dir, f"{self.name}__{self.splitting_type}.h5")
        sel = _select_and_trim_to_time(filename, t_start)

        # optional safety: check shapes
        if sel["params"].shape != tuple(self.params_old.shape):
            raise ValueError(
                f"Saved params shape {sel['params'].shape} != current {tuple(self.params_old.shape)}"
            )
        if sel["coeffs"].shape != tuple(self.coeffs_old.shape):
            raise ValueError(
                f"Saved coeffs shape {sel['coeffs'].shape} != current {tuple(self.coeffs_old.shape)}"
            )

        self.params_old = jnp.asarray(sel["params"])
        self.coeffs_old = jnp.asarray(sel["coeffs"])
        self.t = np.asarray(sel["t"]).item()
        self.params_oldold = None  # reset history across resume
        return sel  # contains t, idx, trimmed flag


def setUpRotheErrorAndGradient_jit(splitting_type):
    """Build the Rothe objective and its JAX gradient.

    Args:
        splitting_type: "none" uses (H, H2); otherwise uses (V, V2) in SHH2.

    Returns:
        (rothe_error, rothe_vg_jit)

                - `rothe_error(params_new, params_old, coeffs_old, SHH2, t, dt, return_cnew=False)`
          returns a real scalar (squared Rothe error) and optionally `c_new`.

        - `rothe_vg_jit` is `jax.value_and_grad(rothe_error, argnums=0)` jitted.

    Notes:
        `SHH2` is passed as a *static* arg to JIT; it should be a stable callable.
    """

    def rothe_error(
        params_new, params_old, coeffs_old, SHH2, t, dt, return_cnew=False, lambda_=1e-10
    ):
        splitting_eff = "none" if splitting_type in (None, "none") else splitting_type
        ngo = params_old.shape[0]
        dtsq4 = dt**2 / 4
        params_concat = jnp.concatenate(
            (params_old, params_new), axis=0
        )  # New shape: (2*n,N_b)
        S_full, H_full, H2_full = SHH2(
            t=t + dt / 2, params=params_concat, splitting_type=splitting_eff
        )
        S_tilde_full = S_full + dtsq4 * H2_full
        rho_mat = (
            S_full[:ngo, ngo:] - dtsq4 * H2_full[:ngo, ngo:] + 1j * dt * H_full[:ngo, ngo:]
        )
        rho_vec = jnp.conj(rho_mat).T @ coeffs_old
        # S_tilde_full = _symh(S_tilde_full)  # Ensure Hermiticity
        S_reg = (
            S_tilde_full + jnp.eye(S_tilde_full.shape[0]) * lambda_
        )  # Regularization term to avoid singular matrices.
        # L = jnp.linalg.cholesky(S_reg[ngo:, ngo:])  # Hermitian PD
        # c_new = jsp.linalg.solve_triangular(L, rho_vec, lower=True, trans="N")
        # c_new = jsp.linalg.solve_triangular(L.conj().T, c_new, lower=False, trans="N")
        # if jnp.abs(dt-jnp.real(dt))>1e-12:
        #    coeffs_new = solve(S_reg[ngo:, ngo:], rho_vec) #Hermiticity not assumed for complex time steps
        # else:
        coeffs_new = solve(S_reg[ngo:, ngo:], rho_vec)
        overlap_term = jnp.conj(coeffs_old) @ S_tilde_full[:ngo, :ngo] @ coeffs_old
        projection_term = jnp.real(jnp.conj(rho_vec).T @ coeffs_new)
        rothe_error = jnp.real(overlap_term - projection_term)
        if return_cnew:
            return rothe_error, coeffs_new
        else:
            return rothe_error

    rothe_vg_jit = jax.jit(
        jax.value_and_grad(rothe_error, argnums=0), static_argnames=("SHH2",)
    )
    return rothe_error, rothe_vg_jit
