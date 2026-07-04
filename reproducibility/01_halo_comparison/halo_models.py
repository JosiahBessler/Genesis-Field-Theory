
import numpy as np
from dataclasses import dataclass
from typing import Callable, Tuple

@dataclass
class HaloModel:
    name: str
    # vfunc(r, params) -> velocity array
    vfunc: Callable[[np.ndarray, Tuple[float, float]], np.ndarray]
    n_params: int = 2

# ---------------- GenFT Template C ----------------
# Uses your exported GenFT halo template:
#   gft_modelA/gft_halo_template_C.npz
#
# That file must contain:
#   xi   : dimensionless radius grid
#   f_xi : dimensionless halo shape f_C(xi), normalized so max(f_C) = 1

import os
import numpy as np
from scipy.interpolate import interp1d

from dataclasses import dataclass
from typing import Callable, Tuple

@dataclass
class HaloModel:
    name: str
    # vfunc(r, params) -> velocity array
    vfunc: Callable[[np.ndarray, Tuple[float, float]], np.ndarray]
    n_params: int = 2

# Path to Template C file (relative to this script)
TEMPLATE_PATH = os.path.join(
    os.path.dirname(__file__),
    "gft_modelA",
    "gft_halo_template_C.npz"
)

try:
    _templateC = np.load(TEMPLATE_PATH)
    _xi_grid = _templateC["xi"]      # dimensionless radii
    _fC_grid = _templateC["f_xi"]    # dimensionless velocities (already normalized)
except Exception as e:
    raise FileNotFoundError(f"Could not load Template C at {TEMPLATE_PATH}: {e}")

# Interpolator for f_C(ξ)
_fC_interp = interp1d(
    _xi_grid,
    _fC_grid,
    kind="cubic",
    fill_value="extrapolate"
)

def f_genft_dimensionless(x: np.ndarray) -> np.ndarray:
    """
    GenFT Template C dimensionless halo shape f_C(ξ),
    where ξ = r / r0 and max f_C = 1.
    """
    x = np.asarray(x, dtype=float)
    return _fC_interp(x)

def v_genft(r: np.ndarray, params):
    """
    v_GenFT(r) = V0 * f_C(r / r0)
    with two free parameters per galaxy:
      r0 : radial scale
      V0 : velocity scale
    """
    r0, V0 = params
    x = r / r0
    return V0 * f_genft_dimensionless(x)


# ---------------- NFW ----------------

def g_nfw(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x = np.where(x == 0.0, 1e-8, x)
    num = np.log(1.0 + x) - x / (1.0 + x)
    den = x
    g2 = num / den
    g2 = np.maximum(g2, 0.0)
    return np.sqrt(g2)

def v_nfw(r: np.ndarray, params):
    r_s, V0 = params
    x = r / r_s
    return V0 * g_nfw(x)

# ---------------- Pseudo-isothermal sphere ----------------

def g_pseudo_iso(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x = np.where(x == 0.0, 1e-8, x)
    num = x - np.arctan(x)
    den = x
    g2 = num / den
    g2 = np.maximum(g2, 0.0)
    return np.sqrt(g2)

def v_pseudo_iso(r: np.ndarray, params):
    r_c, V0 = params
    x = r / r_c
    return V0 * g_pseudo_iso(x)

# ---------------- Burkert profile ----------------

def g_burkert(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x = np.where(x == 0.0, 1e-8, x)
    term = 0.5 * np.log(1.0 + x**2) + np.log(1.0 + x) - np.arctan(x)
    g2 = term / x
    g2 = np.maximum(g2, 0.0)
    return np.sqrt(g2)

def v_burkert(r: np.ndarray, params):
    r0, V0 = params
    x = r / r0
    return V0 * g_burkert(x)

# ---------------- Einasto profile ----------------

def v_einasto(r: np.ndarray, params, alpha: float = 0.18) -> np.ndarray:
    """Einasto profile with fixed alpha (2-parameter version: r_s, V0).

    We compute a dimensionless shape g_E(x) numerically via mass integration.
    This is slower than the analytic forms but fine for modest samples.
    """
    r_s, V0 = params
    r = np.asarray(r, dtype=float)
    if np.all(r == 0):
        return np.zeros_like(r)
    x = r / r_s
    xmax = max(10.0, float(x.max()) * 1.2)
    x_grid = np.linspace(1e-4, xmax, 1200)
    rho = np.exp(-(2.0 / alpha) * (x_grid**alpha - 1.0))
    dx = x_grid[1] - x_grid[0]
    M = np.cumsum(rho * x_grid**2) * dx
    M_interp = np.interp(x, x_grid, M)
    x_safe = np.where(x == 0.0, 1e-8, x)
    g2 = M_interp / x_safe
    g2 = np.maximum(g2, 0.0)
    g = np.sqrt(g2)
    max_g = g.max()
    if max_g > 0:
        g = g / max_g
    return V0 * g

def make_einasto_model(alpha: float = 0.18) -> HaloModel:
    def _v(r, params):
        return v_einasto(r, params, alpha=alpha)
    return HaloModel(name=f"Einasto_alpha{alpha:.2f}", vfunc=_v)

# Convenience list of standard models (not including Einasto by default)

def get_default_models(include_einasto: bool = True):
    models = [
        HaloModel(name="GenFT", vfunc=v_genft),
        HaloModel(name="NFW", vfunc=v_nfw),
        HaloModel(name="Burkert", vfunc=v_burkert),
        HaloModel(name="PseudoIso", vfunc=v_pseudo_iso),
    ]
    if include_einasto:
        models.append(make_einasto_model(alpha=0.18))
    return models
