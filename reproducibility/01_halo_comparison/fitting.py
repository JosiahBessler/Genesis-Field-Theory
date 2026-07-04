
import numpy as np
from typing import Tuple, Dict
from scipy.optimize import minimize
from halo_models import HaloModel

def chi2_model(R: np.ndarray,
               V_obs: np.ndarray,
               sigma_V: np.ndarray,
               model: HaloModel,
               params: Tuple[float, float]) -> Tuple[float, int]:
    V_model = model.vfunc(R, params)
    chi2 = np.sum(((V_obs - V_model) / sigma_V)**2)
    dof = len(R) - model.n_params
    return float(chi2), int(dof)

def fit_halo_model(R: np.ndarray,
                   V_obs: np.ndarray,
                   sigma_V: np.ndarray,
                   model: HaloModel,
                   r_init: float,
                   V_init: float) -> Dict:
    def objective(p):
        r_s, V0 = p
        if r_s <= 0.0 or V0 <= 0.0:
            return 1e30
        chi2, _ = chi2_model(R, V_obs, sigma_V, model, (r_s, V0))
        return chi2

    res = minimize(objective,
                   x0=np.array([r_init, V_init]),
                   method="Nelder-Mead",
                   options={"maxiter": 2000, "xatol": 1e-4, "fatol": 1e-4})

    best_r, best_V = res.x
    chi2, dof = chi2_model(R, V_obs, sigma_V, model, (best_r, best_V))
    chi2_dof = chi2 / dof if dof > 0 else np.nan
    return {
        "model": model.name,
        "r_scale": float(best_r),
        "V0": float(best_V),
        "chi2": float(chi2),
        "dof": int(dof),
        "chi2_dof": float(chi2_dof),
        "success": bool(res.success),
        "nfev": int(res.nfev),
    }
