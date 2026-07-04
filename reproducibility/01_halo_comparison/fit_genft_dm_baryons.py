import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize

from halo_models import v_genft  # uses your Template C

BASE_DIR = os.path.dirname(__file__)
ROT_DIR = os.path.join(BASE_DIR, "Rotmod_LTG")
RESULTS_IN = os.path.join(BASE_DIR, "halo_fit_results.csv")
RESULTS_OUT = os.path.join(BASE_DIR, "genft_dm_baryons_results.csv")


def load_rotmod_file(galaxy_name):
    """
    Load a SPARC Rotmod_LTG file for a given galaxy name (without .dat).

    Returns:
        R       : radii [kpc]
        Vobs    : observed rotation speed [km/s]
        Verr    : velocity uncertainty [km/s]
        Vgas    : gas contribution [km/s]
        Vdisk   : stellar disk contribution [km/s]
        Vbul    : bulge contribution [km/s]
    """
    fname = f"{galaxy_name}.dat"
    path = os.path.join(ROT_DIR, fname)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Could not find rotmod file: {path}")

    R_list = []
    Vobs_list = []
    Verr_list = []
    Vgas_list = []
    Vdisk_list = []
    Vbul_list = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            R_list.append(float(parts[0]))
            Vobs_list.append(float(parts[1]))
            Verr_list.append(abs(float(parts[2])) if float(parts[2]) != 0 else 5.0)
            Vgas_list.append(float(parts[3]))
            Vdisk_list.append(float(parts[4]))
            Vbul_list.append(float(parts[5]))

    R = np.array(R_list)
    Vobs = np.array(Vobs_list)
    Verr = np.array(Verr_list)
    Vgas = np.array(Vgas_list)
    Vdisk = np.array(Vdisk_list)
    Vbul = np.array(Vbul_list)

    # Guard against zero errors
    if np.any(Verr <= 0):
        positive = Verr[Verr > 0]
        default_err = np.median(positive) if positive.size > 0 else 5.0
        Verr[Verr <= 0] = default_err

    return R, Vobs, Verr, Vgas, Vdisk, Vbul


def chi2_genft_dm_baryons(R, Vobs, Verr, Vgas, Vdisk, Vbul, params):
    """
    Compute chi^2 for GenFT DM + baryons model.

    params = (r0, V0, Yd)
      r0 : DM halo scale radius
      V0 : DM halo velocity scale
      Yd : mass-to-light–like scaling applied to both disk and bulge
    """
    r0, V0, Yd = params
    # Enforce positive params
    if r0 <= 0 or V0 <= 0 or Yd < 0:
        return 1e30

    # DM halo
    Vdm = v_genft(R, (r0, V0))

    # Baryonic contribution (gas fixed, disk+bulge scaled by Yd)
    Vbar2 = Vgas**2 + (Yd * Vdisk)**2 + (Yd * Vbul)**2

    Vmodel = np.sqrt(Vdm**2 + Vbar2)

    chi2 = np.sum(((Vobs - Vmodel) / Verr) ** 2)
    return float(chi2)


def fit_genft_dm_baryons_for_galaxy(galaxy_name):
    """
    Fit GenFT DM + baryons (with free Yd) for a single galaxy.

    Returns a dict with best-fit params and chi2/dof.
    """
    R, Vobs, Verr, Vgas, Vdisk, Vbul = load_rotmod_file(galaxy_name)

    if R.size < 4:
        return None  # Not enough data points

    # Initial guesses
    r0_init = np.median(R)
    V0_init = np.max(Vobs)
    Yd_init = 1.0

    def objective(p):
        return chi2_genft_dm_baryons(R, Vobs, Verr, Vgas, Vdisk, Vbul, p)

    res = minimize(
        objective,
        x0=np.array([r0_init, V0_init, Yd_init]),
        method="Nelder-Mead",
        options={"maxiter": 5000, "xatol": 1e-4, "fatol": 1e-4},
    )

    r0_best, V0_best, Yd_best = res.x
    chi2 = objective(res.x)
    dof = max(int(R.size) - 3, 1)
    chi2_dof = chi2 / dof

    return {
        "galaxy": galaxy_name,
        "r0": float(r0_best),
        "V0": float(V0_best),
        "Yd": float(Yd_best),
        "chi2": float(chi2),
        "dof": int(dof),
        "chi2_dof": float(chi2_dof),
        "success": bool(res.success),
        "nfev": int(res.nfev),
    }


def main():
    # Load halo_fit_results.csv and get ALL galaxies
    df = pd.read_csv(RESULTS_IN)
    galaxies = sorted(df["galaxy"].unique())
    print(f"Found {len(galaxies)} galaxies total.")

    results = []

    for i, gal in enumerate(galaxies, 1):
        print(f"[{i}/{len(galaxies)}] Fitting GenFT+bar for {gal} ...")
        try:
            res = fit_genft_dm_baryons_for_galaxy(gal)
            if res is not None:
                results.append(res)
                print(f"    chi2/dof = {res['chi2_dof']:.3f}, Yd = {res['Yd']:.3f}")
            else:
                print("    Not enough data points, skipping.")
        except Exception as e:
            print(f"    ERROR fitting {gal}: {e}")

    if not results:
        print("No fits produced.")
        return

    out_df = pd.DataFrame(results)
    out_df.to_csv(RESULTS_OUT, index=False)
    print(f"Saved GenFT DM + baryons fit results to: {RESULTS_OUT}")



if __name__ == "__main__":
    main()
