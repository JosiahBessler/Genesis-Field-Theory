import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import matplotlib.pyplot as plt

from halo_models import get_default_models  # GenFT, NFW, Burkert, PseudoIso, Einasto

BASE_DIR = os.path.dirname(__file__)
ROT_DIR = os.path.join(BASE_DIR, "Rotmod_LTG")
RESULTS_OUT = os.path.join(BASE_DIR, "halo_dm_baryons_results.csv")
PLOTS_DIR = os.path.join(BASE_DIR, "plots_dm_baryons")


def load_rotmod_file(galaxy_name: str):
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
    # Accept either "<name>.dat" or "<name>_rotmod.dat"
    candidates = [
        f"{galaxy_name}.dat",
        f"{galaxy_name}_rotmod.dat",
    ]
    path = None
    for fn in candidates:
        p = os.path.join(ROT_DIR, fn)
        if os.path.exists(p):
            path = p
            break
    if path is None:
        raise FileNotFoundError(f"Could not find rotmod file for '{galaxy_name}'")

    R_list, Vobs_list, Verr_list = [], [], []
    Vgas_list, Vdisk_list, Vbul_list = [], [], []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            # Expect at least: R  Vobs  Verr  Vgas  Vdisk  Vbul
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

    # Guard against zero or negative errors
    if np.any(Verr <= 0):
        positive = Verr[Verr > 0]
        default_err = np.median(positive) if positive.size > 0 else 5.0
        Verr[Verr <= 0] = default_err

    return R, Vobs, Verr, Vgas, Vdisk, Vbul


def chi2_dm_baryons(
    R,
    Vobs,
    Verr,
    Vgas,
    Vdisk,
    Vbul,
    vfunc,
    params,
):
    """
    Chi^2 for a generic DM halo + baryons model.

    params = (r_s, V0, Yd)
      r_s : halo scale radius
      V0  : halo velocity scale
      Yd  : mass-to-light–like scaling applied to both disk and bulge
    """
    r_s, V0, Yd = params
    if r_s <= 0 or V0 <= 0 or Yd < 0:
        return 1e30

    # DM halo
    Vdm = vfunc(R, (r_s, V0))

    # Baryons (gas fixed, disk+bulge scaled by Yd)
    Vbar2 = Vgas**2 + (Yd * Vdisk)**2 + (Yd * Vbul)**2
    Vmodel = np.sqrt(Vdm**2 + Vbar2)

    chi2 = np.sum(((Vobs - Vmodel) / Verr) ** 2)
    return float(chi2)


def fit_dm_baryons_for_galaxy_and_model(
    galaxy_name: str,
    model,
    R,
    Vobs,
    Verr,
    Vgas,
    Vdisk,
    Vbul,
):
    """
    Fit DM + baryons (with free Yd) for a single galaxy and halo model.

    model: a HaloModel from halo_models.get_default_models()

    Returns a dict with best-fit params and chi2/dof, or None on failure.
    """
    if R.size < 4:
        return None

    # Initial guesses
    r_s_init = np.median(R)
    V0_init = np.max(Vobs)
    Yd_init = 1.0

    def objective(p):
        return chi2_dm_baryons(
            R, Vobs, Verr, Vgas, Vdisk, Vbul, model.vfunc, p
        )

    res = minimize(
        objective,
        x0=np.array([r_s_init, V0_init, Yd_init]),
        method="Nelder-Mead",
        options={"maxiter": 5000, "xatol": 1e-4, "fatol": 1e-4},
    )

    r_s_best, V0_best, Yd_best = res.x
    chi2 = objective(res.x)
    dof = max(int(R.size) - 3, 1)  # 3 free params: r_s, V0, Yd
    chi2_dof = chi2 / dof

    return {
        "galaxy": galaxy_name,
        "model": model.name,
        "r_scale": float(r_s_best),
        "V0": float(V0_best),
        "Yd": float(Yd_best),
        "chi2": float(chi2),
        "dof": int(dof),
        "chi2_dof": float(chi2_dof),
        "success": bool(res.success),
        "nfev": int(res.nfev),
    }


def make_plot(
    galaxy_name,
    R,
    Vobs,
    Verr,
    Vgas,
    Vdisk,
    Vbul,
    model_results,
):
    """
    Create a PNG plot for one galaxy with all halo models overlaid.
    model_results: list of dicts with keys:
        'model', 'r_scale', 'V0', 'Yd', plus HaloModel object in 'model_obj'.
    """

    if R.size == 0:
        return

    # Make a fine radius grid for smooth curves
    R_plot = np.linspace(R.min(), R.max(), 200)

    plt.figure(figsize=(6, 4.5))

    # Plot data + error bars
    plt.errorbar(
        R,
        Vobs,
        yerr=Verr,
        fmt="o",
        markersize=4,
        linewidth=1,
        capsize=2,
        label="data",
    )

    # Some simple color/linestyle choices
    color_map = {
        "GenFT": "tab:orange",
        "NFW": "tab:green",
        "Burkert": "tab:red",
        "PseudoIso": "tab:purple",
    }
    # Einasto may have name like 'Einasto_alpha0.18'
    default_color = "tab:blue"

    for res in model_results:
        model_name = res["model"]
        model_obj = res["model_obj"]
        r_s = res["r_scale"]
        V0 = res["V0"]
        Yd = res["Yd"]

        # Recompute model curve on R_plot
        Vdm = model_obj.vfunc(R_plot, (r_s, V0))
        Vbar2 = (
            np.interp(R_plot, R, Vgas) ** 2
            + (Yd * np.interp(R_plot, R, Vdisk)) ** 2
            + (Yd * np.interp(R_plot, R, Vbul)) ** 2
        )
        Vmodel = np.sqrt(Vdm**2 + Vbar2)

        color = default_color
        for key, c in color_map.items():
            if model_name.startswith(key):
                color = c
                break

        linestyle = "-"
        if "NFW" in model_name:
            linestyle = "--"
        elif "Burkert" in model_name:
            linestyle = "-."
        elif "PseudoIso" in model_name:
            linestyle = ":"

        plt.plot(
            R_plot,
            Vmodel,
            linestyle=linestyle,
            color=color,
            label=model_name,
            linewidth=1.5,
        )

    plt.xlabel("radius [kpc]")
    plt.ylabel("V_c(r) [km/s]")
    plt.title(f"{galaxy_name}: DM + baryons (all halo models)")
    plt.legend(fontsize=8)
    plt.tight_layout()

    os.makedirs(PLOTS_DIR, exist_ok=True)
    out_path = os.path.join(PLOTS_DIR, f"{galaxy_name}_dm_baryons.png")
    plt.savefig(out_path, dpi=150)
    plt.close()


def main():
    os.makedirs(PLOTS_DIR, exist_ok=True)

    # Build galaxy list from files in Rotmod_LTG
    files = [
        fn
        for fn in os.listdir(ROT_DIR)
        if fn.endswith("_rotmod.dat") or fn.endswith(".dat")
    ]
    # Strip extensions
    galaxies = sorted(set(os.path.splitext(fn)[0].replace("_rotmod", "") for fn in files))
    print(f"Found {len(galaxies)} galaxies in Rotmod_LTG.")

    models = get_default_models(include_einasto=True)
    results = []

    total_tasks = len(galaxies) * len(models)
    task_idx = 0

    for gal in galaxies:
        # Load data once per galaxy
        try:
            R, Vobs, Verr, Vgas, Vdisk, Vbul = load_rotmod_file(gal)
        except Exception as e:
            print(f"ERROR loading {gal}: {e}")
            continue

        if R.size < 4:
            print(f"{gal}: not enough data points, skipping.")
            continue

        model_results_for_plot = []

        for model in models:
            task_idx += 1
            print(f"[{task_idx}/{total_tasks}] {gal} - {model.name}")
            try:
                res = fit_dm_baryons_for_galaxy_and_model(
                    gal, model, R, Vobs, Verr, Vgas, Vdisk, Vbul
                )
                if res is not None:
                    res["model_obj"] = model  # attach for plotting
                    results.append({k: v for k, v in res.items() if k != "model_obj"})
                    model_results_for_plot.append(res)
                    print(
                        f"    chi2/dof={res['chi2_dof']:.3f}, "
                        f"Yd={res['Yd']:.3f}, "
                        f"success={res['success']}"
                    )
                else:
                    print("    Not enough data points, skipping this model.")
            except Exception as e:
                print(f"    ERROR fitting {gal} with {model.name}: {e}")

        # After fitting all models for this galaxy, make plot
        try:
            if model_results_for_plot:
                make_plot(
                    gal,
                    R,
                    Vobs,
                    Verr,
                    Vgas,
                    Vdisk,
                    Vbul,
                    model_results_for_plot,
                )
        except Exception as e:
            print(f"    ERROR plotting {gal}: {e}")

    if not results:
        print("No fits produced.")
        return

    out_df = pd.DataFrame(results)
    out_df.to_csv(RESULTS_OUT, index=False)
    print(f"Saved DM + baryons fit results for all models to: {RESULTS_OUT}")
    print(f"PNG plots saved to: {PLOTS_DIR}")


if __name__ == "__main__":
    main()
