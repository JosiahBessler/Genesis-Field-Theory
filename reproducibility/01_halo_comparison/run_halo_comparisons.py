
import os
import numpy as np
import pandas as pd
from halo_models import get_default_models
from fitting import fit_halo_model

DATA_DIR = os.path.join(os.path.dirname(__file__), 'Rotmod_LTG')
OUTPUT_CSV = os.path.join(os.path.dirname(__file__), 'halo_fit_results.csv')

def load_galaxy(filepath):
    """Load a SPARC Rotmod_LTG .dat file.

    Returns R, V_obs, sigma_V as numpy arrays.
    """
    R_vals = []
    V_obs_vals = []
    err_vals = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            R_vals.append(float(parts[0]))
            V_obs_vals.append(float(parts[1]))
            err_vals.append(abs(float(parts[2])) if float(parts[2]) != 0 else 5.0)
    R = np.array(R_vals)
    V_obs = np.array(V_obs_vals)
    sigma_V = np.array(err_vals)
    # Guard against zero or tiny errors
    sigma_V[sigma_V <= 0] = max(1.0, np.median(sigma_V))
    return R, V_obs, sigma_V

def fit_all_for_galaxy(filename):
    filepath = os.path.join(DATA_DIR, filename)
    R, V_obs, sigma_V = load_galaxy(filepath)
    if len(R) <= 3:
        return []
    r_guess = np.median(R)
    V_guess = np.max(V_obs)
    models = get_default_models(include_einasto=True)
    results = []
    for model in models:
        res = fit_halo_model(R, V_obs, sigma_V, model, r_guess, V_guess)
        res['galaxy'] = os.path.splitext(filename)[0]
        results.append(res)
    return results

def main():
    files = sorted([fn for fn in os.listdir(DATA_DIR) if fn.endswith('_rotmod.dat')])
    all_results = []
    for i, fn in enumerate(files, 1):
        print(f"Fitting {fn} ({i}/{len(files)})...")
        try:
            galaxy_results = fit_all_for_galaxy(fn)
            all_results.extend(galaxy_results)
        except Exception as e:
            print(f"  Error fitting {fn}: {e}")
    if not all_results:
        print("No results produced.")
        return
    df = pd.DataFrame(all_results)
    # Determine best model per galaxy by lowest chi2
    df['rank'] = df.groupby('galaxy')['chi2'].rank(method='first')
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"Saved results to {OUTPUT_CSV}")

if __name__ == '__main__':
    main()
