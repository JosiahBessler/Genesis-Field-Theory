"""Generate plausible dummy GenesisFT inputs for testing the pipeline."""
from __future__ import annotations

import pathlib
import numpy as np
import pandas as pd


def main() -> None:
    rng = np.random.default_rng(42)
    n = 600
    lambdas = np.sort(rng.lognormal(mean=0.0, sigma=0.6, size=n))
    mus = rng.lognormal(mean=0.0, sigma=0.3, size=n)

    data_dir = pathlib.Path("data")
    data_dir.mkdir(exist_ok=True)

    np.savez(data_dir / "dummy_genft.npz", lambdas=lambdas, mus=mus)
    df = pd.DataFrame({"lambda": lambdas, "mu": mus})
    df.to_csv(data_dir / "dummy_genft.csv", index=False)
    print(f"Wrote {data_dir / 'dummy_genft.npz'} and {data_dir / 'dummy_genft.csv'}")


if __name__ == "__main__":
    main()
