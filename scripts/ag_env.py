"""Shared setup for AlphaGenome scripts. Import at the top of each one.
These scripts run in the `alphagenome` conda env, NOT `klc2`."""

import os, sys

REQUIRED_ENV = "alphagenome"


def check_env():
    active = os.environ.get("CONDA_DEFAULT_ENV", "<none>")
    if active != REQUIRED_ENV:
        sys.exit(
            f"ERROR: this script must run in the '{REQUIRED_ENV}' conda env "
            f"(currently '{active}').\n"
            f"  conda activate {REQUIRED_ENV}\n"
            f"  python3 {sys.argv[0]}\n"
            f"  conda activate klc2      # switch back when done"
        )
    if "ALPHA_GENOME_API_KEY" not in os.environ:
        sys.exit(
            "ERROR: ALPHA_GENOME_API_KEY is not set.\n"
            "  echo 'export ALPHA_GENOME_API_KEY=\"your_key\"' >> ~/.bashrc\n"
            "  source ~/.bashrc"
        )
    return os.environ["ALPHA_GENOME_API_KEY"]
