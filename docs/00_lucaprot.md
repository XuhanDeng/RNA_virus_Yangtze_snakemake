# 00_lucaprot

## Purpose

Installs the LucaProt deep-learning RdRP predictor from GitHub and downloads/bootstraps its pre-trained model database by running a single test prediction. Must be run on the login node (requires internet access).

## Rules

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `install_lucaprot` | Clones the LucaProt repository if not already present; creates a marker file on success | — | `config["lucaprot"]["marker_install"]` |
| `download_lucaprot_db` | Bootstraps the model weights by running `predict_one_sample.py` on a dummy protein sequence, which triggers PyTorch Hub to download the required model files | `marker_install` | `config["lucaprot"]["marker_db"]` |

## Key Config Parameters

| Config key | Description |
|-----------|-------------|
| `lucaprot.marker_install` | Marker file path confirming repo clone is complete |
| `lucaprot.marker_db` | Marker file path confirming database/weights download is complete |
| `lucaprot.install_dir` | Directory where LucaProt repository is cloned |
| `lucaprot.repo_url` | Git URL for the LucaProt repository |
| `lucaprot.db_dir` | Directory where model weights (PyTorch Hub cache) are stored |

## Dependencies / Tools Used

- `git` — clone the LucaProt repository
- `python` (`predict_one_sample.py`) — triggers model weight download via PyTorch Hub
- Conda environment: `../envs/lucaprot.yaml`

## Notes

- This workflow is designed to run on the **login node** only; compute nodes typically lack internet access.
- The download step works by running a real prediction on a hardcoded dummy amino-acid sequence. There is no dedicated `download` command; the model files are fetched on first use.
- If the repo directory already contains `src/predict_one_sample.py`, the clone step is silently skipped (idempotent).
- The `marker_db` touch file is created using absolute paths because the script must `cd` into the install directory before running Python.
