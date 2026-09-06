# 00_setup_RNA_database

## Purpose

One-time setup workflow that downloads and prepares all external databases required by the Yangtze RNA virus pipeline: SortMeRNA rRNA reference databases, RdRpCATCH, palm_annot, RVMT, ICTV Riboviria sequences, ESViritu, and Kraken2. Must be run on the login node (requires internet access).

## Rules

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `download_sortmerna_db` | Downloads and gunzips four SortMeRNA reference database FASTA files (default, fast, sensitive, RFAM) | URLs from config | Four FASTA files at `config["sortmerna"]["dbs"][0-3]` |
| `download_rdrpcatch` | Uses the `rdrpcatch databases` CLI to fetch the RdRpCATCH HMM database | — | `config["rdrp_catch"]["db_dir"]` (directory) |
| `install_palm_annot` | Clones the palm_annot repository from GitHub and makes binaries executable | URL | `config["palm_annot_install"]["install_dir"]` (directory) |
| `download_RVMT` | Downloads the RVMT zip archive | URL | `config["RVMT"]["zip"]` |
| `extract_RVMT` | Unzips the RVMT archive | `RVMT.zip` | `config["RVMT"]["dir"]` (directory) |
| `download_ICTV_vmr` | Downloads the ICTV Virus Metadata Resource Excel file | URL | `config["ICTV"]["vmr_xlsx"]` |
| `extract_ICTV_accessions` | Parses the VMR spreadsheet and filters Riboviria accessions by realm and coverage | `vmr_xlsx` | `config["ICTV"]["accessions_txt"]` |
| `download_ICTV_sequences` | Downloads GenBank sequences for the extracted accessions using Entrez Direct | `accessions_txt` | `config["ICTV"]["fasta"]` |
| `download_esviritu_db` | Downloads and extracts the ESViritu virus database tarball | URL | `config["esviritu_db"]["dir"]` (directory) |
| `download_kraken2_db` | Downloads and extracts a pre-built Kraken2 database tarball | URL | `config["kraken2_db"]["dir"]` (directory) |

## Key Config Parameters

| Config key | Description |
|-----------|-------------|
| `sortmerna["dbs"]` | List of four SortMeRNA FASTA database paths |
| `sortmerna_db["smr_default_url"]` / `smr_fast_url` / `smr_sensitive_url` / `smr_rfam_url` | Download URLs for SortMeRNA databases |
| `rdrp_catch["db_dir"]` | RdRpCATCH database destination directory |
| `palm_annot_install["install_dir"]` / `["url"]` | palm_annot install directory and Git URL |
| `RVMT["url"]` / `["zip"]` / `["dir"]` | RVMT download URL, zip path, and extraction directory |
| `ICTV["vmr_url"]` / `["vmr_xlsx"]` / `["vmr_sheet"]` / `["realm_filter"]` / `["coverage_filter"]` / `["accessions_txt"]` / `["fasta"]` | ICTV VMR download parameters and output paths |
| `esviritu_db["url"]` / `["tar"]` / `["dir"]` | ESViritu DB download parameters |
| `kraken2_db["url"]` / `["tar"]` / `["dir"]` | Kraken2 DB download parameters |
| `setup_log_dir` | Root directory for setup log files |

## Dependencies / Tools Used

- `wget` — all HTTP downloads
- `gunzip` / `unzip` / `tar` — archive extraction
- `git` + `chmod` — palm_annot installation
- `rdrpcatch databases` — RdRpCATCH CLI (conda: `rdrp_catch.yaml`)
- Python script `extract_ICTV_accessions.py` with `openpyxl` (conda: `openpyxl.yaml`)
- Python script `download_ICTV_sequences.py` with Entrez Direct (conda: `entrez_direct.yaml`)

## Notes

- Run on the **login node** only; compute nodes on the HPC cluster lack internet access.
- `wget -c` is used throughout so interrupted downloads can be resumed.
- The `extract_ICTV_accessions` rule applies both a `realm_filter` (Riboviria) and a `coverage_filter` before writing accession IDs; the exact filter logic is in `scripts/00_setup_RNA_database/extract_ICTV_accessions.py`.
- palm_annot uses precompiled Linux binaries (not conda); will not work on macOS.
- The ESViritu tarball is extracted with `--strip-components=1` so the top-level archive directory is removed.
- No Kraken2 database *building* rule is present in this file (the comment placeholder is empty); downloading a prebuilt tarball is assumed.
