"""
Download sequences from NCBI in batches.
Uses efetch -db nuccore -id directly so WGS contig accessions (e.g. CASDWV010000007)
are fetched by accession rather than term-search.

Batch size is kept at 50 to match efetch's internal sub-batch limit, reducing 502 errors.
Empty/failed batches fall back to per-accession downloads.
"""
import subprocess
import sys
import time

sys.stdout = open(snakemake.log.log, "w")
sys.stderr = open(snakemake.log.err, "w")

BATCH_SIZE  = 50
MAX_RETRIES = 5
SLEEP_SEC   = 0.4   # NCBI rate limit: 3 req/s without API key

accessions_file = snakemake.input[0]
out_fasta       = snakemake.output[0]

with open(accessions_file) as f:
    accessions = [line.strip() for line in f if line.strip()]

batches = [accessions[i:i + BATCH_SIZE] for i in range(0, len(accessions), BATCH_SIZE)]
print(f"Total accessions: {len(accessions)}, batches: {len(batches)}", flush=True)

failed = []


def fetch(ids: list[str]) -> str:
    result = subprocess.run(
        ["efetch", "-db", "nuccore", "-id", ",".join(ids), "-format", "fasta"],
        capture_output=True, text=True, check=True
    )
    output = result.stdout.strip()
    if not output or "EMPTY RESULT" in output or not output.startswith(">"):
        raise ValueError(f"efetch returned no sequences for {len(ids)} ids")
    return result.stdout


with open(out_fasta, "w") as out:
    for idx, batch in enumerate(batches, 1):
        fetched = False
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                out.write(fetch(batch))
                print(f"Batch {idx}/{len(batches)} done ({len(batch)} accessions)", flush=True)
                fetched = True
                break
            except (subprocess.CalledProcessError, ValueError) as e:
                wait = min(2 ** attempt, 60)
                print(
                    f"Batch {idx} attempt {attempt} failed: {e} — retrying in {wait}s",
                    flush=True
                )
                time.sleep(wait)

        if not fetched:
            print(f"Batch {idx} exhausted retries — falling back to per-accession", flush=True)
            for acc in batch:
                for attempt in range(1, MAX_RETRIES + 1):
                    try:
                        out.write(fetch([acc]))
                        break
                    except (subprocess.CalledProcessError, ValueError):
                        if attempt == MAX_RETRIES:
                            print(f"SKIP: {acc} could not be fetched", file=sys.stderr)
                            failed.append(acc)
                        else:
                            time.sleep(min(2 ** attempt, 60))
                time.sleep(SLEEP_SEC)

        time.sleep(SLEEP_SEC)

if failed:
    print(f"\nFailed accessions ({len(failed)}):", file=sys.stderr)
    for acc in failed:
        print(f"  {acc}", file=sys.stderr)

print(f"Done. Output written to {out_fasta}", flush=True)
