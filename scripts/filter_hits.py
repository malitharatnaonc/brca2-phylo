"""
Collapse per-domain hits to per-protein records, filter by domain coverage
(better signal than raw E-value for a protein as divergent as BRCA2), cap
the hit count, and fetch full sequences for whatever survives.

Run locally (laptop).
"""
import os
import sys
import time
import requests
import pandas as pd

UNIPROT_FASTA_API = "https://rest.uniprot.org/uniprotkb/{}.fasta"
# UniParc entries (ID starts with "UPI", e.g. a UniRef cluster whose
# representative isn't in UniProtKB) live in a separate database with their
# own endpoint - confirmed against UniProt's own OpenAPI spec
# (rest.uniprot.org/uniparc/api/docs lists "/uniparc/{upi}" supporting
# "text/plain;format=fasta"), not guessed. The uniprotkb endpoint will
# always 404 for these, no matter how many retries - 335 of the 648
# failures in this project were UPI IDs, which is why retries alone didn't
# help them.
UNIPARC_FASTA_API = "https://rest.uniprot.org/uniparc/{}.fasta"
UNIPROTKB_JSON_API = "https://rest.uniprot.org/uniprotkb/{}.json"

# Core conserved domains of human BRCA2 (P51587, 3418 aa), from UniProt/
# InterPro feature annotations: 8 BRC repeats (RAD51-binding) + the
# C-terminal DNA-binding domain (helical + OB-fold domains + Tower).
# ~989 of 3418 residues (~29%) - the rest is the intrinsically disordered
# majority of the protein that motivated this whole project. A hit's
# `approx_coverage` (span relative to the FULL 3418-residue query) penalises
# any hit that doesn't also align those disordered regions - which, for a
# genuinely divergent ortholog, it correctly won't, since disordered linkers
# don't carry alignable sequence signal even in a real ortholog. That
# systematically discards exactly the distant, structurally-interesting
# homologs this project is trying to find. `core_domain_coverage` instead
# measures what fraction of the CONSERVED, FUNCTIONALLY MEANINGFUL residues
# a hit actually covers - the metric that should drive filtering.
CORE_DOMAINS = [
    (1002, 1036), (1212, 1246), (1421, 1454), (1517, 1550),
    (1664, 1696), (1837, 1869), (1972, 2005), (2051, 2084),  # BRC1-8
    (2474, 3190),  # DBD (helical + OB1 + Tower + OB2/OB3)
]
CORE_TOTAL = sum(end - start + 1 for start, end in CORE_DOMAINS)


def _interval_union(intervals: list[tuple[int, int]]) -> list[list[int]]:
    """Merge overlapping/adjacent (start, end) intervals into a minimal set."""
    if not intervals:
        return []
    merged = [list(intervals[0])]
    for start, end in sorted(intervals[1:]):
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


def _overlap_length(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def _core_domain_coverage(hit_intervals: list[tuple[int, int]]) -> float:
    """Fraction of the CORE_TOTAL conserved-domain residues actually covered
    by this hit's aligned segments (one (hmm_from, hmm_to) pair per domain
    hit row for this protein - BRCA2's own multi-domain architecture means a
    real hit often legitimately has several separate aligned segments, not
    one contiguous span)."""
    hit_union = _interval_union(hit_intervals)
    covered = sum(
        _overlap_length(core_start, core_end, h_start, h_end)
        for core_start, core_end in CORE_DOMAINS
        for h_start, h_end in hit_union
    )
    return covered / CORE_TOTAL


def aggregate_coverage(df: pd.DataFrame) -> pd.DataFrame:
    """Per uniprot_id: both the naive whole-length span coverage (kept for
    reference/debugging) and the core-domain coverage (the one filtering
    should actually use - see CORE_DOMAINS comment above)."""
    grouped = df.groupby("uniprot_id").agg(
        organism_name=("organism_name", "first"),
        taxon_id=("taxon_id", "first"),
        lineage=("lineage", "first"),
        sequence_length=("sequence_length", "first"),
        n_domains=("dom_num", "count"),
        best_evalue=("i_evalue", "min"),
        min_hmm_from=("hmm_from", "min"),
        max_hmm_to=("hmm_to", "max"),
        qlen=("qlen", "first"),
    ).reset_index()
    grouped["approx_coverage"] = (
        (grouped["max_hmm_to"] - grouped["min_hmm_from"] + 1) / grouped["qlen"]
    ).clip(upper=1.0)

    core_cov = (
        df.groupby("uniprot_id")
        .apply(lambda g: _core_domain_coverage(list(zip(g["hmm_from"], g["hmm_to"]))),
               include_groups=False)
        .rename("core_domain_coverage")
    )
    grouped = grouped.merge(core_cov, on="uniprot_id")
    return grouped


def _fetch_fasta_once(url: str, timeout: int = 15) -> tuple[str | None, str | None]:
    """Single raw attempt against a fasta-format endpoint. Returns
    (sequence, error_reason); error_reason is None only on real success."""
    try:
        resp = requests.get(url, timeout=timeout)
        if resp.status_code == 404:
            return None, "404 not found"
        resp.raise_for_status()
        lines = resp.text.strip().splitlines()
        if not lines:
            return None, "empty response"
        return "".join(lines[1:]), None  # drop header
    except requests.exceptions.HTTPError as e:
        if e.response is not None and e.response.status_code == 429:
            return None, "429 rate limited"
        return None, f"HTTP {e.response.status_code if e.response is not None else '?'}"
    except requests.exceptions.Timeout:
        return None, "timeout"
    except Exception as e:
        return None, f"{type(e).__name__}: {e}"


def _resolve_uniparc_id(accession: str) -> str | None:
    """A UniProtKB accession can return HTTP 200 with an empty fasta body
    when the entry has been deleted/demerged from UniProtKB's curated set
    (e.g. 'Not part of a reference proteome') - confirmed on a real failing
    accession, whose JSON entry showed entryType":"Inactive" plus
    extraAttributes.uniParcId. The underlying sequence is NOT gone - UniParc
    (UniProt's permanent sequence archive) never deletes anything, even
    when UniProtKB prunes an entry from its reference set. This looks up
    that UniParc ID so the sequence can still be fetched."""
    try:
        resp = requests.get(UNIPROTKB_JSON_API.format(accession), timeout=15)
        resp.raise_for_status()
        return resp.json().get("extraAttributes", {}).get("uniParcId")
    except Exception:
        return None


def fetch_sequence(accession: str, retries: int = 3, backoff: float = 10.0) -> tuple[str | None, str | None]:
    """Retries transient failures before giving up, with backoff long enough
    to plausibly clear a rate-limit window (10/20/30s, not an earlier
    2/4/6s that saw a 100% retry failure rate on a batch where a manual
    curl immediately after succeeded cleanly - consistent with the whole
    batch happening inside one still-active rate-limited window).

    Also falls back to UniParc when UniProtKB returns an empty body for an
    inactive/deleted entry (see _resolve_uniparc_id) - this fallback itself
    is retried through the same loop, so a transient failure while
    resolving the UniParc ID doesn't permanently give up either.

    Returns (sequence, error_reason) - error_reason is None on success,
    otherwise a short string identifying what actually went wrong, so a
    failure log is self-diagnosing instead of requiring a manual curl to
    find out why."""
    url = UNIPARC_FASTA_API.format(accession) if accession.startswith("UPI") \
        else UNIPROT_FASTA_API.format(accession)
    last_error = None
    for attempt in range(retries):
        seq, error = _fetch_fasta_once(url)
        if seq:
            return seq, None
        last_error = error

        if error == "empty response" and not accession.startswith("UPI"):
            uniparc_id = _resolve_uniparc_id(accession)
            if uniparc_id:
                seq2, error2 = _fetch_fasta_once(UNIPARC_FASTA_API.format(uniparc_id))
                if seq2:
                    return seq2, None
                last_error = f"empty response, UniParc fallback also failed: {error2}"

        if attempt < retries - 1:
            time.sleep(backoff * (attempt + 1))
    return None, last_error


def load_existing_sequences(fasta_path: str) -> dict[str, tuple[str, str]]:
    """Read any sequences a previous (partial/failed) run already fetched
    successfully, keyed by accession -> (organism_name, sequence). Lets a
    re-run skip everything that already succeeded and only retry what's
    actually missing, rather than re-fetching ~1000 sequences to recover 648
    failures."""
    existing: dict[str, tuple[str, str]] = {}
    if not os.path.exists(fasta_path):
        return existing
    acc, organism, seq_lines = None, "", []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if acc:
                    existing[acc] = (organism, "".join(seq_lines))
                header = line[1:].split(" ", 1)
                acc = header[0]
                organism = header[1] if len(header) > 1 else ""
                seq_lines = []
            else:
                seq_lines.append(line)
        if acc:
            existing[acc] = (organism, "".join(seq_lines))
    return existing


def select_top_per_species(agg: pd.DataFrame, per_species_n: int) -> pd.DataFrame:
    """Within each species (taxon_id), keep only the top `per_species_n` hits
    by coverage (best_evalue as tiebreak). Without this, a global coverage
    ranking can let a cluster of near-identical close relatives (e.g. many
    primates all with near-complete, high-coverage BRCA2) crowd out weaker
    -but-real hits from genuinely divergent, taxonomically interesting
    species - exactly the ones this project cares most about. Hits with no
    resolved taxon_id (missing taxonomy lookup) are grouped together as one
    "unknown" bucket and capped the same way, as a safety net rather than an
    ideal outcome - worth spot-checking that bucket if it ends up large."""
    agg = agg.sort_values(["core_domain_coverage", "best_evalue"], ascending=[False, True])
    return agg.groupby("taxon_id", dropna=False, group_keys=False).head(per_species_n)


def filter_hits(df: pd.DataFrame, min_cov: float, max_hits: int,
                 per_species_n: int | None = None) -> pd.DataFrame:
    agg = aggregate_coverage(df)
    filtered = agg[agg["core_domain_coverage"] >= min_cov].copy()
    if per_species_n:
        filtered = select_top_per_species(filtered, per_species_n)
    filtered = filtered.sort_values("core_domain_coverage", ascending=False)
    filtered = filtered.head(max_hits)
    return filtered


if __name__ == "__main__":
    if "snakemake" in globals():
        in_path = snakemake.input["tax_table"]
        out_fasta = snakemake.output["filtered"]
        out_table = snakemake.output["filtered_table"]
        min_cov = snakemake.params["min_cov"]
        max_hits = snakemake.params["max_hits"]
        per_species_n = snakemake.params.get("per_species_n")
    else:
        in_path, out_fasta, out_table, min_cov, max_hits = (
            sys.argv[1], sys.argv[2], sys.argv[3],
            float(sys.argv[4]), int(sys.argv[5])
        )
        per_species_n = int(sys.argv[6]) if len(sys.argv) > 6 else None

    df = pd.read_csv(in_path, sep="\t")
    filtered = filter_hits(df, min_cov, max_hits, per_species_n)

    existing = load_existing_sequences(out_fasta)
    if existing:
        print(f"Found {len(existing)} already-fetched sequences from a previous "
              f"run at {out_fasta} - skipping those, only fetching what's missing.")

    failed: list[tuple[str, str]] = []  # (accession, reason)
    newly_fetched: dict[str, tuple[str, str]] = {}
    to_fetch = [row for _, row in filtered.iterrows() if row["uniprot_id"] not in existing]
    for i, row in enumerate(to_fetch):
        acc = row["uniprot_id"]
        seq, error = fetch_sequence(acc)
        if seq:
            newly_fetched[acc] = (row["organism_name"], seq)
        else:
            failed.append((acc, error or "unknown"))
        if i % 50 == 0:
            print(f"  fetched sequence {i}/{len(to_fetch)}")
        # 1s between requests, not the original 0.1s - a batch of 648 at 0.1s
        # saw a 100% retry failure rate that a manual curl immediately after
        # succeeded cleanly on, consistent with rate-limiting rather than
        # missing data. This is a cautious, unverified-exact-limit choice,
        # not a number UniProt publishes - if failures persist even at this
        # pace, that's worth treating as new evidence, not assuming this
        # guess was sufficient.
        time.sleep(1.0)

    with open(out_fasta, "w") as fh:
        n_written = 0
        for _, row in filtered.iterrows():
            acc = row["uniprot_id"]
            if acc in existing:
                organism, seq = existing[acc]
            elif acc in newly_fetched:
                organism, seq = newly_fetched[acc]
            else:
                continue
            fh.write(f">{acc} {organism}\n{seq}\n")
            n_written += 1

    if failed:
        failed_path = out_fasta.rsplit(".", 1)[0] + "_failed_accessions.txt"
        with open(failed_path, "w") as f:
            f.write("\n".join(f"{acc}\t{reason}" for acc, reason in failed))
        reason_counts = pd.Series([r for _, r in failed]).value_counts()
        print(f"WARNING: {len(failed)} accessions failed to fetch after retries - "
              f"listed in {failed_path} (accession + reason, tab-separated)")
        print(f"  Reasons: {dict(reason_counts)}")
        print(f"  Re-run this exact command to retry just these - already-"
              f"successful fetches are skipped automatically, not redone.")
    else:
        print("All sequences fetched successfully - no failures.")

    filtered.to_csv(out_table, sep="\t", index=False)
    print(f"Final: {n_written}/{len(filtered)} sequences written -> {out_fasta}")
