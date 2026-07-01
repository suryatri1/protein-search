#!/usr/bin/env python3
"""
Protein Presence/Absence Survey
================================
Given query protein(s) and a bacterial species, determine how many GTDB
strains carry each gene using DIAMOND blastp and (optionally) HMMER.

Usage:
    python protein_survey.py --query my_proteins.faa --species "Phocaeicola vulgatus" --outdir results/
    python protein_survey.py --query my_proteins.faa --species "Phocaeicola vulgatus" --pfam-db Pfam-A.hmm
    python protein_survey.py --query my_proteins.faa --species "Phocaeicola vulgatus" --hmm-profile custom.hmm

Thresholds default to Price et al. 2024 (fast.genomics): >=30% identity,
>=50% query coverage, E <= 1e-3.
"""

import argparse
import csv
import gzip
import json
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile
from collections import defaultdict
from pathlib import Path

__version__ = "0.1.0"

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYXBZJUOacdefghiklmnpqrstvwyxbzjuo*")


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def log(msg, level="INFO"):
    ts = time.strftime("%H:%M:%S")
    print(f"[{ts}] {level}: {msg}", file=sys.stderr, flush=True)


def run_cmd(cmd, desc=None, check=True):
    """Run a shell command, stream stderr, return CompletedProcess."""
    if desc:
        log(desc)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if check and result.returncode != 0:
        sys.exit(f"ERROR: {' '.join(cmd[:3])}... failed:\n{result.stderr.strip()}")
    return result


def check_tool(name, flag="--version"):
    path = shutil.which(name)
    if not path:
        return None
    try:
        r = subprocess.run([path, flag], capture_output=True, text=True, timeout=10)
        version = (r.stdout.strip() or r.stderr.strip()).split("\n")[0]
    except Exception:
        version = "unknown"
    return version


def parse_fasta(path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header, seq_parts = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    if header is not None:
        yield header, "".join(seq_parts)


def validate_query(path):
    """Parse and validate query FASTA. Returns list of (id, length) tuples."""
    queries = []
    for header, seq in parse_fasta(path):
        qid = header.split()[0]
        non_aa = set(seq) - AMINO_ACIDS
        if non_aa:
            sys.exit(f"ERROR: Non-amino-acid characters in {qid}: {non_aa}")
        if len(seq) < 10:
            log(f"WARNING: {qid} is only {len(seq)} aa — very short query", "WARN")
        queries.append((qid, len(seq)))
    if not queries:
        sys.exit("ERROR: No sequences found in query FASTA")
    log(f"Query: {len(queries)} protein(s), lengths {min(l for _,l in queries)}-{max(l for _,l in queries)} aa")
    return queries


# ---------------------------------------------------------------------------
# Step 1: GTDB metadata filtering
# ---------------------------------------------------------------------------

def normalize_species(name):
    """Convert user input to GTDB s__ format."""
    name = name.strip().strip('"').strip("'")
    name = re.sub(r"\s+", " ", name)
    if not name.startswith("s__"):
        name = f"s__{name}"
    return name


def filter_gtdb_metadata(metadata_path, species, max_genomes=None, include_genbank=False):
    """
    Stream-parse GTDB metadata TSV, filter by species.
    Returns list of dicts with genome info.
    """
    species_tag = normalize_species(species)
    log(f"Filtering GTDB metadata for {species_tag}")

    opener = gzip.open if str(metadata_path).endswith(".gz") else open
    genomes = []
    col_map = {}

    with opener(metadata_path, "rt") as f:
        for i, line in enumerate(f):
            fields = line.rstrip("\n").split("\t")
            if i == 0:
                for j, col in enumerate(fields):
                    col_map[col] = j
                required = ["accession", "gtdb_taxonomy"]
                for r in required:
                    if r not in col_map:
                        sys.exit(f"ERROR: GTDB metadata missing column '{r}'")
                continue

            taxonomy = fields[col_map["gtdb_taxonomy"]]
            if species_tag not in taxonomy:
                continue

            accession = fields[col_map["accession"]]
            prefix = ""
            if accession.startswith("RS_") or accession.startswith("GB_"):
                prefix = accession[:3]
                bare = accession[3:]
            else:
                bare = accession

            if not include_genbank and prefix == "GB_":
                continue

            def safe_get(col, default=""):
                return fields[col_map[col]] if col in col_map and col_map[col] < len(fields) else default

            genomes.append({
                "accession": bare,
                "gtdb_accession": accession,
                "source": "RefSeq" if prefix == "RS_" else "GenBank",
                "assembly_level": safe_get("ncbi_assembly_level", "unknown"),
                "checkm2_completeness": safe_get("checkm2_completeness", ""),
                "checkm2_contamination": safe_get("checkm2_contamination", ""),
                "protein_count": safe_get("protein_count", ""),
                "genome_size": safe_get("genome_size", ""),
            })

    if not genomes:
        log(f"No genomes found for {species_tag}. Trying without s__ prefix...", "WARN")
        alt = species_tag.replace("s__", "")
        parts = alt.split("_")
        if len(parts) >= 2:
            suggestions = []
            with opener(metadata_path, "rt") as f:
                for i, line in enumerate(f):
                    if i == 0:
                        continue
                    fields = line.rstrip("\n").split("\t")
                    taxonomy = fields[col_map["gtdb_taxonomy"]].lower()
                    if parts[0].lower() in taxonomy and parts[1].lower() in taxonomy:
                        match = re.search(r"s__(\S+)", fields[col_map["gtdb_taxonomy"]])
                        if match and match.group(0) not in suggestions:
                            suggestions.append(match.group(0))
            if suggestions:
                log(f"Did you mean one of: {', '.join(suggestions[:5])}?", "WARN")
        sys.exit(f"ERROR: No genomes found for species '{species}'")

    log(f"Found {len(genomes)} genomes ({sum(1 for g in genomes if g['source']=='RefSeq')} RefSeq, "
        f"{sum(1 for g in genomes if g['source']=='GenBank')} GenBank)")

    if max_genomes and len(genomes) > max_genomes:
        log(f"Sampling {max_genomes} of {len(genomes)} genomes (--max-genomes)")
        random.seed(42)
        genomes = random.sample(genomes, max_genomes)

    return genomes


# ---------------------------------------------------------------------------
# Step 2: Fetch proteomes from NCBI
# ---------------------------------------------------------------------------

def fetch_proteomes(genomes, outdir, threads=1):
    """Download protein FASTAs from NCBI using datasets CLI."""
    dl_dir = Path(outdir) / "ncbi_download"
    dl_dir.mkdir(parents=True, exist_ok=True)

    accession_file = dl_dir / "accessions.txt"
    with open(accession_file, "w") as f:
        for g in genomes:
            f.write(g["accession"] + "\n")

    zip_path = dl_dir / "proteomes.zip"

    cmd = [
        "datasets", "download", "genome", "accession",
        "--inputfile", str(accession_file),
        "--include", "protein",
        "--filename", str(zip_path),
    ]
    run_cmd(cmd, desc=f"Downloading proteomes for {len(genomes)} genomes from NCBI")

    log("Extracting proteomes")
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(dl_dir)

    data_dir = dl_dir / "ncbi_dataset" / "data"
    found = 0
    missing = []
    for g in genomes:
        faa = data_dir / g["accession"] / "protein.faa"
        if faa.exists():
            g["faa_path"] = str(faa)
            found += 1
        else:
            g["faa_path"] = None
            missing.append(g["accession"])

    log(f"Proteomes: {found}/{len(genomes)} downloaded, {len(missing)} missing")
    if missing:
        log(f"Missing accessions (first 10): {', '.join(missing[:10])}", "WARN")

    return genomes


# ---------------------------------------------------------------------------
# Step 3: Build unified target
# ---------------------------------------------------------------------------

def build_target_db(genomes, outdir, threads=8):
    """Concatenate proteomes with {accession}|{locus} headers, build DIAMOND DB."""
    target_dir = Path(outdir) / "target_db"
    target_dir.mkdir(parents=True, exist_ok=True)
    all_faa = target_dir / "all_proteins.faa"
    manifest_path = Path(outdir) / "manifest.tsv"

    total_proteins = 0
    with open(all_faa, "w") as out, open(manifest_path, "w", newline="") as mf:
        writer = csv.writer(mf, delimiter="\t")
        writer.writerow([
            "accession", "source", "assembly_level",
            "checkm2_completeness", "checkm2_contamination", "protein_count",
        ])

        for g in genomes:
            if not g["faa_path"]:
                writer.writerow([
                    g["accession"], g["source"], g["assembly_level"],
                    g["checkm2_completeness"], g["checkm2_contamination"], 0,
                ])
                continue

            pcount = 0
            for header, seq in parse_fasta(g["faa_path"]):
                locus = header.split()[0]
                desc = header[len(locus):].strip()
                out.write(f">{g['accession']}|{locus} {desc}\n{seq}\n")
                pcount += 1

            total_proteins += pcount
            writer.writerow([
                g["accession"], g["source"], g["assembly_level"],
                g["checkm2_completeness"], g["checkm2_contamination"], pcount,
            ])

    log(f"Concatenated {total_proteins:,} proteins from {sum(1 for g in genomes if g['faa_path'])} genomes")

    dmnd_path = target_dir / "target_db"
    run_cmd(
        ["diamond", "makedb", "--in", str(all_faa), "--db", str(dmnd_path),
         "--threads", str(threads)],
        desc="Building DIAMOND database",
    )

    return str(all_faa), str(dmnd_path), str(manifest_path)


# ---------------------------------------------------------------------------
# Step 4: DIAMOND blastp
# ---------------------------------------------------------------------------

def run_diamond(query_faa, db_path, outdir, evalue, threads):
    """Run DIAMOND blastp, return path to raw results."""
    blast_out = Path(outdir) / "blast_raw.tsv"
    cmd = [
        "diamond", "blastp",
        "--query", query_faa,
        "--db", db_path,
        "--out", str(blast_out),
        "--outfmt", "6",
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "qlen", "slen", "qcovhsp",
        "--evalue", str(evalue),
        "--max-target-seqs", "0",
        "--sensitive",
        "--threads", str(threads),
    ]
    run_cmd(cmd, desc="Running DIAMOND blastp")
    return str(blast_out)


def parse_blast_results(blast_tsv, min_pident, min_qcov):
    """Parse DIAMOND output, apply thresholds, return best hit per (query, accession)."""
    hits = defaultdict(list)
    total = 0
    passed = 0

    with open(blast_tsv) as f:
        for line in f:
            total += 1
            fields = line.rstrip("\n").split("\t")
            qseqid = fields[0]
            sseqid = fields[1]
            pident = float(fields[2])
            length = int(fields[3])
            evalue = float(fields[10])
            bitscore = float(fields[11])
            qlen = int(fields[12])
            slen = int(fields[13])
            qcovhsp = float(fields[14])

            accession = sseqid.split("|")[0]
            locus = sseqid.split("|")[1].split()[0] if "|" in sseqid else sseqid

            if pident < min_pident or qcovhsp < min_qcov:
                continue
            passed += 1

            hits[(qseqid, accession)].append({
                "query_id": qseqid,
                "accession": accession,
                "hit_locus": locus,
                "pident": pident,
                "evalue": evalue,
                "bitscore": bitscore,
                "qcov": qcovhsp,
                "length": length,
                "qlen": qlen,
                "slen": slen,
            })

    best_hits = {}
    for key, hit_list in hits.items():
        best_hits[key] = max(hit_list, key=lambda h: h["bitscore"])

    log(f"BLAST: {total} raw hits → {passed} passed thresholds → {len(best_hits)} best hits (query×genome)")
    return best_hits


# ---------------------------------------------------------------------------
# Step 5: HMM search
# ---------------------------------------------------------------------------

def run_hmmscan(query_faa, pfam_db, outdir, threads):
    """hmmscan query against Pfam, return dict of query_id -> best Pfam hit."""
    domtbl = Path(outdir) / "pfam_scan.domtbl"
    cmd = [
        "hmmscan", "--cpu", str(threads),
        "--domtblout", str(domtbl),
        "--noali",
        str(pfam_db), query_faa,
    ]
    run_cmd(cmd, desc="Scanning query proteins against Pfam")
    return parse_domtblout(str(domtbl), key_field="query")


def parse_domtblout(domtbl_path, key_field="query"):
    """
    Parse HMMER domtblout format.
    key_field="query": group by query name (for hmmscan results)
    key_field="target": group by target name (for hmmsearch results)
    """
    hits = defaultdict(list)
    with open(domtbl_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 22:
                continue
            target_name = parts[0]
            target_acc = parts[1]
            query_name = parts[3]
            query_acc = parts[4]
            full_evalue = float(parts[6])
            full_score = float(parts[7])
            dom_evalue = float(parts[12])
            dom_score = float(parts[13])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])
            env_from = int(parts[19])
            env_to = int(parts[20])
            tlen = int(parts[2])
            qlen = int(parts[5])

            hit = {
                "target_name": target_name,
                "target_acc": target_acc if target_acc != "-" else target_name,
                "query_name": query_name,
                "full_evalue": full_evalue,
                "full_score": full_score,
                "dom_evalue": dom_evalue,
                "dom_score": dom_score,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "hmm_coverage": (hmm_to - hmm_from + 1) / tlen if key_field == "query" else 0,
                "ali_from": ali_from,
                "ali_to": ali_to,
            }

            if key_field == "query":
                hits[query_name].append(hit)
            else:
                hits[target_name].append(hit)

    return dict(hits)


def run_hmmsearch(pfam_db, pfam_acc, all_proteins_faa, outdir, threads):
    """Extract a Pfam profile and search against target proteomes."""
    hmm_dir = Path(outdir) / "hmm_profiles"
    hmm_dir.mkdir(exist_ok=True)

    profile_path = hmm_dir / f"{pfam_acc}.hmm"
    run_cmd(
        ["hmmfetch", "-o", str(profile_path), str(pfam_db), pfam_acc],
        desc=f"Extracting HMM profile {pfam_acc}",
    )

    if not profile_path.exists() or profile_path.stat().st_size == 0:
        log(f"Could not extract profile for {pfam_acc} — skipping", "WARN")
        return {}

    domtbl = Path(outdir) / f"hmmsearch_{pfam_acc}.domtbl"
    run_cmd(
        ["hmmsearch", "--cpu", str(threads),
         "--domtblout", str(domtbl),
         "--noali",
         str(profile_path), all_proteins_faa],
        desc=f"Searching {pfam_acc} against target proteomes",
    )

    return parse_domtblout(str(domtbl), key_field="target")


def do_hmm_step(query_faa, queries, pfam_db, hmm_profile, all_proteins_faa, outdir, threads):
    """
    Orchestrate the HMM search step with multi-domain architecture matching.
    Returns (hmm_results, domain_fingerprints, per_domain_counts) where:
      hmm_results: query_id -> {accession -> hit_dict}
      domain_fingerprints: query_id -> list of {pfam_acc, name, evalue, coverage}
      per_domain_counts: query_id -> {pfam_acc -> n_genomes_with_domain}
    """
    hmm_results = {}
    domain_fingerprints = {}
    per_domain_counts = {}

    if hmm_profile:
        log(f"Using user-supplied HMM profile: {hmm_profile}")
        domtbl = Path(outdir) / "hmmsearch_custom.domtbl"
        run_cmd(
            ["hmmsearch", "--cpu", str(threads),
             "--domtblout", str(domtbl),
             "--noali",
             hmm_profile, all_proteins_faa],
            desc="Searching custom HMM profile against targets",
        )
        target_hits = parse_domtblout(str(domtbl), key_field="target")
        for qid, _ in queries:
            hmm_results[qid] = _group_hmm_by_accession_single(target_hits, "custom_profile")
        return hmm_results, {}, {}

    if not pfam_db:
        log("No Pfam database provided — skipping HMM search")
        return {}, {}, {}

    pfam_hits = run_hmmscan(query_faa, pfam_db, outdir, threads)

    for qid, _ in queries:
        qhits = pfam_hits.get(qid, [])
        significant = [h for h in qhits if h["dom_evalue"] < 1e-3 and h["hmm_coverage"] > 0.5]

        if not significant:
            log(f"  {qid}: no significant Pfam domain found — BLAST-only")
            hmm_results[qid] = None
            domain_fingerprints[qid] = []
            per_domain_counts[qid] = {}
            continue

        domains_by_acc = {}
        for h in significant:
            acc = h["target_acc"]
            if acc not in domains_by_acc or h["dom_score"] > domains_by_acc[acc]["dom_score"]:
                domains_by_acc[acc] = h

        fingerprint = sorted(domains_by_acc.values(), key=lambda h: h["ali_from"])
        domain_fingerprints[qid] = [
            {"pfam_acc": d["target_acc"], "name": d["target_name"],
             "evalue": d["dom_evalue"], "coverage": d["hmm_coverage"]}
            for d in fingerprint
        ]

        domain_labels = ", ".join(
            f"{d['target_acc']} ({d['target_name']})" for d in fingerprint
        )
        log(f"  {qid}: {len(fingerprint)} Pfam domains: {domain_labels}")

        all_domain_hits = {}
        domain_genome_counts = {}
        for d in fingerprint:
            pfam_acc = d["target_acc"]
            target_hits = run_hmmsearch(pfam_db, pfam_acc, all_proteins_faa, outdir, threads)

            locus_hits = {}
            genome_set = set()
            for target_name, hit_list in target_hits.items():
                if "|" not in target_name:
                    continue
                accession = target_name.split("|")[0]
                locus = target_name.split("|")[1].split()[0]
                genome_set.add(accession)
                best_hit = max(hit_list, key=lambda h: h["dom_score"])
                key = f"{accession}|{locus}"
                if key not in locus_hits or best_hit["dom_score"] > locus_hits[key]["dom_score"]:
                    locus_hits[key] = best_hit

            all_domain_hits[pfam_acc] = locus_hits
            domain_genome_counts[pfam_acc] = len(genome_set)
            log(f"    {pfam_acc}: found in {len(genome_set)} genomes")

        per_domain_counts[qid] = domain_genome_counts

        hmm_results[qid] = _build_architecture_results(
            all_domain_hits, fingerprint
        )

    return hmm_results, domain_fingerprints, per_domain_counts


def _build_architecture_results(all_domain_hits, fingerprint):
    """
    For each genome, find the ORF with the most domain matches.
    Require ALL domains in the same ORF for architecture match.
    """
    required_domains = {d["target_acc"] for d in fingerprint}
    n_required = len(required_domains)

    locus_domain_map = defaultdict(dict)
    for pfam_acc, locus_hits in all_domain_hits.items():
        for key, hit in locus_hits.items():
            locus_domain_map[key][pfam_acc] = hit

    by_accession = defaultdict(list)
    for key, domains_found in locus_domain_map.items():
        accession, locus = key.split("|", 1)
        matched_set = set(domains_found.keys()) & required_domains
        n_matched = len(matched_set)

        scores = [domains_found[d]["dom_score"] for d in matched_set]
        evalues = [domains_found[d]["dom_evalue"] for d in matched_set]

        by_accession[accession].append({
            "locus": locus,
            "domains_matched": matched_set,
            "n_matched": n_matched,
            "total_score": sum(scores),
            "worst_evalue": max(evalues),
            "best_evalue": min(evalues),
        })

    results = {}
    for acc, orf_list in by_accession.items():
        best_orf = max(orf_list, key=lambda o: (o["n_matched"], o["total_score"]))
        is_arch_match = best_orf["n_matched"] == n_required

        results[acc] = {
            "accession": acc,
            "hit_locus": best_orf["locus"],
            "hmm_domains_matched": best_orf["n_matched"],
            "hmm_domains_required": n_required,
            "hmm_architecture_match": is_arch_match,
            "hmm_domain": ",".join(sorted(best_orf["domains_matched"])),
            "hmm_evalue": best_orf["worst_evalue"],
            "hmm_score": best_orf["total_score"],
        }

    return results


def _group_hmm_by_accession_single(target_hits, domain_name):
    """Group hmmsearch hits by genome accession (single-domain / custom HMM path)."""
    by_accession = defaultdict(list)
    for target_name, hit_list in target_hits.items():
        if "|" not in target_name:
            continue
        accession = target_name.split("|")[0]
        locus = target_name.split("|")[1].split()[0]
        for h in hit_list:
            by_accession[accession].append({
                "accession": accession,
                "hit_locus": locus,
                "hmm_domain": domain_name,
                "hmm_evalue": h["dom_evalue"],
                "hmm_score": h["dom_score"],
                "hmm_domains_matched": 1,
                "hmm_domains_required": 1,
                "hmm_architecture_match": True,
            })

    best = {}
    for acc, hits in by_accession.items():
        best[acc] = max(hits, key=lambda h: h["hmm_score"])
    return best


def _load_genomes_from_manifest(manifest_path, target_dir):
    """Reconstruct genome list from manifest.tsv for --reuse-target mode."""
    genomes = []
    with open(manifest_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc = row["accession"]
            faa = target_dir.parent / "ncbi_download" / "ncbi_dataset" / "data" / acc / "protein.faa"
            genomes.append({
                "accession": acc,
                "gtdb_accession": acc,
                "source": row.get("source", "RefSeq"),
                "assembly_level": row.get("assembly_level", "unknown"),
                "checkm2_completeness": row.get("checkm2_completeness", ""),
                "checkm2_contamination": row.get("checkm2_contamination", ""),
                "protein_count": row.get("protein_count", ""),
                "genome_size": "",
                "faa_path": str(faa) if faa.exists() else None,
            })
    return genomes


# ---------------------------------------------------------------------------
# Step 6: Merge & report
# ---------------------------------------------------------------------------

def merge_and_report(queries, genomes, blast_hits, hmm_results,
                     domain_fingerprints, per_domain_counts, outdir, args):
    """Produce presence_absence.tsv, summary.txt, and run_log.json."""
    outdir = Path(outdir)

    pa_path = outdir / "presence_absence.tsv"
    with open(pa_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        header = [
            "query_id", "accession", "source", "assembly_level", "completeness",
            "blast_hit", "blast_pident", "blast_evalue", "blast_bitscore",
            "blast_qcov", "blast_locus",
            "hmm_hit", "hmm_evalue", "hmm_score", "hmm_domain",
            "hmm_domains_matched", "hmm_arch_match",
            "call",
        ]
        w.writerow(header)

        for qid, qlen in queries:
            for g in genomes:
                acc = g["accession"]
                if not g["faa_path"]:
                    continue

                bkey = (qid, acc)
                bhit = blast_hits.get(bkey)

                hhit = None
                if qid in hmm_results and hmm_results[qid] is not None:
                    hhit = hmm_results[qid].get(acc)

                blast_present = bhit is not None
                hmm_arch_match = (hhit is not None
                                  and hhit.get("hmm_architecture_match", False))

                completeness = g.get("checkm2_completeness", "")
                try:
                    comp_val = float(completeness)
                except (ValueError, TypeError):
                    comp_val = None

                if blast_present or hmm_arch_match:
                    call = "present"
                elif comp_val is not None and comp_val < 90:
                    call = "absent*"
                else:
                    call = "absent"

                n_matched = hhit["hmm_domains_matched"] if hhit else 0
                n_required = hhit["hmm_domains_required"] if hhit else 0

                w.writerow([
                    qid, acc, g["source"], g["assembly_level"], completeness,
                    "yes" if blast_present else "no",
                    f"{bhit['pident']:.1f}" if bhit else "-",
                    f"{bhit['evalue']:.1e}" if bhit else "-",
                    f"{bhit['bitscore']:.1f}" if bhit else "-",
                    f"{bhit['qcov']:.1f}" if bhit else "-",
                    bhit["hit_locus"] if bhit else "-",
                    "yes" if hhit else "no",
                    f"{hhit['hmm_evalue']:.1e}" if hhit else "-",
                    f"{hhit['hmm_score']:.1f}" if hhit else "-",
                    hhit["hmm_domain"] if hhit else "-",
                    f"{n_matched}/{n_required}" if hhit else "-",
                    "yes" if hmm_arch_match else "no",
                    call,
                ])

    log(f"Wrote {pa_path}")

    write_summary(queries, genomes, blast_hits, hmm_results,
                  domain_fingerprints, per_domain_counts, outdir, args)
    write_run_log(outdir, args, domain_fingerprints)


def write_summary(queries, genomes, blast_hits, hmm_results,
                  domain_fingerprints, per_domain_counts, outdir, args):
    """Write human-readable summary.txt."""
    summary_path = outdir / "summary.txt"
    searchable = [g for g in genomes if g["faa_path"]]
    n_searched = len(searchable)

    level_counts = defaultdict(int)
    for g in searchable:
        level_counts[g["assembly_level"]] += 1

    with open(summary_path, "w") as f:
        f.write(f"Protein Presence/Absence Survey\n")
        f.write(f"{'=' * 40}\n")
        f.write(f"Species: {args.species}\n")
        f.write(f"Genomes in GTDB: {len(genomes)}\n")
        f.write(f"Proteomes searched: {n_searched} "
                f"({len(genomes) - n_searched} skipped — no protein annotation)\n")
        f.write(f"Thresholds: pident >= {args.min_pident}%, qcov >= {args.min_qcov}%, "
                f"E-value <= {args.evalue}\n\n")

        for qid, qlen in queries:
            f.write(f"=== {qid} ({qlen} aa) ===\n\n")

            blast_present = []
            blast_absent = []
            low_comp_absent = 0
            pidents = []
            qcovs = []

            for g in searchable:
                bkey = (qid, g["accession"])
                if bkey in blast_hits:
                    blast_present.append(g)
                    pidents.append(blast_hits[bkey]["pident"])
                    qcovs.append(blast_hits[bkey]["qcov"])
                else:
                    blast_absent.append(g)
                    try:
                        if float(g.get("checkm2_completeness", 100)) < 90:
                            low_comp_absent += 1
                    except (ValueError, TypeError):
                        pass

            np_blast = len(blast_present)
            f.write(f"DIAMOND blastp:\n")
            f.write(f"  Present: {np_blast}/{n_searched} ({100*np_blast/n_searched:.1f}%)\n")
            f.write(f"  Absent:  {n_searched-np_blast}/{n_searched} ({100*(n_searched-np_blast)/n_searched:.1f}%)\n")
            if low_comp_absent:
                f.write(f"    Of which {low_comp_absent} are low-completeness (<90%) "
                        f"— true absence uncertain\n")
            if pidents:
                mean_pi = sum(pidents) / len(pidents)
                std_pi = (sum((x - mean_pi) ** 2 for x in pidents) / len(pidents)) ** 0.5
                mean_qc = sum(qcovs) / len(qcovs)
                f.write(f"  Pident (hits): mean {mean_pi:.1f}% +/- {std_pi:.1f}  |  "
                        f"range {min(pidents):.1f}-{max(pidents):.1f}%\n")
                f.write(f"  Qcov (hits):   mean {mean_qc:.1f}%\n")
            f.write("\n")

            if qid in hmm_results:
                qhmm = hmm_results[qid]
                if qhmm is None:
                    f.write(f"HMM: no Pfam domain found — skipped\n\n")
                else:
                    fp = domain_fingerprints.get(qid, [])
                    dc = per_domain_counts.get(qid, {})

                    if fp:
                        f.write(f"Domain fingerprint ({len(fp)} Pfam domains):\n")
                        for d in fp:
                            n_with = dc.get(d["pfam_acc"], 0)
                            f.write(f"  {d['pfam_acc']} ({d['name']}): "
                                    f"E={d['evalue']:.1e}, cov={d['coverage']:.0%} "
                                    f"— found in {n_with}/{n_searched} genomes "
                                    f"({100*n_with/n_searched:.0f}%)\n")
                        f.write("\n")

                    all_accs = {g["accession"] for g in searchable}
                    blast_present_set = {g["accession"] for g in blast_present}

                    arch_match_set = {
                        acc for acc, hit in qhmm.items()
                        if hit.get("hmm_architecture_match", False) and acc in all_accs
                    }
                    any_domain_set = {
                        acc for acc in qhmm if acc in all_accs
                    }

                    n_arch = len(arch_match_set)
                    n_any = len(any_domain_set)

                    f.write(f"HMM domain architecture (all {len(fp)} domains in same ORF):\n")
                    f.write(f"  Any domain hit:    {n_any}/{n_searched} ({100*n_any/n_searched:.1f}%)\n")
                    f.write(f"  Architecture match: {n_arch}/{n_searched} ({100*n_arch/n_searched:.1f}%)\n\n")

                    both = len(blast_present_set & arch_match_set)
                    blast_only = len(blast_present_set - arch_match_set)
                    hmm_arch_only = len(arch_match_set - blast_present_set)
                    neither = n_searched - len(blast_present_set | arch_match_set)

                    f.write(f"Concordance (BLAST vs HMM architecture):\n")
                    f.write(f"  BLAST+ HMM-arch+: {both}\n")
                    f.write(f"  BLAST+ HMM-arch-: {blast_only}\n")
                    f.write(f"  BLAST- HMM-arch+: {hmm_arch_only}\n")
                    f.write(f"  BLAST- HMM-arch-: {neither}\n\n")

                    if hmm_arch_only > 0.25 * n_searched:
                        warn_msg = (
                            f"WARNING: HMM architecture matches "
                            f"{100*hmm_arch_only/n_searched:.0f}% of genomes "
                            f"without BLAST support — domain architecture may "
                            f"overlap with a different protein family. "
                            f"BLAST-only calling recommended."
                        )
                        f.write(f"  {warn_msg}\n\n")
                        log(warn_msg, "WARN")

            f.write(f"By assembly level:\n")
            for level in sorted(level_counts.keys()):
                level_genomes = [g for g in searchable if g["assembly_level"] == level]
                level_present = sum(1 for g in level_genomes if (qid, g["accession"]) in blast_hits)
                n_level = len(level_genomes)
                pct = 100 * level_present / n_level if n_level else 0
                note = "  <- lower rate may reflect assembly gaps" if level == "Contig" and pct < 85 else ""
                f.write(f"  {level} ({n_level}): {level_present}/{n_level} ({pct:.1f}%){note}\n")
            f.write("\n")

    log(f"Wrote {summary_path}")


def write_run_log(outdir, args, domain_fingerprints=None):
    """Write run parameters and tool versions as JSON."""
    log_path = outdir / "run_log.json"
    run_log = {
        "tool": "protein_survey",
        "version": __version__,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "args": {
            "species": args.species,
            "evalue": args.evalue,
            "min_pident": args.min_pident,
            "min_qcov": args.min_qcov,
            "threads": args.threads,
            "max_genomes": args.max_genomes,
            "include_genbank": args.include_genbank,
            "pfam_db": args.pfam_db,
            "hmm_profile": args.hmm_profile,
        },
        "tool_versions": {
            "diamond": check_tool("diamond") or "not found",
            "hmmscan": check_tool("hmmscan", "-h") or "not found",
            "datasets": check_tool("datasets") or "not found",
        },
    }
    if domain_fingerprints:
        run_log["hmm_domains"] = domain_fingerprints
    with open(log_path, "w") as f:
        json.dump(run_log, f, indent=2)
    log(f"Wrote {log_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Protein presence/absence survey across GTDB strains",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--query", required=True,
                        help="Query protein FASTA (1 or more sequences)")
    parser.add_argument("--species", required=True,
                        help="Species name (e.g. 'Phocaeicola vulgatus')")
    parser.add_argument("--outdir", default="protein_survey_results",
                        help="Output directory (default: protein_survey_results)")
    parser.add_argument("--threads", type=int, default=8,
                        help="CPU threads for DIAMOND/HMMER (default: 8)")
    parser.add_argument("--evalue", type=float, default=1e-3,
                        help="E-value threshold (default: 1e-3, per Price et al. 2024)")
    parser.add_argument("--min-pident", type=float, default=30.0,
                        help="Minimum percent identity (default: 30)")
    parser.add_argument("--min-qcov", type=float, default=50.0,
                        help="Minimum query coverage %% (default: 50)")
    parser.add_argument("--pfam-db",
                        help="Path to Pfam-A.hmm (enables HMM search)")
    parser.add_argument("--hmm-profile",
                        help="Path to a custom HMM profile (bypasses Pfam scan)")
    parser.add_argument("--gtdb-metadata", required=True,
                        help="Path to GTDB bac120_metadata_r*.tsv.gz")
    parser.add_argument("--max-genomes", type=int,
                        help="Cap number of genomes (random sample, seed=42)")
    parser.add_argument("--include-genbank", action="store_true",
                        help="Include GenBank-only genomes (default: RefSeq only)")
    parser.add_argument("--prepare-only", action="store_true",
                        help="Download proteomes and build target DB, then exit")
    parser.add_argument("--reuse-target",
                        help="Path to pre-built target dir (skip download/build steps)")

    args = parser.parse_args()

    # Validate inputs
    if not Path(args.query).exists():
        sys.exit(f"ERROR: Query file not found: {args.query}")
    if not Path(args.gtdb_metadata).exists():
        sys.exit(f"ERROR: GTDB metadata not found: {args.gtdb_metadata}")
    if args.pfam_db and not Path(args.pfam_db).exists():
        sys.exit(f"ERROR: Pfam database not found: {args.pfam_db}")
    if args.hmm_profile and not Path(args.hmm_profile).exists():
        sys.exit(f"ERROR: HMM profile not found: {args.hmm_profile}")

    # Check required tools
    for tool in ["diamond", "datasets"]:
        if not shutil.which(tool):
            sys.exit(f"ERROR: '{tool}' not found in PATH. Install it first.")
    if (args.pfam_db or args.hmm_profile) and not shutil.which("hmmscan"):
        sys.exit("ERROR: 'hmmscan' not found in PATH. Install HMMER or remove --pfam-db/--hmm-profile.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    t_start = time.time()

    # Step 0: Validate query
    queries = validate_query(args.query)

    if args.reuse_target:
        target_dir = Path(args.reuse_target)
        all_faa = str(target_dir / "all_proteins.faa")
        db_path = str(target_dir / "target_db")
        manifest_path = str(target_dir.parent / "manifest.tsv")

        if not Path(all_faa).exists():
            sys.exit(f"ERROR: {all_faa} not found in --reuse-target dir")
        if not Path(db_path + ".dmnd").exists():
            sys.exit(f"ERROR: {db_path}.dmnd not found in --reuse-target dir")

        genomes = _load_genomes_from_manifest(manifest_path, target_dir)
        log(f"Reusing target DB from {target_dir} ({len(genomes)} genomes)")
    else:
        # Step 1: Filter GTDB metadata
        genomes = filter_gtdb_metadata(
            args.gtdb_metadata, args.species,
            max_genomes=args.max_genomes,
            include_genbank=args.include_genbank,
        )

        # Step 2: Fetch proteomes
        genomes = fetch_proteomes(genomes, str(outdir), threads=args.threads)

        # Step 3: Build unified target DB
        all_faa, db_path, manifest_path = build_target_db(genomes, str(outdir), threads=args.threads)

        if args.prepare_only:
            elapsed = time.time() - t_start
            log(f"Prepare-only done in {elapsed:.0f}s. Target DB in {outdir}/target_db/")
            return

    # Step 4: DIAMOND blastp
    blast_tsv = run_diamond(args.query, db_path, str(outdir), args.evalue, args.threads)
    blast_hits = parse_blast_results(blast_tsv, args.min_pident, args.min_qcov)

    # Step 5: HMM search (conditional)
    hmm_results, domain_fingerprints, per_domain_counts = do_hmm_step(
        args.query, queries,
        args.pfam_db, args.hmm_profile,
        all_faa, str(outdir), args.threads,
    )

    # Step 6: Merge & report
    merge_and_report(queries, genomes, blast_hits, hmm_results,
                     domain_fingerprints, per_domain_counts, outdir, args)

    elapsed = time.time() - t_start
    log(f"Done in {elapsed:.0f}s. Results in {outdir}/")


if __name__ == "__main__":
    main()
