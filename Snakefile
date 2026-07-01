"""
Protein Presence/Absence Survey — Snakemake Workflow
=====================================================
Searches for query proteins across all GTDB strains of a bacterial species
using DIAMOND blastp and HMMER domain architecture matching.

Usage:
    snakemake --cores 4
    snakemake --cores 4 -n          # dry-run
    snakemake --cores 4 --forceall  # re-run everything
"""

from pathlib import Path

configfile: "config.yaml"

SPECIES = config["species"]
OUTDIR = config["outdir"]
QUERY_DIR = config["query_dir"]
THREADS = config.get("threads", 4)

QUERIES = glob_wildcards(f"{QUERY_DIR}/{{name}}.faa").name

if not QUERIES:
    raise ValueError(f"No .faa files found in {QUERY_DIR}/")


rule all:
    input:
        expand(f"{OUTDIR}/{{name}}/summary.txt", name=QUERIES),
        f"{OUTDIR}/combined_presence_absence.tsv",


rule prepare_target:
    """Download proteomes and build DIAMOND database (once per species)."""
    input:
        gtdb=config["gtdb_metadata"],
    output:
        faa=f"{OUTDIR}/target/target_db/all_proteins.faa",
        dmnd=f"{OUTDIR}/target/target_db/target_db.dmnd",
        manifest=f"{OUTDIR}/target/manifest.tsv",
    params:
        species=SPECIES,
        outdir=f"{OUTDIR}/target",
        genbank="--include-genbank" if config.get("include_genbank") else "",
        max_genomes=f"--max-genomes {config['max_genomes']}" if config.get("max_genomes") else "",
    threads: THREADS
    shell:
        """
        python protein_survey.py \
            --query {QUERY_DIR}/{QUERIES[0]}.faa \
            --species "{params.species}" \
            --gtdb-metadata {input.gtdb} \
            --outdir {params.outdir} \
            --threads {threads} \
            --prepare-only \
            {params.genbank} {params.max_genomes}
        """


rule search_protein:
    """Run BLAST + HMM search for each query protein."""
    input:
        query=f"{QUERY_DIR}/{{name}}.faa",
        faa=rules.prepare_target.output.faa,
        dmnd=rules.prepare_target.output.dmnd,
        manifest=rules.prepare_target.output.manifest,
    output:
        summary=f"{OUTDIR}/{{name}}/summary.txt",
        tsv=f"{OUTDIR}/{{name}}/presence_absence.tsv",
        log_json=f"{OUTDIR}/{{name}}/run_log.json",
    params:
        species=SPECIES,
        target_dir=f"{OUTDIR}/target/target_db",
        pfam=f"--pfam-db {config['pfam_db']}" if config.get("pfam_db") else "",
        evalue=config.get("evalue", 1e-3),
        min_pident=config.get("min_pident", 30.0),
        min_qcov=config.get("min_qcov", 50.0),
    threads: THREADS
    shell:
        """
        python protein_survey.py \
            --query {input.query} \
            --species "{params.species}" \
            --gtdb-metadata {config[gtdb_metadata]} \
            --reuse-target {params.target_dir} \
            --outdir {OUTDIR}/{wildcards.name} \
            --threads {threads} \
            --evalue {params.evalue} \
            --min-pident {params.min_pident} \
            --min-qcov {params.min_qcov} \
            {params.pfam}
        """


rule combine_results:
    """Merge all per-query presence_absence.tsv files into one."""
    input:
        expand(f"{OUTDIR}/{{name}}/presence_absence.tsv", name=QUERIES),
    output:
        f"{OUTDIR}/combined_presence_absence.tsv",
    run:
        import csv
        header_written = False
        with open(output[0], "w", newline="") as out:
            writer = None
            for tsv in input:
                with open(tsv) as f:
                    reader = csv.reader(f, delimiter="\t")
                    hdr = next(reader)
                    if not header_written:
                        writer = csv.writer(out, delimiter="\t")
                        writer.writerow(hdr)
                        header_written = True
                    for row in reader:
                        writer.writerow(row)
