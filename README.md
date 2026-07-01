# Protein Presence/Absence Survey

Search for query proteins across all GTDB strains of a bacterial species.
Reports per-genome presence/absence using DIAMOND blastp and (optionally) HMMER domain architecture matching.

## Setup

```bash
conda env create -f environment.yml
conda activate protein_survey
```

### Required databases (download once)

**GTDB metadata** (~500 MB compressed):
```bash
curl -fSLO https://data.gtdb.ecogenomic.org/releases/release220/220.0/bac120_metadata_r220.tsv.gz
```

**Pfam-A** (optional, ~1.4 GB compressed — enables HMM search):
```bash
curl -fSLO https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

## Quick start (Snakemake)

1. Place query protein FASTAs (one per file) in `query_proteins/`
2. Edit `config.yaml` — set species, database paths, thresholds
3. Run:

```bash
snakemake --cores 4
```

The workflow downloads proteomes once, then runs BLAST + HMM per query protein.

### Snakemake output

```
results/
├── target/                        # Shared across queries (built once)
│   ├── target_db/all_proteins.faa
│   ├── target_db/target_db.dmnd
│   └── manifest.tsv
├── pucA/                          # Per-query results
│   ├── presence_absence.tsv
│   ├── summary.txt
│   └── run_log.json
├── pucB/
│   └── ...
└── combined_presence_absence.tsv  # All queries merged
```

## Quick start (standalone)

```bash
python protein_survey.py \
  --query my_protein.faa \
  --species "Phocaeicola vulgatus" \
  --gtdb-metadata bac120_metadata_r220.tsv.gz \
  --pfam-db Pfam-A.hmm \
  --outdir results/
```

## How it works

1. **Filter GTDB metadata** for all strains of the target species (GTDB r220 taxonomy, 95% ANI species boundaries)
2. **Download proteomes** from NCBI using `datasets` CLI
3. **Build unified DIAMOND database** from all target proteomes
4. **DIAMOND blastp** query protein(s) against the combined database
5. **(Optional) HMM domain architecture matching**:
   - `hmmscan` query against Pfam to identify all significant domains (E < 1e-3, coverage > 50%)
   - `hmmsearch` each domain against all target proteomes
   - Require all query domains to co-occur in the same ORF for a positive call
   - Report concordance between BLAST and HMM results
6. **Report** per-genome presence/absence with metrics

## Default thresholds

Following Price et al. 2024 ([fast.genomics](https://doi.org/10.1371/journal.pone.0301871)):

| Parameter | Default | Flag |
|---|---|---|
| Percent identity | >= 30% | `--min-pident` |
| Query coverage | >= 50% | `--min-qcov` |
| E-value | <= 1e-3 | `--evalue` |

## Output columns (presence_absence.tsv)

| Column | Description |
|---|---|
| `query_id` | Query protein identifier |
| `accession` | NCBI assembly accession |
| `source` | RefSeq or GenBank |
| `assembly_level` | Chromosome / Complete Genome / Scaffold / Contig |
| `completeness` | CheckM2 completeness (%) |
| `blast_hit` | yes/no |
| `blast_pident` | Percent identity of best BLAST hit |
| `blast_evalue` | E-value of best BLAST hit |
| `blast_bitscore` | Bitscore of best BLAST hit |
| `blast_qcov` | Query coverage (%) |
| `blast_locus` | Locus tag of best BLAST hit |
| `hmm_hit` | yes/no (if HMM search was run) |
| `hmm_evalue` | Worst domain E-value among matched domains |
| `hmm_score` | Sum of domain scores in best ORF |
| `hmm_domain` | Matched Pfam domain(s) |
| `hmm_domains_matched` | N matched / N required |
| `hmm_arch_match` | yes/no — all domains in same ORF |
| `call` | `present`, `absent`, or `absent*` |

`absent*` = genome completeness < 90%; absence may reflect incomplete assembly rather than true gene loss.

## HMM domain architecture matching

When `--pfam-db` is provided, the tool identifies the complete Pfam domain fingerprint of each query protein, then checks whether the same domain architecture exists in each target genome. This catches distant homologs that BLAST may miss.

A concordance analysis compares BLAST and HMM results. When > 25% of genomes have HMM architecture matches without BLAST support, the tool warns that the domain architecture may overlap with a different protein family (e.g., TetQ shares all domains with EF-G). In such cases, BLAST-only calling is more reliable.

## CLI options

```
--query FASTA          Query protein FASTA (required)
--species NAME         Species name, e.g. "Phocaeicola vulgatus" (required)
--gtdb-metadata FILE   Path to bac120_metadata_r220.tsv.gz (required)
--outdir DIR           Output directory (default: protein_survey_results)
--threads N            CPU threads (default: 8)
--evalue FLOAT         E-value threshold (default: 1e-3)
--min-pident FLOAT     Min percent identity (default: 30)
--min-qcov FLOAT       Min query coverage % (default: 50)
--pfam-db FILE         Path to Pfam-A.hmm (enables HMM search)
--hmm-profile FILE     Custom HMM profile (bypasses Pfam scan)
--max-genomes N        Cap genome count (random sample, seed=42)
--include-genbank      Include GenBank-only genomes (default: RefSeq only)
--prepare-only         Download proteomes and build target DB, then exit
--reuse-target DIR     Path to pre-built target_db/ (skip download/build)
```

## Dependencies

| Tool | Version | Purpose |
|---|---|---|
| [DIAMOND](https://github.com/bbuchfink/diamond) | >= 2.0 | Protein search |
| [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) | >= 16 | Proteome download |
| [HMMER](http://hmmer.org/) | >= 3.3 | HMM domain search (optional) |
| [Snakemake](https://snakemake.readthedocs.io/) | >= 8.0 | Workflow orchestration (optional) |
| Python | >= 3.9 | No pip packages required (stdlib only) |

## Notes

- **Species names**: Accepts common formats ("Phocaeicola vulgatus", "P. vulgatus"). Suggests corrections if no match is found.
- **RefSeq vs GenBank**: Defaults to RefSeq-only for higher assembly quality. Use `--include-genbank` for broader strain coverage.
- **Large species**: For species with thousands of genomes (e.g., *E. coli*), use `--max-genomes` to sample.
- **Assembly completeness**: Contig-level assemblies may show lower gene presence due to assembly gaps. The `absent*` flag highlights this.

## Citation

Default thresholds are based on:

> Price MN, Arkin AP. A fast comparative genome browser for diverse bacteria and archaea. *PLoS ONE* 19(4): e0301871 (2024). https://doi.org/10.1371/journal.pone.0301871
