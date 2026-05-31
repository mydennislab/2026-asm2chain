# asm2chain

Bidirectional UCSC liftOver chain files from two genome assemblies. Given
FASTA `A` and FASTA `B` (same or closely related species, any naming
convention, gzipped or plain), `asm2chain` aligns them with `minimap2` in
both directions, converts each PAF to a UCSC chain with `transanno`,
sorts and gzips with UCSC `chainSort`, and emits per-target-sequence
chain statistics for each direction. The output is suitable for `liftOver`
in either direction.

## What you get

For configured labels `A` and `B`:

| Path                                  | Description                                                       |
| ------------------------------------- | ----------------------------------------------------------------- |
| `results/<A>_to_<B>.chain.gz`         | UCSC chain with target=A; lifts A coordinates onto B              |
| `results/<B>_to_<A>.chain.gz`         | UCSC chain with target=B; lifts B coordinates onto A              |
| `results/<A>_chain_stats.tsv`         | per-target-sequence summary for the A-target chain                |
| `results/<B>_chain_stats.tsv`         | per-target-sequence summary for the B-target chain                |
| `results/align/target_<A>.paf`        | raw minimap2 PAF, target=A (kept for inspection / re-chaining)    |
| `results/align/target_<B>.paf`        | raw minimap2 PAF, target=B                                        |
| `staged/<A>.fa{,.fai}`                | decompressed + indexed input A                                    |
| `staged/<B>.fa{,.fai}`                | decompressed + indexed input B                                    |

Use as:

```bash
liftOver input.bed results/<A>_to_<B>.chain.gz output.bed unmapped.bed
```

`input.bed` is in A coordinates; `output.bed` is in B coordinates.

## Required tools

| Tool        | Tested version | Notes                                              |
| ----------- | -------------- | -------------------------------------------------- |
| `minimap2`  | >= 2.30         | https://github.com/lh3/minimap2                    |
| `transanno` | 0.4.5          | https://github.com/informationsea/transanno        |
| `chainSort` | UCSC tools     | https://hgdownload.soe.ucsc.edu/admin/exe/         |
| `samtools`  | >= 1.19         | https://github.com/samtools/samtools               |
| `snakemake` | >= 8            | https://snakemake.readthedocs.io                   |
| `python3`   | >= 3.8          | for the chain-stats script                         |

The pipeline calls each tool by name on `$PATH` by default; override any
binary path in `config.yaml` if needed (`minimap2:`, `transanno:`,
`chainsort:`, `samtools:`).

## Quick start

```bash
git clone https://github.com/mr-eyes/2026-asm2chain.git
cd 2026-asm2chain
cp config.example.yaml config.yaml
# edit config.yaml: paths to assembly_a, assembly_b, labels, preset
snakemake -s Snakefile --configfile config.yaml -n           # dry-run
snakemake -s Snakefile --configfile config.yaml --cores 32   # local
```

SLURM:

```bash
snakemake -s Snakefile --configfile config.yaml \
    --executor slurm --jobs 8 \
    --default-resources slurm_account=<account> slurm_partition=<partition>
```

## Configuration reference

See `config.example.yaml` for an annotated template.

| Key                  | Required | Default        | Description                                                                 |
| -------------------- | -------- | -------------- | --------------------------------------------------------------------------- |
| `label_a`            | yes      | -              | short identifier for assembly A (`[A-Za-z0-9._-]`)                          |
| `label_b`            | yes      | -              | short identifier for assembly B (`[A-Za-z0-9._-]`, must differ from `label_a`) |
| `assembly_a`         | yes      | -              | path to assembly A FASTA (plain or `.gz`)                                   |
| `assembly_b`         | yes      | -              | path to assembly B FASTA (plain or `.gz`)                                   |
| `minimap2_preset`    | no       | `asm5`         | one of `asm5`, `asm10`, `asm20` (see below)                                 |
| `align_threads`      | no       | `32`           | `minimap2 -t` value                                                          |
| `align_mem_mb`       | no       | `98304`        | memory request (MB) for each alignment job                                  |
| `align_runtime_min`  | no       | `720`          | runtime budget (min) for each alignment job                                 |
| `chain_mem_mb`       | no       | `16384`        | memory request for chain / stats jobs                                       |
| `chain_runtime_min`  | no       | `120`          | runtime budget for chain / stats jobs                                       |
| `minimap2`           | no       | `minimap2`     | binary path override                                                        |
| `transanno`          | no       | `transanno`    | binary path override                                                        |
| `chainsort`          | no       | `chainSort`    | binary path override                                                        |
| `samtools`           | no       | `samtools`     | binary path override                                                        |
| `results_dir`        | no       | `results`      | output root                                                                  |
| `logs_dir`           | no       | `logs`         | log root                                                                     |
| `stage_dir`          | no       | `staged`       | staged-FASTA root                                                            |

### minimap2 preset choice

Pick the preset to match the divergence of your two assemblies. As a rule
of thumb: `asm5` for same-species curated assemblies (~0.1% divergence),
`asm10` for closely related species (~1%), and `asm20` for more divergent
genomes (~5%) [REF]. The choice materially affects chain coverage and
fragmentation; document the preset alongside published chains.

[REF] https://lh3.github.io/minimap2/minimap2.html

## Output format

### Chain files

Standard UCSC chain format (gzipped) as defined at
https://genome.ucsc.edu/goldenPath/help/chain.html. Each chain header
encodes target (`tName`/`tSize`/`tStart`/`tEnd`/`tStrand`) and query
(`qName`/`qSize`/`qStart`/`qEnd`/`qStrand`); blocks follow as
`size dt dq` lines (last block has only `size`). Chains are sorted by
score with `chainSort` before gzipping.

### Chain stats TSV

One row per target sequence (i.e. column 1 from `samtools faidx` of the
target FASTA, in fai order):

| Column                    | Description                                                                  |
| ------------------------- | ---------------------------------------------------------------------------- |
| `target_seq`              | name of the target sequence                                                  |
| `target_size_bp`          | total length of the target sequence (from `samtools faidx`)                  |
| `n_chains`                | number of chain records whose target is this sequence                        |
| `n_unique_query_seqs`     | number of distinct query sequences hit                                       |
| `total_target_aligned_bp` | sum of ungapped block sizes on the target side                               |
| `target_aligned_fraction` | `total_target_aligned_bp / target_size_bp`                                   |
| `top_query_seq`           | query sequence contributing the most aligned bp to this target               |
| `top_query_aligned_bp`    | aligned bp contributed by `top_query_seq`                                    |

Target sequences with zero chains are still emitted (all-zero row) so the
TSV mirrors the target FASTA directly.

## Method

The pipeline is a thin wrapper around four standard steps:

```text
stage_assembly        zcat / symlink + samtools faidx
align_to              minimap2 -cx <preset> --cs target.fa query.fa > paf
paf_to_chain          transanno minimap2chain paf -o unsorted.chain
sort_and_gzip_chain   chainSort unsorted.chain stdout | gzip > sorted.chain.gz
chain_stats           scripts/chain_stats.py
```

Both directions (`target=A` and `target=B`) are produced. The pipeline
performs two independent alignments rather than `chainSwap`-ing a single
PAF, because minimap2 anchors differ asymmetrically between the two
directions and a freshly-aligned chain is usually preferred for downstream
liftOver work.

## Citation

If you use this tool, please cite:

> Abuelanin M, Kaya G, Lake JA, Lambert C, Wu MV, Berendzen KM, Krasheninnikova K, Wood JMD, Solomon NG, Donaldson ZR, Bales KL, Howe K, Korlach J, Manoli D, Tollkuhn J, Dennis MY. Single-library chromosome-scale diploid assemblies of vole genomes resolve a species-specific duplication implicated in pair bonding. bioRxiv 2026.03.13.711624. https://doi.org/10.64898/2026.03.13.711624

## License

MIT. See `LICENSE`.
