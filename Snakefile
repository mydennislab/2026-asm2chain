"""
asm2chain
=========
Bidirectional UCSC liftOver chains from two assemblies.

Pipeline:
  1. Stage both assemblies (decompress if gzipped) and index with samtools faidx.
  2. Run minimap2 -cx <preset> --cs in both directions (A as target, B as target).
  3. Convert each PAF to a UCSC chain with `transanno minimap2chain`.
  4. Sort with UCSC `chainSort` and gzip.
  5. Emit per-target-sequence chain stats (count, unique query partners,
     aligned bp on target, top partner) for each direction.

Outputs (under results/):
  {labelA}_to_{labelB}.chain.gz   target=A; lifts A coords -> B coords
  {labelB}_to_{labelA}.chain.gz   target=B; lifts B coords -> A coords
  {labelA}_chain_stats.tsv        per-target-seq stats for chain with A as target
  {labelB}_chain_stats.tsv        per-target-seq stats for chain with B as target

Run:
  snakemake -s Snakefile --configfile config.yaml -n        # dry-run
  snakemake -s Snakefile --configfile config.yaml --cores 32
"""

import os
import sys

# Configuration
# No `configfile:` directive: pass --configfile explicitly so a stray
# config.yaml in the working dir is never picked up silently.

# Required keys (validated below)
REQUIRED_KEYS = ("assembly_a", "assembly_b", "label_a", "label_b")
for k in REQUIRED_KEYS:
    if k not in config:
        sys.exit(
            f"[asm2chain] config error: missing required key '{k}'. "
            "See config.example.yaml."
        )

LABEL_A = str(config["label_a"]).strip()
LABEL_B = str(config["label_b"]).strip()
if not LABEL_A or not LABEL_B:
    sys.exit("[asm2chain] config error: label_a and label_b must be non-empty.")
if LABEL_A == LABEL_B:
    sys.exit("[asm2chain] config error: label_a and label_b must differ.")
for lab in (LABEL_A, LABEL_B):
    bad = [c for c in lab if not (c.isalnum() or c in "._-")]
    if bad:
        sys.exit(
            f"[asm2chain] config error: label '{lab}' contains invalid characters "
            f"({''.join(sorted(set(bad)))}). Use [A-Za-z0-9._-] only."
        )

ASM_A = os.path.expanduser(str(config["assembly_a"]))
ASM_B = os.path.expanduser(str(config["assembly_b"]))

PRESET = str(config.get("minimap2_preset", "asm5"))
ALIGN_THREADS = int(config.get("align_threads", 32))
ALIGN_MEM_MB = int(config.get("align_mem_mb", 96 * 1024))
ALIGN_RUNTIME_MIN = int(config.get("align_runtime_min", 12 * 60))

CHAIN_MEM_MB = int(config.get("chain_mem_mb", 16 * 1024))
CHAIN_RUNTIME_MIN = int(config.get("chain_runtime_min", 2 * 60))

# Tool binaries (allow override; fall back to PATH lookup).
MINIMAP2 = str(config.get("minimap2", "minimap2"))
TRANSANNO = str(config.get("transanno", "transanno"))
CHAINSORT = str(config.get("chainsort", "chainSort"))
SAMTOOLS = str(config.get("samtools", "samtools"))

RESULTS_DIR = str(config.get("results_dir", "results"))
LOGS_DIR = str(config.get("logs_dir", "logs"))
STAGE_DIR = str(config.get("stage_dir", "staged"))


def asm_path(label):
    return ASM_A if label == LABEL_A else ASM_B


# Targets
rule all:
    input:
        f"{RESULTS_DIR}/{LABEL_A}_to_{LABEL_B}.chain.gz",
        f"{RESULTS_DIR}/{LABEL_B}_to_{LABEL_A}.chain.gz",
        f"{RESULTS_DIR}/{LABEL_A}_chain_stats.tsv",
        f"{RESULTS_DIR}/{LABEL_B}_chain_stats.tsv",


# 1. Stage assemblies (decompress + index)
rule stage_assembly:
    """
    Stage an input assembly into a plain (uncompressed) FASTA and index it
    with samtools faidx. Handles .gz transparently; symlinks otherwise to
    keep IO cheap.
    """
    input:
        src=lambda wc: asm_path(wc.label),
    output:
        fa=f"{STAGE_DIR}/{{label}}.fa",
        fai=f"{STAGE_DIR}/{{label}}.fa.fai",
    log:
        f"{LOGS_DIR}/stage_{{label}}.log",
    threads: 1
    resources:
        runtime=60,
        mem_mb=4 * 1024,
    shell:
        r"""
        mkdir -p "$(dirname {output.fa})" "$(dirname {log})"
        : > {log}
        if [ ! -s "{input.src}" ]; then
            echo "[asm2chain] input missing or empty: {input.src}" >&2
            exit 2
        fi
        case "{input.src}" in
            *.gz)
                echo "[asm2chain] decompressing {input.src} -> {output.fa}" >> {log}
                zcat "{input.src}" > "{output.fa}" 2>> {log}
                ;;
            *)
                echo "[asm2chain] symlinking {input.src} -> {output.fa}" >> {log}
                ln -sf "$(readlink -f {input.src})" "{output.fa}" 2>> {log}
                ;;
        esac
        {SAMTOOLS} faidx "{output.fa}" 2>> {log}
        """


# 2. Alignment: minimap2 in both directions
# We name the PAF by its target. align_to_{label} produces a PAF whose target
# (column 6) is the {label} assembly and whose query is the other one.
rule align_to:
    """
    minimap2 -cx <preset> --cs target.fa query.fa > paf
    The target here is {label}; the query is the other assembly.
    """
    input:
        target=f"{STAGE_DIR}/{{label}}.fa",
        target_fai=f"{STAGE_DIR}/{{label}}.fa.fai",
        query=lambda wc: f"{STAGE_DIR}/{LABEL_B if wc.label == LABEL_A else LABEL_A}.fa",
    output:
        paf=f"{RESULTS_DIR}/align/target_{{label}}.paf",
    log:
        f"{LOGS_DIR}/align_target_{{label}}.log",
    threads: ALIGN_THREADS
    params:
        preset=PRESET,
    resources:
        runtime=ALIGN_RUNTIME_MIN,
        mem_mb=ALIGN_MEM_MB,
    shell:
        r"""
        mkdir -p "$(dirname {output.paf})" "$(dirname {log})"
        {MINIMAP2} -t {threads} -cx {params.preset} --cs \
            {input.target} {input.query} > {output.paf} 2> {log}
        """


# 3. PAF -> UCSC chain (transanno)
rule paf_to_chain:
    """
    Convert a target-anchored PAF to an unsorted UCSC chain via
    `transanno minimap2chain`. The chain inherits the target/query roles
    from the PAF, so target_{label}.paf -> chain whose target is {label}.
    """
    input:
        paf=f"{RESULTS_DIR}/align/target_{{label}}.paf",
    output:
        chain=temp(f"{RESULTS_DIR}/chain_unsorted/target_{{label}}.unsorted.chain"),
    log:
        f"{LOGS_DIR}/paf2chain_target_{{label}}.log",
    threads: 1
    resources:
        runtime=CHAIN_RUNTIME_MIN,
        mem_mb=CHAIN_MEM_MB,
    shell:
        r"""
        mkdir -p "$(dirname {output.chain})" "$(dirname {log})"
        {TRANSANNO} minimap2chain {input.paf} --output {output.chain} 2> {log}
        """


# 4. Sort + gzip the chain
# Naming: a chain whose target is X lifts X coordinates onto the query, which
# is the *other* assembly. So target_A.chain lifts A -> B and is published as
# {labelA}_to_{labelB}.chain.gz.
def _sorted_chain_target(wc):
    # source label is the target of the chain
    return wc.src


rule sort_and_gzip_chain:
    """
    Sort the unsorted chain by score with UCSC chainSort and gzip the result.
    Output name reflects coordinate-lift direction:
      <target_label>_to_<query_label>.chain.gz lifts target -> query coords.
    """
    input:
        chain=f"{RESULTS_DIR}/chain_unsorted/target_{{src}}.unsorted.chain",
    output:
        chain_gz=f"{RESULTS_DIR}/{{src}}_to_{{dst}}.chain.gz",
    log:
        f"{LOGS_DIR}/sort_chain_{{src}}_to_{{dst}}.log",
    threads: 1
    resources:
        runtime=CHAIN_RUNTIME_MIN,
        mem_mb=CHAIN_MEM_MB,
    wildcard_constraints:
        # Restrict to the configured (src, dst) pairs only. Two valid combos.
        src=f"({LABEL_A}|{LABEL_B})",
        dst=f"({LABEL_A}|{LABEL_B})",
    shell:
        r"""
        mkdir -p "$(dirname {output.chain_gz})" "$(dirname {log})"
        {CHAINSORT} {input.chain} stdout 2> {log} | gzip -c > {output.chain_gz}
        """


# 5. Per-target chain stats
rule chain_stats:
    """
    Emit one row per target sequence with chain count, unique query partners
    hit, total aligned bp on target side, and the dominant query partner.
    See scripts/chain_stats.py for column definitions.
    """
    input:
        chain=lambda wc: (
            f"{RESULTS_DIR}/{wc.label}_to_"
            f"{LABEL_B if wc.label == LABEL_A else LABEL_A}.chain.gz"
        ),
        target_fai=f"{STAGE_DIR}/{{label}}.fa.fai",
        query_fai=lambda wc: (
            f"{STAGE_DIR}/{LABEL_B if wc.label == LABEL_A else LABEL_A}.fa.fai"
        ),
    output:
        tsv=f"{RESULTS_DIR}/{{label}}_chain_stats.tsv",
    log:
        f"{LOGS_DIR}/chain_stats_{{label}}.log",
    threads: 1
    resources:
        runtime=30,
        mem_mb=8 * 1024,
    shell:
        r"""
        mkdir -p "$(dirname {output.tsv})" "$(dirname {log})"
        python3 scripts/chain_stats.py \
            --chain {input.chain} \
            --target-fai {input.target_fai} \
            --query-fai {input.query_fai} \
            --output {output.tsv} 2> {log}
        """
