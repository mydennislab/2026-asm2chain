#!/usr/bin/env python3
"""
Per-target-sequence statistics for a UCSC chain file.

For each target sequence in the chain, report:
  target_seq               name of the target sequence (column 6 in the
                           chain header line "chain ... tName ...")
  target_size_bp           total length of that target sequence (from the
                           --target-fai)
  n_chains                 number of chain records whose target is this seq
  n_unique_query_seqs      number of distinct query sequences that have at
                           least one chain to this target
  total_target_aligned_bp  sum of ungapped block sizes (column 1 of each
                           chain block line) for chains on this target
  target_aligned_fraction  total_target_aligned_bp / target_size_bp
  top_query_seq            the query sequence contributing the most
                           aligned bp to this target
  top_query_aligned_bp     that contribution

Accepts plain or gzipped chain. Writes TSV to --output.
"""

import argparse
import gzip
import sys
from collections import defaultdict


def parse_fai(path):
    sizes = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            sizes[parts[0]] = int(parts[1])
    return sizes


def iter_chains(path):
    """Yield (header_dict, total_aligned_bp_on_target) per chain.

    Chain format (UCSC):
        chain score tName tSize tStrand tStart tEnd qName qSize qStrand
            qStart qEnd id
        size dt dq    (one block per line, last block has just `size`)
        <blank line>
    Some chain producers omit the blank separator and emit headers back-to-
    back; handle both. The `size` field (block length) is identical on
    target and query sides by definition, so we accumulate it as the
    aligned-bp contribution.
    """
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        header = None
        aligned = 0
        for raw in f:
            line = raw.rstrip("\n")
            if not line:
                if header is not None:
                    yield header, aligned
                    header = None
                    aligned = 0
                continue
            if line.startswith("chain"):
                # Defensive: chains may not be blank-separated.
                if header is not None:
                    yield header, aligned
                fields = line.split()
                # Indices:
                #  0       1     2      3      4        5       6     7      8      9        10      11    12
                # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
                header = {
                    "score":   int(fields[1]),
                    "tName":   fields[2],
                    "tSize":   int(fields[3]),
                    "tStrand": fields[4],
                    "tStart":  int(fields[5]),
                    "tEnd":    int(fields[6]),
                    "qName":   fields[7],
                    "qSize":   int(fields[8]),
                    "qStrand": fields[9],
                    "qStart":  int(fields[10]),
                    "qEnd":    int(fields[11]),
                    "id":      fields[12] if len(fields) > 12 else "",
                }
                aligned = 0
            else:
                parts = line.split()
                aligned += int(parts[0])
        if header is not None:
            yield header, aligned


def main():
    ap = argparse.ArgumentParser(
        description="Per-target-seq chain statistics for a UCSC chain file."
    )
    ap.add_argument("--chain", required=True, help="UCSC chain (plain or .gz).")
    ap.add_argument(
        "--target-fai",
        required=True,
        help="samtools faidx of the target FASTA (provides target_size_bp).",
    )
    ap.add_argument(
        "--query-fai",
        required=True,
        help=(
            "samtools faidx of the query FASTA. Currently used only as a "
            "sanity check that query names in the chain are recognized; the "
            "stats themselves are target-anchored."
        ),
    )
    ap.add_argument("--output", required=True, help="Output TSV path.")
    args = ap.parse_args()

    target_sizes = parse_fai(args.target_fai)
    query_sizes = parse_fai(args.query_fai)

    n_chains = defaultdict(int)
    queries_hit = defaultdict(set)
    aligned_total = defaultdict(int)
    aligned_by_query = defaultdict(lambda: defaultdict(int))

    unknown_q = 0
    n_total = 0
    for header, aligned in iter_chains(args.chain):
        t = header["tName"]
        q = header["qName"]
        n_chains[t] += 1
        queries_hit[t].add(q)
        aligned_total[t] += aligned
        aligned_by_query[t][q] += aligned
        if q not in query_sizes:
            unknown_q += 1
        n_total += 1

    print(
        f"[chain_stats] parsed {n_total} chains across {len(n_chains)} "
        f"target seqs ({unknown_q} chains referenced query seqs not in "
        f"--query-fai)",
        file=sys.stderr,
    )

    cols = [
        "target_seq",
        "target_size_bp",
        "n_chains",
        "n_unique_query_seqs",
        "total_target_aligned_bp",
        "target_aligned_fraction",
        "top_query_seq",
        "top_query_aligned_bp",
    ]
    with open(args.output, "w") as out:
        out.write("\t".join(cols) + "\n")
        # Emit rows in fai (target) order so users see chromosomes first.
        ordered = list(target_sizes.keys())
        # Append any target names that appeared in chains but not the fai
        # (shouldn't happen with consistent inputs).
        for t in n_chains:
            if t not in target_sizes:
                ordered.append(t)
        for t in ordered:
            size = target_sizes.get(t, 0)
            n = n_chains.get(t, 0)
            if n == 0:
                out.write(
                    "\t".join(
                        [t, str(size), "0", "0", "0", "0.0000", "", "0"]
                    )
                    + "\n"
                )
                continue
            tot = aligned_total[t]
            frac = (tot / size) if size else 0.0
            by_q = aligned_by_query[t]
            top_q, top_bp = max(by_q.items(), key=lambda kv: kv[1])
            out.write(
                "\t".join(
                    [
                        t,
                        str(size),
                        str(n),
                        str(len(queries_hit[t])),
                        str(tot),
                        f"{frac:.4f}",
                        top_q,
                        str(top_bp),
                    ]
                )
                + "\n"
            )


if __name__ == "__main__":
    main()
