#!/usr/bin/env python3
"""Convert a Kraken-style report (e.g. Bracken's *.bracken.report.txt) into a
MetaPhlAn-style lineage table for HUMAnN translated-only mode.

The 6 columns of a Kraken report are: pct, reads_clade, reads_at, rank_code,
taxid, indented_name. We walk it depth-first, maintain a per-rank lineage
stack, and emit a MetaPhlAn line per taxon at ranks D|P|C|O|F|G|S (mapped to
MetaPhlAn's k__|p__|c__|o__|f__|g__|s__). Sub-ranks (S1, S2, K, R, R1, ...)
are skipped — they have no MetaPhlAn equivalent.

Output goes to stdout (or -o). Default header tag is the one HUMAnN 4 alpha
expects; override with --header if you need a different vintage.

NOTE: HUMAnN translated-only ignores the profile *content* — only the header
line matters. This converter exists for audit trail / SUSHI app aesthetics.
"""

import sys
import re
import argparse

RANK_MAP   = {"D": "k__", "P": "p__", "C": "c__", "O": "o__",
              "F": "f__", "G": "g__", "S": "s__"}
RANK_ORDER = ["D", "P", "C", "O", "F", "G", "S"]
# Some Kraken DBs (GTDB-derived) already prefix names with d__/p__/.../s__.
# Strip any such prefix before we attach our own MetaPhlAn one.
GTDB_PREFIX_RE = re.compile(r"^[dpcofgs]__")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("report", help="Bracken/Kraken report file")
    ap.add_argument("-o", "--out", default="-", help="output path (default stdout)")
    ap.add_argument("--header", default="#mpa_vOct22_CHOCOPhlAnSGB_202403",
                    help="MetaPhlAn DB header tag (default: vOct22)")
    args = ap.parse_args()

    lineage = {}                                       # rank_code -> (name, taxid)
    out_lines = [args.header]

    with open(args.report) as fh:
        for raw in fh:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            pct, _rc, _ra, rank, taxid, raw_name = parts[:6]
            if rank not in RANK_MAP:
                continue
            name = GTDB_PREFIX_RE.sub("", raw_name.strip()).replace(" ", "_")
            # Clear all deeper ranks before recording this one
            idx = RANK_ORDER.index(rank)
            for r in RANK_ORDER[idx + 1:]:
                lineage.pop(r, None)
            lineage[rank] = (name, taxid)
            try:
                pct_val = float(pct)
            except ValueError:
                continue
            if pct_val <= 0:
                continue
            tags, taxids = [], []
            for r in RANK_ORDER:
                if r in lineage:
                    nm, ti = lineage[r]
                    tags.append(f"{RANK_MAP[r]}{nm}")
                    taxids.append(ti)
            out_lines.append("\t".join(["|".join(tags), "|".join(taxids),
                                        f"{pct_val:.5f}"]))

    blob = "\n".join(out_lines) + "\n"
    if args.out == "-":
        sys.stdout.write(blob)
    else:
        with open(args.out, "w") as f:
            f.write(blob)


if __name__ == "__main__":
    main()
