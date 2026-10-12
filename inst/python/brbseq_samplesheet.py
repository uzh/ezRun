#!/usr/bin/env python3
"""Build the BRB-seq sample sheet for an order from its B-Fabric library plate.

The sheet is the input contract of the BRBseq SUSHI/ezRun app: one row per
plate well with the 14 nt BRB well barcode, the well, the library name and its
role. The app finds it by convention at

    /srv/GT/gstore/projects/p<project>/o<order>_metaData/o<order>_BRB_SampleSheet.tsv

What B-Fabric holds (checked on order 42799, plate 7248):
  * Each well is a "Library - Illumina" sample on a library plate; its
    `multiplexid` is the 14 nt BRB barcode (plain text), its `tubeid` points to
    the biological sample (e.g. 42799/1#5), the description holds the input
    ("input 150ng in 10ul").
  * FGCZ control and negative wells live in the PROJECT container, with tubeids
    FGCZ_control_RNA_v3_<n> / FGCZ_control_negative_<n>.
  * Library names (Plate_<id>_<n>) are NOT unique across containers, so rows
    are named after the biological sample (controls: after the tube).
  * Which wells went into which sequencing pool is NOT recorded -> `Pool` is
    left empty and the app assigns pools from the data.

The plate version (V5A..V5D) is detected by matching each (well, barcode) pair
against the kit file; V5A and V5D share sequences but not well positions.

Usage (bfabricPy, e.g. /usr/local/ngseq/bin/python):
    brbseq_samplesheet.py --order 42799 [--plate 7248] [--out sheet.tsv] [--deliver]
--deliver copies the sheet to the gstore convention path with `g-req copynow -f`.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from collections import Counter

KIT_FILE = "/srv/GT/databases/BRB/whiteList_SampleBarcodes/BRB_V5_barcodes.tsv"
GSTORE = "/srv/gstore/projects"
BARCODE_LEN = 14
COLUMNS = ["Name", "Barcode", "Well", "Role", "Pool", "Input [Factor]",
           "Condition [Factor]", "Sample Id [B-Fabric]", "Tube Id", "Library Id",
           "Plate Id", "Plate Version"]


def norm_well(pos):
    """'A1' / 'A01' -> 'A01'."""
    m = re.fullmatch(r"([A-Pa-p])0*(\d{1,2})", pos.strip())
    if not m:
        raise ValueError(f"unexpected plate position: {pos!r}")
    return f"{m.group(1).upper()}{int(m.group(2)):02d}"


def well_key(well):
    """Column-major plate order: A01, B01, ..., H01, A02, ..."""
    return (int(well[1:]), well[0])


def safe_name(x):
    return re.sub(r"[^A-Za-z0-9_-]+", "_", x).strip("_")


def read_kit(path):
    kit = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            kit.setdefault(row["Plate"], {})[row["Well"]] = row["Barcode"]
    return kit


def detect_plate_version(wells, kit):
    """Return (version, n_matching) for the kit plate whose well->barcode map fits best."""
    scores = {v: sum(kit[v].get(w) == bc for w, bc in wells) for v in kit}
    best = max(scores, key=scores.get)
    return best, scores[best], scores


def read_many(bf, endpoint, ids, chunk=100):
    out = []
    ids = list(ids)
    for i in range(0, len(ids), chunk):
        out += bf.read(endpoint, {"id": ids[i:i + chunk]}, max_results=None).to_list_dict()
    return out


def find_plates(bf, order_id):
    """Library plates holding the order's 14 nt BRB well libraries."""
    libs = bf.read("sample", {"containerid": order_id, "type": "Library - Illumina"},
                   max_results=None).to_list_dict()
    brb = [l for l in libs if len(l.get("multiplexid") or "") == BARCODE_LEN]
    if not brb:
        sys.exit(f"order {order_id}: no 'Library - Illumina' samples with a {BARCODE_LEN} nt barcode")
    plates, covered = {}, set()
    for l in brb:                      # one plate lookup per plate, not per well
        if int(l["id"]) in covered:
            continue
        for p in bf.read("plate", {"sampleid": int(l["id"])}).to_list_dict():
            plates[int(p["id"])] = p["name"]
            covered |= {int(s["id"]) for s in p.get("sample", [])}
    return plates


def build_sheet(bf, order_id, plate_ids, kit):
    rows = []
    for plate_id in plate_ids:
        plate = bf.read("plate", {"id": plate_id}).to_list_dict()[0]
        pos = {int(s["id"]): norm_well(s["_gridposition"]) for s in plate["sample"]}
        libs = read_many(bf, "sample", pos)
        tubes = sorted({l.get("tubeid") for l in libs if "#" in (l.get("tubeid") or "")})
        bio = {}
        for i in range(0, len(tubes), 100):
            for b in bf.read("sample", {"tubeid": tubes[i:i + 100]}, max_results=None).to_list_dict():
                if b.get("type", "").startswith("Biological"):
                    bio[b["tubeid"]] = b
        wells = []
        for l in libs:
            bc = (l.get("multiplexid") or "").strip().upper()
            if len(bc) != BARCODE_LEN:
                continue                       # not a BRB well library (e.g. a pool)
            tube = l.get("tubeid") or ""
            b = bio.get(tube)
            if tube.startswith("FGCZ_control_negative"):
                role = "neg"
            elif tube.startswith("FGCZ_control"):
                role = "control"
            else:
                role = "sample"
            m = re.search(r"input\s*([\d.]+)\s*ng", l.get("description") or "", re.I)
            wells.append({
                "Name": safe_name(b["name"] if b else tube or l["name"]),
                "Barcode": bc,
                "Well": pos[int(l["id"])],
                "Role": role,
                "Pool": "",
                "Input [Factor]": f"{m.group(1)}ng" if m else "",
                "Condition [Factor]": "",
                "Sample Id [B-Fabric]": b["id"] if b else "",
                "Tube Id": tube,
                "Library Id": l["id"],
                "Plate Id": plate_id,
            })
        version, n_ok, scores = detect_plate_version([(w["Well"], w["Barcode"]) for w in wells], kit)
        if n_ok < len(wells):
            bad = [w for w in wells if kit[version].get(w["Well"]) != w["Barcode"]]
            print(f"WARNING plate {plate_id}: {len(bad)} wells do not match kit {version} "
                  f"(scores {scores}): " + ", ".join(f"{w['Well']}={w['Barcode']}" for w in bad[:8]),
                  file=sys.stderr)
        for w in wells:
            w["Plate Version"] = version
        rows += wells
        print(f"plate {plate_id} '{plate['name']}': {len(wells)} BRB wells, kit {version} "
              f"({n_ok}/{len(wells)} well/barcode pairs match)", file=sys.stderr)

    # unique names: suffix the well only where a name repeats
    dup = {n for n, k in Counter(r["Name"] for r in rows).items() if k > 1}
    for r in rows:
        if r["Name"] in dup:
            r["Name"] = f"{r['Name']}_{r['Well']}"
    names = Counter(r["Name"] for r in rows)
    bcs = Counter(r["Barcode"] for r in rows)
    if max(names.values()) > 1:
        sys.exit("duplicate names after well suffix: " + ", ".join(n for n, k in names.items() if k > 1))
    if max(bcs.values()) > 1:
        print("WARNING: a barcode is used in more than one well (several plates of the same "
              "version?) -> set the Pool column before running the app", file=sys.stderr)
    return sorted(rows, key=lambda r: (str(r["Plate Id"]), well_key(r["Well"])))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--order", type=int, required=True)
    ap.add_argument("--plate", type=int, action="append", help="library plate id (repeatable); default: found from the order")
    ap.add_argument("--out", help="output TSV (default: ./o<order>_BRB_SampleSheet.tsv)")
    ap.add_argument("--kit", default=KIT_FILE)
    ap.add_argument("--deliver", action="store_true", help="copy to gstore o<order>_metaData/ with g-req")
    ap.add_argument("--config-env", default="PRODUCTION")
    a = ap.parse_args()

    from bfabric import Bfabric
    bf = Bfabric.from_config(config_env=a.config_env)
    kit = read_kit(a.kit)

    plate_ids = a.plate or sorted(find_plates(bf, a.order))
    if not plate_ids:
        sys.exit(f"order {a.order}: no library plate found")
    rows = build_sheet(bf, a.order, plate_ids, kit)

    out = a.out or f"o{a.order}_BRB_SampleSheet.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    roles = Counter(r["Role"] for r in rows)
    print(f"wrote {out}: {len(rows)} wells ({', '.join(f'{k} {v}' for k, v in sorted(roles.items()))})",
          file=sys.stderr)

    if a.deliver:
        order = bf.read("order", {"id": a.order}).to_list_dict()[0]
        project = int(order["project"]["id"])
        dest = os.path.join(GSTORE, f"p{project}", f"o{a.order}_metaData")
        target = os.path.join(dest, f"o{a.order}_BRB_SampleSheet.tsv")
        if os.path.basename(out) != os.path.basename(target):
            sys.exit(f"--deliver needs the file to be named {os.path.basename(target)}")
        print(f"g-req copynow -f {out} {dest}", file=sys.stderr)
        subprocess.run(["g-req", "-w", "copynow", "-f", os.path.abspath(out), dest], check=True)


if __name__ == "__main__":
    main()
