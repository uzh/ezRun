#!/usr/bin/env python3
"""Register a SUSHI dataset's files in B-Fabric as one workunit + resources + dataset.

Self-contained replacement for btools' register_sushi_dataset_into_bfabric, using
bfabricPy only. Meant for apps that deliver raw data into short-term storage
(application 372, e.g. EnaApp, MergeRunData), where the app runs BEFORE the
files are copied to gStore. Sizes and md5 checksums are therefore taken from the
local copies in --local-dir (the job's working directory). A gStore copy is used
as a fallback, so the script also works for back-filling delivered data.

What gets created (same layout btools produces):
  workunit  name=description=<name>, application, container, status
            processing -> available
  resource  one per unique "[File]" path: name=basename, relativepath, size,
            filechecksum, storage, status available; expirationdate +5y if > 1 GiB
  dataset   name=<name>, one attribute per column ("[File]" columns typed Resource),
            one item per row, linked to the workunit
If anything fails after the workunit exists, the workunit is set to failed and
deleted (that cascades to its resources and dataset) and the script exits non-zero.

stdout: one JSON line {"workunit_id", "dataset_id", "n_resources", "n_rows"}.
"""

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

GSTORE_ROOT = "/srv/gstore/projects"
CHUNK = 100  # B-Fabric accepts at most 100 objects per save call
EXCLUDED_COLUMNS = {"ID", "Project"}
GIB = 1024**3


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", required=True, help="SUSHI dataset.tsv (gStore-relative [File] paths)")
    p.add_argument("--container-id", required=True, type=int, help="B-Fabric project or order id")
    p.add_argument("--name", required=True, help="workunit and dataset name")
    p.add_argument("--application", type=int, default=372, help="B-Fabric application id (default 372)")
    p.add_argument("--storage-id", type=int, default=5, help="B-Fabric storage id (default 5 = gStore)")
    p.add_argument("--local-dir", default=".", help="directory holding the files before delivery")
    p.add_argument("--input-dataset-id", type=int, default=0, help="B-Fabric input dataset id (optional)")
    p.add_argument("--config-file", default="~/.bfabricpy.yml")
    p.add_argument("--config-env", default=None, help="config environment, e.g. PRODUCTION or TEST")
    p.add_argument("--threads", type=int, default=8, help="parallel md5 workers")
    p.add_argument("--no-md5", action="store_true", help="do not compute checksums")
    p.add_argument("--dry-run", action="store_true", help="resolve files and print the plan; write nothing")
    return p.parse_args()


def read_dataset(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    rows = [[c.rstrip("\r\n").strip() for c in r] for r in rows if any(c.strip() for c in r)]
    if len(rows) < 2:
        raise ValueError(f"{path}: no data rows")
    header, body = rows[0], rows[1:]
    for i, r in enumerate(body, start=2):
        if len(r) != len(header):
            raise ValueError(f"{path}: line {i} has {len(r)} fields, header has {len(header)}")
    return header, body


def is_file_column(col):
    return "[File]" in col


def clean_name(col):
    return re.sub(r"\s*\[[^\]]*\]\s*$", "", col).strip()


def resolve_local(relpath, local_dir):
    """Where the bytes are right now: job working dir first, then gStore."""
    for cand in (
        os.path.join(local_dir, os.path.basename(relpath)),
        os.path.join(local_dir, relpath),
        os.path.join(GSTORE_ROOT, relpath),
    ):
        if os.path.exists(cand):
            return cand
    return None


def md5_of(path):
    if os.path.isdir(path):
        return "folder"
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def size_of(path):
    if os.path.isdir(path):
        total = 0
        for root, _, files in os.walk(path):
            total += sum(os.path.getsize(os.path.join(root, f)) for f in files)
        return total
    return os.path.getsize(path)


def collect_resources(header, body, local_dir):
    file_cols = [i for i, c in enumerate(header) if is_file_column(c)]
    if not file_cols:
        raise ValueError("dataset has no '[File]' column")
    seen, resources, missing = set(), [], []
    for row in body:
        for i in file_cols:
            for relpath in (p.strip() for p in row[i].split(",")):
                if not relpath or relpath in seen:
                    continue  # duplicate paths would violate B-Fabric's (storage, path, workunit) key
                seen.add(relpath)
                local = resolve_local(relpath, local_dir)
                if local is None:
                    missing.append(relpath)
                resources.append({"relativepath": relpath, "local": local})
    if missing:
        raise FileNotFoundError(f"{len(missing)} file(s) not found locally or on gStore, e.g. {missing[:3]}")
    return resources


def dataset_payload(header, body):
    keep = [i for i, c in enumerate(header) if c not in EXCLUDED_COLUMNS]
    attributes = [
        {"name": clean_name(header[i]), "position": pos, "type": "Resource" if is_file_column(header[i]) else "String"}
        for pos, i in enumerate(keep, start=1)
    ]
    items = [
        {"position": r_pos, "field": [
            {"attributeposition": pos, "value": row[i] if row[i] else "-"} for pos, i in enumerate(keep, start=1)
        ]}
        for r_pos, row in enumerate(body, start=1)
    ]
    return attributes, items


def chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def main():
    a = parse_args()
    header, body = read_dataset(a.dataset)
    resources = collect_resources(header, body, a.local_dir)

    for r in resources:
        r["size"] = size_of(r["local"])
    if a.no_md5:
        for r in resources:
            r["md5"] = ""
    else:
        with ProcessPoolExecutor(max_workers=max(1, a.threads)) as pool:
            for r, md5 in zip(resources, pool.map(md5_of, [r["local"] for r in resources])):
                r["md5"] = md5
    attributes, items = dataset_payload(header, body)

    if a.dry_run:
        print(json.dumps({
            "dry_run": True, "n_rows": len(body), "n_resources": len(resources),
            "total_bytes": sum(r["size"] for r in resources),
            "attributes": attributes, "first_resources": resources[:3],
        }, indent=1))
        return 0

    from bfabric import Bfabric  # imported late so --dry-run works without credentials

    env = a.config_env or os.environ.get("BFABRICPY_CONFIG_ENV") or "default"
    bf = Bfabric.connect(config_file_path=os.path.expanduser(a.config_file), config_file_env=env)

    wu = {"name": a.name, "description": a.name, "applicationid": a.application,
          "containerid": a.container_id, "status": "processing"}
    if a.input_dataset_id > 0:
        wu["inputdatasetid"] = a.input_dataset_id
    workunit_id = bf.save("workunit", wu).to_list_dict()[0]["id"]

    try:
        expiry = (dt.date.today().replace(year=dt.date.today().year + 5)).isoformat()
        objs = []
        for r in resources:
            o = {"name": os.path.basename(r["relativepath"].rstrip("/")), "relativepath": r["relativepath"],
                 "size": r["size"], "filechecksum": r["md5"], "storageid": a.storage_id,
                 "status": "available", "workunitid": workunit_id}
            if r["size"] > GIB:
                o["expirationdate"] = expiry
            objs.append(o)
        n_created = 0
        for part in chunks(objs, CHUNK):
            n_created += len(bf.save("resource", part).to_list_dict())
        if n_created != len(objs):
            raise RuntimeError(f"created {n_created} of {len(objs)} resources")

        ds = bf.save("dataset", {"name": a.name, "containerid": a.container_id, "workunitid": workunit_id,
                                 "attribute": attributes, "item": items}).to_list_dict()[0]
        bf.save("workunit", {"id": workunit_id, "status": "available"})
    except BaseException:
        try:
            # B-Fabric refuses to delete a workunit while it is 'processing'
            bf.save("workunit", {"id": workunit_id, "status": "failed"})
            bf.delete("workunit", workunit_id)
            print(f"registration failed; deleted workunit {workunit_id}", file=sys.stderr)
        except Exception as e:  # leave a trace so it can be cleaned up by hand
            print(f"registration failed; could NOT delete workunit {workunit_id}: {e}", file=sys.stderr)
        raise

    print(json.dumps({"workunit_id": int(workunit_id), "dataset_id": int(ds["id"]),
                      "n_resources": n_created, "n_rows": len(body)}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
