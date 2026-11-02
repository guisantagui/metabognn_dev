################################################################################
####################----------------------------------------####################
##########----------            MetaboReprogram             ----------##########
# Gene ID mapping and memory-light reading of genes x samples matrices.        #
# Used by expr_symb_to_ensembl.py and by the MetChange preparation scripts.    #
#                                                                              #
################################################################################

import csv
from collections import defaultdict

import numpy as np
import pandas as pd


def gem_symbol_overrides(gem):
    """Symbol -> GEM gene id from the model's own hgnc.symbol annotation. HGNC and
    Human-GEM use different Ensembl ids for a few genes (SOD2, SLC37A4, LRTOMT,
    PPIAL4C); without this they would not match the GEM."""
    return {g.annotation["hgnc.symbol"]: g.id for g in gem.genes if g.annotation.get("hgnc.symbol")}


def symbol_to_ensembl_map(hgnc, overrides=None):
    """Approved symbol -> Ensembl id; aliases and previous symbols only when
    they point to a single gene and are not an approved symbol themselves.
    overrides (e.g. gem_symbol_overrides) take precedence."""
    hgnc = hgnc[hgnc["ensembl_gene_id"].notna() & (hgnc["ensembl_gene_id"] != "")]
    mapping = dict(zip(hgnc["symbol"], hgnc["ensembl_gene_id"]))
    for col in ("prev_symbol", "alias_symbol"):
        hits = defaultdict(set)
        for ens, val in zip(hgnc["ensembl_gene_id"], hgnc[col].fillna("")):
            for s in str(val).split("|"):
                if s:
                    hits[s].add(ens)
        for s, ens in hits.items():
            if s not in mapping and len(ens) == 1:
                mapping[s] = next(iter(ens))
    mapping.update(overrides or {})
    return mapping


def read_header(path):
    with open(path) as fh:
        return next(csv.reader(fh))


def merge_rows(rows, is_log2, how):
    """Merge rows of one gene: in linear space if the data are log2, result back in log2."""
    x = np.vstack(rows).astype(np.float64)
    if is_log2:
        x = np.exp2(x)
    x = x.sum(0) if how == "sum" else x.mean(0)
    return np.log2(x) if is_log2 else x


def read_matrix(path, id_col=None, mapping=None, keep=None, is_log2=True, agg="sum", chunk_rows=100):
    """Read a genes x samples CSV in row chunks into one float32 array.

    id_col   name of the gene id column (default: the first column)
    mapping  dict old id -> new id; rows without a mapping are dropped
    keep     set of (new) ids to retain, e.g. the GEM genes; others are skipped while reading
    Rows that end up with the same id are merged (merge_rows). Memory is the
    output array (n_genes x n_samples x 4 bytes) plus one chunk.
    Returns (ids, samples, array, report) with report = one row per input row.
    """
    header = read_header(path)
    idx = 0 if id_col is None else header.index(id_col)
    samples = [h for i, h in enumerate(header) if i != idx]
    ids_in = pd.read_csv(path, usecols=[idx], dtype=str, keep_default_na=False).iloc[:, 0].tolist()
    new = [mapping.get(g, "") for g in ids_in] if mapping is not None else list(ids_in)
    status = ["kept" if n else "unmapped" for n in new]
    if keep is not None:
        for k, n in enumerate(new):
            if n and n not in keep:
                new[k], status[k] = "", "not_kept"
    out_row, count = {}, defaultdict(int)
    for n in new:
        if n:
            out_row.setdefault(n, len(out_row))
            count[n] += 1
    for k, n in enumerate(new):
        if n and count[n] > 1:
            status[k] = "merged"
    X = np.empty((len(out_row), len(samples)), np.float32)
    pending = defaultdict(list)
    row = 0
    for chunk in pd.read_csv(path, chunksize=chunk_rows, index_col=idx, dtype={header[idx]: str}):
        values = chunk.to_numpy(dtype=np.float64)
        for k in range(len(chunk)):
            n = new[row]
            row += 1
            if not n:
                continue
            if count[n] == 1:
                X[out_row[n]] = values[k]
            else:
                pending[n].append(values[k])
                if len(pending[n]) == count[n]:
                    X[out_row[n]] = merge_rows(pending.pop(n), is_log2, agg)
    assert not pending, "incomplete duplicate groups: %s" % list(pending)[:5]
    report = pd.DataFrame({"input_id": ids_in, "output_id": new, "status": status})
    return list(out_row), samples, X, report
