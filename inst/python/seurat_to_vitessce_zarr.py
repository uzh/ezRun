#!/usr/bin/env python3
"""
Convert h5ad to Vitessce-optimized Zarr format.

This script is designed to be called after Seurat → h5ad conversion,
producing a Zarr store optimized for Vitessce visualization:
- Small VAR_CHUNK_SIZE (10) for fast gene lookups
- CSC sparse format for efficient column access
- Minimal metadata (only essential obs columns)

Usage:
    python seurat_to_vitessce_zarr.py <input.h5ad> <output_dir>

The script will create:
    <output_dir>/vitessce/data.zarr/
    <output_dir>/vitessce/metadata.json
"""

import sys
import json
import logging
from pathlib import Path

import anndata
import numpy as np
import scipy.sparse as sp

# Vitessce-optimized settings
VAR_CHUNK_SIZE = 10  # Critical for fast gene lookups
OBS_CHUNK_SIZE = None  # Use full obs dimension (all cells)

# Essential metadata columns to keep (others are stripped)
ESSENTIAL_OBS_COLS = [
    # Clustering
    'seurat_clusters', 'banksy_cluster', 'cluster', 'Cluster',
    # Cell type annotations
    'RCTD_Main', 'RCTD_Full', 'predicted.celltype', 'celltype',
    'cell_type', 'CellType', 'annotation',
    # QC metrics
    'nCount_Xenium', 'nFeature_Xenium', 'nCount_RNA', 'nFeature_RNA',
    'cell_area', 'nucleus_area', 'percent.mt',
    # Sample info
    'orig.ident', 'Sample', 'sample', 'condition', 'Condition',
    # Spatial
    'fov', 'FOV',
]

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def optimize_anndata_for_vitessce(adata: anndata.AnnData) -> anndata.AnnData:
    """
    Optimize AnnData object for Vitessce visualization.

    1. Convert expression matrix to CSC sparse format (for column/gene access)
    2. Filter obs columns to essential ones only
    3. Remove unnecessary data (layers, uns, varm)

    Parameters
    ----------
    adata : anndata.AnnData
        Input AnnData object (will be modified in place)

    Returns
    -------
    anndata.AnnData
        Optimized AnnData object
    """
    logger.info(f"Input: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    # 1. Optimize expression matrix format
    if adata.X is not None:
        if sp.issparse(adata.X):
            if not sp.isspmatrix_csc(adata.X):
                logger.info("Converting expression matrix to CSC sparse format...")
                adata.X = sp.csc_matrix(adata.X)
            else:
                logger.info("Expression matrix already in CSC format")
        else:
            logger.info("Expression matrix is dense (keeping as-is)")

    # 2. Filter obs columns to essential ones
    original_cols = list(adata.obs.columns)
    keep_cols = [c for c in ESSENTIAL_OBS_COLS if c in adata.obs.columns]

    # Also keep any numeric columns that might be useful
    for col in adata.obs.columns:
        if col not in keep_cols:
            if adata.obs[col].dtype in [np.float64, np.float32, np.int64, np.int32]:
                if col not in ['nCount_Xenium', 'nFeature_Xenium']:  # Skip already included
                    keep_cols.append(col)

    if keep_cols:
        logger.info(f"Keeping {len(keep_cols)}/{len(original_cols)} obs columns: {keep_cols}")
        adata.obs = adata.obs[keep_cols].copy()
    else:
        logger.warning("No essential columns found, keeping all obs columns")

    # 3. Keep only gene names in var (strip other annotations)
    adata.var = adata.var[[]]  # Empty DataFrame, keeps index (gene names)

    # 4. Clear unnecessary data structures
    adata.uns = {}
    adata.varm = {}
    adata.layers = {}
    # Keep obsm (embeddings like X_umap, spatial) - these are essential!

    # 5. Log obsm contents
    if adata.obsm:
        logger.info(f"Available embeddings (obsm): {list(adata.obsm.keys())}")

    return adata


def write_optimized_zarr(adata: anndata.AnnData, zarr_path: str) -> None:
    """
    Write AnnData to Zarr with Vitessce-optimized chunking.

    Uses VAR_CHUNK_SIZE=10 for fast gene lookups (load 10 genes per request).
    This is critical for Vitessce performance when selecting genes.

    Parameters
    ----------
    adata : anndata.AnnData
        Optimized AnnData object
    zarr_path : str
        Output path for Zarr store
    """
    import shutil

    # Remove existing Zarr if present
    zarr_path = Path(zarr_path)
    if zarr_path.exists():
        logger.info(f"Removing existing Zarr: {zarr_path}")
        shutil.rmtree(zarr_path)

    # Calculate optimal chunks
    n_obs, n_var = adata.shape
    obs_chunk = OBS_CHUNK_SIZE if OBS_CHUNK_SIZE else n_obs
    var_chunk = min(VAR_CHUNK_SIZE, n_var)
    chunks = (obs_chunk, var_chunk)

    logger.info(f"Writing Zarr with chunks: {chunks}")
    logger.info(f"  - {n_obs:,} cells in {n_obs // obs_chunk + (1 if n_obs % obs_chunk else 0)} chunk(s)")
    logger.info(f"  - {n_var:,} genes in {n_var // var_chunk + (1 if n_var % var_chunk else 0)} chunk(s)")

    # Write Zarr
    adata.write_zarr(str(zarr_path), chunks=chunks)

    # Calculate approximate size
    zarr_size = sum(f.stat().st_size for f in zarr_path.rglob('*') if f.is_file())
    logger.info(f"Zarr size: {zarr_size / 1024 / 1024:.1f} MB")


def write_metadata(adata: anndata.AnnData, output_dir: Path) -> None:
    """
    Write metadata JSON for Shiny app fast-path detection.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object
    output_dir : Path
        Directory containing the Zarr store
    """
    metadata = {
        'n_cells': int(adata.shape[0]),
        'n_genes': int(adata.shape[1]),
        'available_embeddings': list(adata.obsm.keys()) if adata.obsm else [],
        'available_cell_sets': list(adata.obs.columns),
        'var_chunk_size': VAR_CHUNK_SIZE,
        'is_sparse': sp.issparse(adata.X) if adata.X is not None else False,
    }

    metadata_path = output_dir / 'metadata.json'
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)

    logger.info(f"Metadata written to: {metadata_path}")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        print("Error: Missing arguments")
        print(f"Usage: {sys.argv[0]} <input.h5ad> <output_dir>")
        sys.exit(1)

    h5ad_path = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])

    if not h5ad_path.exists():
        logger.error(f"Input file not found: {h5ad_path}")
        sys.exit(1)

    # Create output directory structure
    vitessce_dir = output_dir / 'vitessce'
    vitessce_dir.mkdir(parents=True, exist_ok=True)
    zarr_path = vitessce_dir / 'data.zarr'

    logger.info(f"Input: {h5ad_path}")
    logger.info(f"Output: {zarr_path}")
    logger.info("=" * 50)

    # Load h5ad
    logger.info("Loading h5ad...")
    adata = anndata.read_h5ad(h5ad_path)

    # Optimize for Vitessce
    logger.info("Optimizing for Vitessce...")
    adata = optimize_anndata_for_vitessce(adata)

    # Write optimized Zarr
    logger.info("Writing optimized Zarr...")
    write_optimized_zarr(adata, zarr_path)

    # Write metadata
    write_metadata(adata, vitessce_dir)

    logger.info("=" * 50)
    logger.info("Done!")


if __name__ == '__main__':
    main()
