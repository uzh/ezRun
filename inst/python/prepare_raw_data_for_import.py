#!/usr/bin/env python3
"""
Prepare raw data for import.

This script prepares raw data files for import into the analysis pipeline.

Usage:
    python prepare_raw_data_for_import.py <input_file> -d <directory> [-o <output_dir>]

Arguments:
    input_file                  Input file to process
    -d, --directory            Base directory for output
    -o, --output              Output directory (optional, defaults to --directory value)
"""

import sys
import argparse
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def prepare_data(input_path: Path, output_dir: Path) -> None:
    """
    Prepare raw data for import.

    Parameters
    ----------
    input_path : Path
        Input file path
    output_dir : Path
        Output directory path
    """
    logger.info(f"Input: {input_path}")
    logger.info(f"Output directory: {output_dir}")
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # TODO: Add data preparation logic here
    logger.info("Data preparation complete!")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare raw data for import into the analysis pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Using -d/--directory flag (output defaults to directory)
    python prepare_raw_data_for_import.py input.txt -d /path/to/output
    
    # Specifying separate output directory
    python prepare_raw_data_for_import.py input.txt -d /path/to/base -o /path/to/output
        """
    )
    
    parser.add_argument(
        'input',
        type=str,
        help='Input file path'
    )
    
    parser.add_argument(
        '-d', '--directory',
        type=str,
        required=True,
        help='Base directory for output (used as default for output directory)'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        help='Output directory (optional, defaults to --directory value if not specified)'
    )
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    
    # Determine output directory: use -o if provided, otherwise use -d
    if args.output:
        output_dir = Path(args.output)
    else:
        output_dir = Path(args.directory)
    
    if not input_path.exists():
        logger.error(f"Input file not found: {input_path}")
        sys.exit(1)
    
    logger.info("=" * 50)
    prepare_data(input_path, output_dir)
    logger.info("=" * 50)
    logger.info("Done!")


if __name__ == '__main__':
    main()
