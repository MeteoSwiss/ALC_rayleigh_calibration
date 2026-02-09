#!/usr/bin/env python3
"""
Automated Rayleigh Calibration Runner.

This script runs Rayleigh calibration for all configured instruments.
It's designed to be run as a cron job or scheduled task.

Usage
-----
    python -m rayleigh_calibration.main [OPTIONS]

Options
-------
    --date YYYYMMDD     Process specific date (default: yesterday)
    --config PATH       Path to options JSON file
    --instruments PATH  Path to instruments JSON file
    --output PATH       Override output directory
    --verbose          Enable verbose logging
    --dry-run          Don't write output files

Example
-------
    python -m rayleigh_calibration.main --date 20240115 --verbose
"""

from __future__ import annotations

import argparse
import logging
import sys
from datetime import date, datetime, timedelta
from pathlib import Path
from typing import Optional, Tuple, List, Dict

from . import (
    calibrate_rayleigh,
    load_instruments,
    CalibrationOptions,
    CalibrationResult,
    InstrumentInfo,
)


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for the calibration run."""
    level = logging.DEBUG if verbose else logging.INFO
    format_str = "%(asctime)s | %(levelname)-8s | %(name)s | %(message)s"
    
    logging.basicConfig(
        level=level,
        format=format_str,
        handlers=[
            logging.StreamHandler(sys.stdout),
        ]
    )


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Run Rayleigh calibration for ceilometers/lidars",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    parser.add_argument(
        "--date",
        type=str,
        default=None,
        help="Date to process (YYYYMMDD format, default: yesterday)",
    )
    
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("E-PROFILE_options_calibration_rayleigh.json"),
        help="Path to calibration options JSON file",
    )
    
    parser.add_argument(
        "--instruments",
        type=Path,
        default=Path("E-PROFILE_instruments_L2_auto.json"),
        help="Path to instruments configuration JSON file",
    )
    
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Override output directory",
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Don't write output files",
    )
    
    parser.add_argument(
        "--station",
        type=str,
        default=None,
        help="Process only this station (WMO ID)",
    )
    
    return parser.parse_args()


def run_calibration_batch(
    date_str: str,
    instruments: List[InstrumentInfo],
    options: CalibrationOptions,
    station_filter: Optional[str] = None,
) -> Tuple[Dict[str, CalibrationResult], List[str], List[str]]:
    """
    Run calibration for all instruments.
    
    Parameters
    ----------
    date_str : str
        Date string (YYYYMMDD).
    instruments : list
        List of instrument configurations.
    options : CalibrationOptions
        Calibration options.
    station_filter : str, optional
        If provided, only process this WMO ID.
        
    Returns
    -------
    tuple
        (results dict, errors list, not_implemented list)
    """
    logger = logging.getLogger(__name__)
    results: Dict[str, CalibrationResult] = {}
    errors: List[str] = []
    not_implemented: List[str] = []
    
    for info in instruments:
        station_key = f"{info.wmo_id}_{info.identifier}"
        
        # Filter by station if requested
        if station_filter and info.wmo_id != station_filter:
            continue
            
        logger.info("=" * 70)
        logger.info(f"Processing: {info.site_name}")
        logger.info("=" * 70)
        
        # Skip if calibration not requested
        if not info.calibrated:
            logger.info("No calibration requested for this instrument")
            continue
            
        # Check if instrument type supports calibration
        if not info.instrument_type.supports_calibration:
            logger.info(f"Calibration not supported for {info.instrument_type.value}")
            continue
            
        try:
            result = calibrate_rayleigh(date_str, info, options)
            results[station_key] = result
            
            if result.is_successful:
                logger.info(f"✓ Success: CL = {result.lidar_constant:.4e}")
            else:
                logger.warning(f"✗ Failed: {result.flag_meaning}")
                
        except NotImplementedError as e:
            logger.warning(f"Not implemented: {e}")
            not_implemented.append(info.instrument_type.value)
            
        except Exception as e:
            logger.error(f"Error processing {info.site_name}: {e}", exc_info=True)
            errors.append(info.site_name)
            
    return results, errors, not_implemented


def print_summary(
    results: Dict[str, CalibrationResult],
    errors: List[str],
    not_implemented: List[str],
    elapsed_time: float,
) -> None:
    """Print summary of calibration run."""
    logger = logging.getLogger(__name__)
    
    logger.info("")
    logger.info("=" * 70)
    logger.info("CALIBRATION SUMMARY")
    logger.info("=" * 70)
    
    # Count successes and failures
    n_success = sum(1 for r in results.values() if r.is_successful)
    n_failed = sum(1 for r in results.values() if not r.is_successful)
    
    logger.info(f"Processed:     {len(results)}")
    logger.info(f"Successful:    {n_success}")
    logger.info(f"Failed:        {n_failed}")
    logger.info(f"Errors:        {len(errors)}")
    logger.info(f"Elapsed time:  {elapsed_time:.1f}s")
    
    if errors:
        logger.warning(f"Stations with errors: {errors}")
        
    if not_implemented:
        unique_types = list(set(not_implemented))
        logger.warning(f"Unimplemented instrument types: {unique_types}")


def main() -> int:
    """Main entry point."""
    args = parse_args()
    setup_logging(args.verbose)
    
    logger = logging.getLogger(__name__)
    start_time = datetime.now()
    
    # Determine date to process
    if args.date:
        date_str = args.date
    else:
        yesterday = date.today() - timedelta(days=1)
        date_str = yesterday.strftime("%Y%m%d")
        
    logger.info(f"Processing date: {date_str}")
    
    # Load configuration
    try:
        instruments = load_instruments(args.instruments)
        logger.info(f"Loaded {len(instruments)} instruments from {args.instruments}")
    except FileNotFoundError:
        logger.error(f"Instruments file not found: {args.instruments}")
        return 1
        
    try:
        options = CalibrationOptions.from_json(args.config)
        logger.info(f"Loaded options from {args.config}")
    except FileNotFoundError:
        logger.error(f"Config file not found: {args.config}")
        return 1
        
    # Override output directory if specified
    if args.output:
        options.folder_output = args.output
        
    # Run calibration
    results, errors, not_impl = run_calibration_batch(
        date_str=date_str,
        instruments=instruments,
        options=options,
        station_filter=args.station,
    )
    
    # Print summary
    elapsed = (datetime.now() - start_time).total_seconds()
    print_summary(results, errors, not_impl, elapsed)
    
    # Return non-zero if there were errors
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
