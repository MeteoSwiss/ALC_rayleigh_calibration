#!/usr/bin/env python3
"""
Automated Rayleigh Calibration Runner.

This script runs Rayleigh calibration for all configured instruments.
It's designed to be run as a cron job or scheduled task.

Usage
-----
    python -m calibration.main [OPTIONS]

Options
-------
    --date YYYYMMDD [YYYYMMDD]  Process one date or a date range (default: yesterday)
    --config PATH               Path to options JSON file
    --instruments PATH          Path to instruments JSON file
    --output PATH               Override output directory
    --verbose                   Enable verbose logging
    --dry-run                   Don't write output files

Example
-------
    python -m calibration.main --date 20240115 --verbose
    python -m calibration.main --date 20240101 20240131 --verbose
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
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

    # Suppress noisy third-party loggers
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    logging.getLogger("PIL").setLevel(logging.WARNING)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Run Rayleigh calibration for ceilometers/lidars",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    parser.add_argument(
        "--date",
        type=str,
        nargs='+',
        default=None,
        help="Date(s) to process: one YYYYMMDD or two for a range (default: yesterday)",
    )
    
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("options.json"),
        help="Path to calibration options JSON file",
    )
    
    parser.add_argument(
        "--instruments",
        type=Path,
        default=Path("instruments.json"),
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
    
    parser.add_argument(
        "--workers", "-j",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1, 0 = number of CPUs)",
    )
    
    return parser.parse_args()


def _calibrate_one(
    date_str: str,
    info: InstrumentInfo,
    options: CalibrationOptions,
) -> Tuple[str, str, Optional[CalibrationResult], Optional[str]]:
    """
    Calibrate a single instrument (designed to run in a worker process).

    Returns
    -------
    tuple
        (station_key, site_name, result_or_None, error_type_or_None)
        error_type is 'not_implemented', 'error', or None on success.
    """
    station_key = f"{info.wmo_id}_{info.identifier}"
    try:
        result = calibrate_rayleigh(date_str, info, options)
        return (station_key, info.site_name, result, None)
    except NotImplementedError:
        return (station_key, info.site_name, None, "not_implemented")
    except Exception as e:
        # Log the full traceback inside the worker process
        logging.getLogger(__name__).error(
            f"Error processing {info.site_name}: {e}", exc_info=True
        )
        return (station_key, info.site_name, None, "error")


def run_calibration_batch(
    date_str: str,
    instruments: List[InstrumentInfo],
    options: CalibrationOptions,
    station_filter: Optional[str] = None,
    max_workers: int = 1,
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
    max_workers : int
        Number of parallel workers (1 = sequential).
        
    Returns
    -------
    tuple
        (results dict, errors list, not_implemented list)
    """
    logger = logging.getLogger(__name__)
    results: Dict[str, CalibrationResult] = {}
    errors: List[str] = []
    not_implemented: List[str] = []

    # Build the list of instruments to process
    to_process: List[InstrumentInfo] = []
    for info in instruments:
        if station_filter and info.wmo_id != station_filter:
            continue
        if not info.calibrated:
            logger.info(f"Skipping {info.site_name}: no calibration requested")
            continue
        if not info.instrument_type.supports_calibration:
            logger.info(f"Skipping {info.site_name}: {info.instrument_type.value} not supported")
            continue
        to_process.append(info)

    logger.info(f"{len(to_process)} instruments to calibrate")

    # --- Sequential path (workers=1): keeps all logging in order ----------
    if max_workers == 1:
        for info in to_process:
            logger.info("=" * 70)
            logger.info(f"Processing: {info.site_name}")
            logger.info("=" * 70)
            key, name, result, err_type = _calibrate_one(date_str, info, options)
            if err_type == "not_implemented":
                logger.warning(f"Not implemented for {info.instrument_type.value}")
                not_implemented.append(info.instrument_type.value)
            elif err_type == "error":
                errors.append(name)
            else:
                assert result is not None
                results[key] = result
                if result.is_successful:
                    logger.info(f"✓ Success: CL = {result.lidar_constant:.4e}")
                else:
                    logger.warning(f"✗ Failed: {result.flag_meaning}")
        return results, errors, not_implemented

    # --- Parallel path ----------------------------------------------------
    logger.info(f"Using {max_workers} parallel workers")

    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(_calibrate_one, date_str, info, options): info
            for info in to_process
        }
        for future in as_completed(futures):
            info = futures[future]
            key, name, result, err_type = future.result()
            if err_type == "not_implemented":
                logger.warning(f"Not implemented: {info.instrument_type.value}")
                not_implemented.append(info.instrument_type.value)
            elif err_type == "error":
                errors.append(name)
            else:
                assert result is not None
                results[key] = result
                if result.is_successful:
                    logger.info(f"✓ {name}: CL = {result.lidar_constant:.4e}")
                else:
                    logger.warning(f"✗ {name}: {result.flag_meaning}")

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


def date_range(start_str: str, end_str: str) -> List[str]:
    """
    Generate a list of YYYYMMDD date strings from *start_str* to *end_str*
    (inclusive).
    """
    start = datetime.strptime(start_str, "%Y%m%d").date()
    end = datetime.strptime(end_str, "%Y%m%d").date()
    if end < start:
        raise ValueError(
            f"End date ({end_str}) must not be before start date ({start_str})"
        )
    dates: List[str] = []
    current = start
    while current <= end:
        dates.append(current.strftime("%Y%m%d"))
        current += timedelta(days=1)
    return dates


def main() -> int:
    """Main entry point."""
    args = parse_args()
    setup_logging(args.verbose)
    
    logger = logging.getLogger(__name__)
    start_time = datetime.now()
    
    # Determine date(s) to process
    if args.date is None:
        yesterday = date.today() - timedelta(days=1)
        dates_to_process = [yesterday.strftime("%Y%m%d")]
    elif len(args.date) == 1:
        dates_to_process = [args.date[0]]
    elif len(args.date) == 2:
        try:
            dates_to_process = date_range(args.date[0], args.date[1])
        except ValueError as e:
            logger.error(str(e))
            return 1
    else:
        logger.error("--date accepts at most two values: START_DATE [END_DATE]")
        return 1

    logger.info(
        f"Processing {len(dates_to_process)} date(s): "
        f"{dates_to_process[0]} → {dates_to_process[-1]}"
    )
    
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

    # Resolve worker count
    max_workers = args.workers
    if max_workers == 0:
        max_workers = os.cpu_count() or 1
    logger.info(f"Workers: {max_workers}")

    # ---- Loop over all requested dates ----
    all_results: Dict[str, CalibrationResult] = {}
    all_errors: List[str] = []
    all_not_impl: List[str] = []

    for date_str in dates_to_process:
        logger.info("")
        logger.info("#" * 70)
        logger.info(f"# DATE: {date_str}")
        logger.info("#" * 70)

        results, errors, not_impl = run_calibration_batch(
            date_str=date_str,
            instruments=instruments,
            options=options,
            station_filter=args.station,
            max_workers=max_workers,
        )

        # Prefix keys with date so multi-day results don't collide
        for key, result in results.items():
            all_results[f"{date_str}/{key}"] = result
        all_errors.extend(errors)
        all_not_impl.extend(not_impl)
    
    # Print overall summary
    elapsed = (datetime.now() - start_time).total_seconds()
    print_summary(all_results, all_errors, all_not_impl, elapsed)
    
    # Return non-zero if there were errors
    return 1 if all_errors else 0


if __name__ == "__main__":
    sys.exit(main())
