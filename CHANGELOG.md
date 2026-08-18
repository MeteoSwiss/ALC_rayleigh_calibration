# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!-- This changelog was seeded from the project's current state; entries for
     versions before 2.0.0 are summarized. Please refine dates and details, and
     add an entry to [Unreleased] with every user-facing change. -->

## [Unreleased]

### Added
- `LICENSE` (MIT), matching the license declared in `pyproject.toml`.
- GitHub Actions CI (`.github/workflows/CI_test.yaml`): ruff, mypy and pytest on
  pushes to `main` and on pull requests.
- Community-health files: `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `AUTHORS`,
  and GitHub issue / pull-request templates.

### Changed
- `requires-python`, ruff `target-version` and mypy `python_version` aligned to
  the version actually used and deployed (3.12).

## [2.0.0]
<!-- release date: TBD -->

The v2 calibration line: a read-once pipeline performing Rayleigh, liquid-cloud
and Kalman calibration for the E-PROFILE automated lidar/ceilometer network.

### Added
- Rayleigh (molecular) calibration against a CAMS-derived molecular profile.
- Liquid-cloud calibration (O'Connor method) with PVC (Hogan 2006) multiple-
  scattering η tables, per Vaisala instrument type.
- Water-vapour transmission correction from CAMS model levels — mandatory for
  910 nm instruments (a no-WV mode is rejected).
- Observation-minus-background (OmB) and sensitivity diagnostics.
- Out-of-domain CAMS routing via per-station lean boxes for affiliates outside
  the Europe+Arctic domain.
- Static monitoring dashboard (`monitoring/` + `scripts/build_dashboard.py`).

### Fixed
- Klett inversion sign error (still present in the operational v1.0 line).

## [1.x]

Earlier operational E-PROFILE calibration. See the git history for detail.

[Unreleased]: https://github.com/MeteoSwiss/ALC_rayleigh_calibration/compare/v2.0.0...HEAD
