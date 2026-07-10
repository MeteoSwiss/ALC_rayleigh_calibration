# Contributing

Thank you for your interest in contributing to the **E-PROFILE ALC calibration**
pipeline — the Rayleigh and liquid-cloud calibration for automated lidars and
ceilometers (ALC) operated across the EUMETNET E-PROFILE network. Contributions,
bug reports and questions are welcome.

## Code of Conduct

This project adheres to the [Contributor Covenant](CODE_OF_CONDUCT.md) Code of
Conduct. By participating, you are expected to uphold it.

## How to contribute

### Reporting bugs

Please check existing issues first to avoid duplicates, then open a bug report
using the **Bug report** issue template. A good report includes:

- what you did, what you expected, and what actually happened;
- a minimal way to reproduce it (a code snippet or the failing command);
- your environment: Python version, OS, and — where relevant — the instrument
  type (CHM15k, CL31, CL51, CL61, MiniMPL/MPL) and station.

### Suggesting features

Open an issue with the **Feature request** template describing the proposed
change and why it would be useful.

### Contributing code

1. Fork the repository and create a branch for your change.
2. Make your changes, following the conventions below.
3. Add or update tests for the behaviour you change.
4. Make sure the test suite and the quality checks pass locally.
5. Open a pull request using the PR template.

## Development setup

Prerequisites: Python 3.10+ (the pipeline is developed and deployed on 3.12) and
git.

```console
$ git clone https://github.com/MeteoSwiss/ALC_rayleigh_calibration.git
$ cd ALC_rayleigh_calibration
$ python -m venv .venv
$ source .venv/bin/activate        # Windows: .venv\Scripts\activate
$ pip install -e ".[dev,dashboard]"
```

Add the `download` extra (`pip install -e ".[dev,dashboard,download]"`) only if
you work on the CAMS auto-download path — it pulls `cdsapi` and `cfgrib`.

## Running the tests

```console
$ pytest
```

Many tests are **integration tests that self-skip** when their inputs are not
present — the full CAMS archive, HITRAN look-up tables, ADS credentials, or the
MATLAB parity fixtures. In a clean checkout they skip automatically and the
remaining unit tests run against the bundled sample data under `examples/data/`.
This is expected; a green run with skips is a pass.

## Code style and quality

The project uses [ruff](https://docs.astral.sh/ruff/) for linting and
formatting, and [mypy](https://mypy-lang.org/) for type checking. Configuration
lives in [`pyproject.toml`](pyproject.toml).

```console
$ ruff check calibration monitoring      # lint
$ ruff format calibration monitoring     # auto-format
$ mypy calibration                       # type check
```

New function definitions should carry type annotations (`disallow_untyped_defs`
is enabled). Prefer clear, readable code and comment the *why*, not just the
*what*. These checks also run in CI (`.github/workflows/CI_test.yaml`); they are
currently advisory and will become blocking once the tree is fully clean.

## Pull request process

- Reference any related issue (e.g. "Closes #123").
- Keep the change focused; unrelated cleanups belong in their own PR.
- If this is your first contribution, add your name to the [`AUTHORS`](AUTHORS)
  file.
- Make sure CI passes.
- Do not commit data files, credentials (`~/.cdsapirc`, tokens) or large
  binaries. See [`.gitignore`](.gitignore) for what is deliberately excluded.

### Commit messages

Use the imperative mood ("Add …", "Fix …"), keep the summary line short, and add
a body explaining the reasoning when the change is non-trivial.

## License

By contributing, you agree that your contributions will be licensed under the
[MIT License](LICENSE).

## Questions

Operational context for maintainers lives in `doc/OPERATIONS.md` and `CLAUDE.md`.
For anything else, open an issue.
