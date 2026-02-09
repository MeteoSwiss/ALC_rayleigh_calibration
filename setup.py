#!/usr/bin/env python3
"""
Setup script for rayleigh-calibration package.

This provides a fallback installation method when Poetry is not available.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
readme_path = Path(__file__).parent / "README.md"
long_description = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""

setup(
    name="rayleigh-calibration",
    version="2.0.0",
    author="E-PROFILE",
    author_email="contact@e-profile.eu",
    description="Rayleigh calibration tools for automated lidars and ceilometers",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://e-profile.eu",
    packages=find_packages(exclude=["tests", "tests.*"]),
    python_requires=">=3.9",
    install_requires=[
        "numpy>=1.20",
        "scipy>=1.7",
        "netCDF4>=1.5",
    ],
    extras_require={
        "plotting": ["matplotlib>=3.4"],
        "dev": ["pytest>=7.0", "pytest-cov>=4.0"],
    },
    entry_points={
        "console_scripts": [
            "rayleigh-calibration=rayleigh_calibration.main:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
    ],
)
