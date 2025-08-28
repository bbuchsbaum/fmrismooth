# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

fmrismooth is an R package providing fast edge-preserving 3D/4D smoothing and denoising for fMRI data. It implements advanced filtering techniques including:
- MP-PCA (Marchenko-Pastur PCA) denoising
- Total Variation (TV) denoising with robust loss functions (Huber, Tukey)
- Bilateral and joint bilateral filters using permutohedral lattice approximation
- Spatio-temporal guided filters
- Variance Stabilizing Transforms (VST)

## Build and Development Commands

### Building the package
```bash
# Install package dependencies and build
R CMD build .
R CMD INSTALL fmrismooth_*.tar.gz

# Or from R:
devtools::load_all()  # Load for development
devtools::install()   # Install locally
```

### Documentation
```bash
# Generate documentation from roxygen2 comments
devtools::document()
```

### Testing
No test framework currently detected. To add tests:
```bash
# Set up testthat
usethis::use_testthat()
# Then run tests with:
devtools::test()
```

## Architecture

### Core Components

1. **C++ Backend** (`src/`): High-performance implementations using Rcpp
   - `bilateral_lattice.cpp`: Permutohedral lattice for fast bilateral filtering
   - `mp_pca4d.cpp`: Marchenko-Pastur PCA denoising
   - `robust_tv_huber_4d.cpp`: Robust TV smoothing with Huber/Tukey loss
   - `st_guided_filter.cpp`: Spatio-temporal guided filtering
   - `tv_denoise_4d.cpp`: Total variation denoising

2. **R Interface** (`R/`):
   - **Main pipelines**: 
     - `default_pipeline.R`: One-liner interface with automatic parameter selection
     - `pipeline.R`: Alternative pipeline implementation
   - **Core algorithms**:
     - `mp_pca4d.R`: MP-PCA denoising wrapper
     - `joint_bilateral.R`: Joint bilateral filtering with multiple guides
     - `bilateral_lattice.R`: Lattice-based bilateral filtering
     - `st_robust_smooth4d.R`: Robust spatio-temporal smoothing
   - **Utilities**:
     - `auto_params.R`: Automatic parameter recommendation based on data characteristics
     - `utils_io.R`, `utils_compat.R`: neuroim2 integration for neuroimaging data
     - `resample_align.R`, `utils_resample.R`: Spatial alignment and resampling

### Key Design Patterns

- **neuroim2 Integration**: Optional dependency for neuroimaging data structures (NeuroVol, NeuroVec)
- **Flexible Input**: Functions accept both plain arrays and neuroim2 objects
- **Automatic Resampling**: Guides and masks are automatically aligned to fMRI grid when neuroim2 is available
- **Parameter Auto-tuning**: `recommend_params()` suggests optimal parameters based on voxel size and target smoothness

### Typical Workflow

1. Load fMRI data as 4D array or neuroim2::NeuroVec
2. Optionally provide anatomical T1 guide and tissue probability maps
3. Call `fmrismooth_default()` with automatic parameter selection
4. Or build custom pipeline combining MP-PCA → robust smoothing → bilateral filtering

## Dependencies

- **Required**: Rcpp, RcppArmadillo (for C++ linear algebra)
- **Optional**: neuroim2 (for neuroimaging data handling and spatial resampling)
- **System**: BLAS/LAPACK libraries (handled by R installation)