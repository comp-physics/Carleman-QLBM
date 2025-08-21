# GitHub Actions CI/CD Workflows

This directory contains automated testing workflows for the CLBM (Carleman Lattice Boltzmann Method) project.

## Workflows

### 1. **Main CI (`ci.yml`)**
- **Triggers**: Push/PR to main/master/develop branches
- **Tests**: Julia 1.9 & 1.10 on Ubuntu & macOS
- **Coverage**: 
  - Sparse vs Dense matrix tests (adaptive on macOS)
  - Main CLBM simulation (full on Linux, minimal on macOS)
  - Unit tests (comprehensive on Linux, minimal on macOS)
  - Integration tests
- **macOS Compatibility**: Uses `--break-system-packages` and adaptive testing to handle Python environment restrictions

### 2. **Quick Tests (`quick-test.yml`)**
- **Triggers**: Changes to `src/CLBM/` files
- **Purpose**: Fast feedback for development
- **Tests**:
  - Syntax validation
  - Minimal sparse vs dense test
  - Function loading verification

### 3. **Performance Tests (`performance-tests.yml`)**
- **Triggers**: Push/PR to main/master + weekly schedule
- **Purpose**: Ensure performance across different problem sizes
- **Tests**:
  - Multiple truncation orders (2, 3, 4)
  - Different matrix sizes
  - Memory usage validation

### 4. **Code Quality (`code-quality.yml`)**
- **Triggers**: Push/PR to any branch
- **Purpose**: Maintain code quality standards
- **Checks**:
  - Syntax validation
  - Common issue detection
  - Memory efficiency tests
  - Documentation coverage

## Local Testing

Before pushing, you can run tests locally:

```bash
# Run comprehensive tests
cd src/CLBM
julia test_sparse_vs_dense.jl

# Run unit tests
julia unit_tests.jl

# Run main simulation
julia clbm_run.jl
```

## Test Structure

### Unit Tests (`unit_tests.jl`)
- Configuration validation
- Sparse Kronecker function tests
- Matrix dimension checks
- LBM setup verification
- Matrix construction tests

### Integration Tests (`test_sparse_vs_dense.jl`)
- Mathematical equivalence
- Performance comparison
- Memory usage analysis
- Sparsity validation

## Requirements

The CI workflows install these Julia packages automatically:
- `SymPy`
- `PyPlot` 
- `HDF5`
- `LaTeXStrings`
- `SparseArrays`
- `LinearAlgebra`
- `Test`

## Environment Variables

Required environment variables:
- `QCFD_HOME`: Repository root path
- `QCFD_SRC`: Source directory path (`$QCFD_HOME/src/`)

## macOS Compatibility

### Python Environment Issues
Modern macOS systems use externally-managed Python environments (PEP 668) which prevent direct `pip install` commands. Our workflows handle this by:

1. **Using `--break-system-packages`**: Safe for CI environments
2. **Adaptive testing**: Falls back to minimal tests if plotting packages fail
3. **Platform-specific logic**: Different test strategies for Linux vs macOS

### Test Adaptation
- **Linux**: Full comprehensive tests with all plotting capabilities
- **macOS**: Adaptive tests that gracefully handle matplotlib issues
- **Fallback**: Minimal unit tests ensure core functionality always works

## Badge Status

Add these badges to your main README.md:

```markdown
[![CI](https://github.com/YOUR_USERNAME/YOUR_REPO/workflows/CI/badge.svg)](https://github.com/YOUR_USERNAME/YOUR_REPO/actions)
[![Performance Tests](https://github.com/YOUR_USERNAME/YOUR_REPO/workflows/Performance%20Tests/badge.svg)](https://github.com/YOUR_USERNAME/YOUR_REPO/actions)
[![Code Quality](https://github.com/YOUR_USERNAME/YOUR_REPO/workflows/Code%20Quality/badge.svg)](https://github.com/YOUR_USERNAME/YOUR_REPO/actions)
```
