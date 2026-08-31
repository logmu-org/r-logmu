# CRAN submission comments — logmu 0.1.0

## Test environments

- Windows 11, R 4.6.0 (local)
- Windows (R-devel), via win-builder
- Windows (R-release), via win-builder
- macOS-latest (R-release), Windows-latest (R-release),
  Ubuntu-latest (R-devel, R-release, R-oldrel-1),
  via GitHub Actions (`r-lib/actions/check-standard`)

## R CMD check results

0 errors | 1 warning | 2 note

* WARNING: "checking pragmas in C/C++ headers and code"

  The flagged pragmas occur only in unmodified, vendored upstream source
  from the EVE (Expressive Vector Engine) header-only SIMD library, under
  `src/eve/`. They are part of EVE's own handling of compiler-specific
  diagnostics around SIMD intrinsics and are not introduced or modified by
  this package. EVE and its dependency SPY are bundled under the Boost Software
  License 1.0; see `LICENSE.note`, `src/eve/LICENSE.txt` and 
  `src/spy/LICENSE.txt` for attribution.

* NOTE: Standard new-submission flag from the CRAN incoming feasibility check.

* NOTE: "checking compilation flags used"

  The non-portable flags (`-mavx2`, `-mfma`, `-mavx512f`, `-mavx512cd`,
  `-mavx512bw`, `-mavx512dq`, `-mavx512vl`) are applied only to a small
  number of isolated translation units (`vec_kernel_avx2.cpp`,
  `vec_kernel_avx512.cpp`) that implement function-multiversioned SIMD
  kernels.

  At runtime the package detects CPU support via CPUID and only dispatches
  to these kernels on hosts whose CPU supports the corresponding
  instruction set (see `cpp_vec_active_tier()` / `cpp_vec_active_lanes()`).
  On all other hosts it falls back to `vec_kernel_baseline.cpp`, which is
  compiled without any of these flags. This is the standard "function
  multiversioning" approach to portable SIMD code, and does not affect the
  portability or safety of the resulting binary.

## Downstream dependencies

There are currently no downstream dependencies for this package.
