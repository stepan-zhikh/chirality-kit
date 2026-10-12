# Changelog

All notable changes to chirality-kit are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.6.3] - 2026-10-12

Metadata, documentation and force-field data release. Apart from Drude
systems containing Na+, Li+ or Ca2+, the same input and seed produce
byte-identical `.pdb` and `.psf` files to 0.6.2.

### Changed

- GitHub Issues is now the primary support route; the email address remains as
  a fallback.
- The third-party material notice moved from `LICENSE` to a new `NOTICE` file,
  so `LICENSE` is the plain MIT text and is detected as MIT.

### Fixed

- The repository URL in `CITATION.cff` and the README clone command now point
  to `github.com/stepan-zhikh/chirality-kit`.
- The force-field README now names the `system.force-field` input key correctly.
- The third-party notice points to `src/force-fields/LICENSE.charmm` instead of
  a per-directory licence file that does not exist.
- In the Drude force field (`charmm36m-ndp`), Na+, Li+ and Ca2+ now carry their
  core-Drude bond. Their Drude particle used to be written without a bond in the
  `.psf` and `.data` output. NAMD was unaffected, because it assigns Drude
  springs from the polarisability, but tools that take the spring from the bond
  list saw a free Drude particle.
- The built-in default input in both layouts recorded `"version": "1.0"` and
  the old filename `chirality_main.py` in its metadata, which runs without
  `--json` archived. It now records `0.6.3` and `chirality_kit.py`.

## [0.6.2] - 2026-09-14

Reference-structure validation release.

### Added

- A Structure Validation section documenting armchair, zigzag, and chiral
  carbon nanotubes against peer-reviewed structural values.
- Separate validation results for additive CHARMM36m and Drude-polarisable
  carbon nanotubes, with Drude comparisons restricted to carbon core sites.

### Internal validation

- Development-only tests generate unsolvated, unfunctionalised nanotubes with
  periodicity disabled and validate atom counts, bonding, coordination,
  diameters, chiral angles, axial lengths, and C-C spacing.
- The internal tests pass for both the standalone and split implementations;
  the manually maintained README tables record the corresponding results.

## [0.6.1] - 2026-09-07

First public release.

### Added

- A basic roadmap section in the README.
