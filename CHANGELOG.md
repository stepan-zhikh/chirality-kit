# Changelog

All notable changes to chirality-kit are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.8.0] - 2026-11-17

Tube-end-wall release: holed graphene sheets at the two ends of the tube, and a
`side` key that places each solvent group in its own part of the box.

The graphene walls (box ends, 0.7.0, and tube ends, 0.8.0) are still being
tested in production simulations. Builds without walls are unchanged; check a
walled build before relying on it, and report anything unexpected through
[GitHub Issues](https://github.com/stepan-zhikh/chirality-kit/issues).

### Added

- `system.tube-sheets`: a holed graphene sheet at each end of the tube, on its
  own or together with the box-end sheets. Each sheet lies in the plane of its
  end's innermost rim carbon (a tube carbon bonded to only two other tube
  carbons), so on a chiral tube the whole rim reaches the sheet or passes it.
  The hole is centred on the tube axis and is R + `hole-buffer` wide; carbons
  left with fewer than two neighbours are removed too. The tube passes through
  the hole and is not bonded to the sheet.
- `sheets` block keys `hole-buffer` (default 2.5 Å) and the residue names
  `tube-z-min` and `tube-z-max` (`MEMN` and `MEMP` by default). Both kinds of
  sheet read the one block, so all four sheets share one lattice and one x and y
  size.
- `rim-offset` in the `sheets` block (default 0 Å) moves both tube-end sheets
  along the tube from the innermost rim carbons: positive outward, as far as the
  outermost rim carbon, negative inward, keeping the sheets more than
  2 × `wall-gap` apart. The log and the `.summary` report the rim carbons and
  the sheet planes, and the log warns when a sheet carbon comes within 2 Å of
  the tube or its groups.
- The annulus is cleared when the system is built: every solvent molecule with
  an atom between the tube-end sheets outside the tube, or within `wall-gap` of
  a sheet carbon, is removed and listed in `<system-name>.carved`. Ions are kept
  out of that slab when they are placed, so none is deleted and the charge is
  unchanged.
- `restraints.ref` holds the tube-end sheets with the tube; the box-end sheets
  stay free.
- `side` for every solvation group, in the same terms as the terminal groups'
  `side`: `"full"` (the whole box, the default), `"-"` or `"+"` (one side of the
  tube), `"both"` (the two sides, nothing around or inside the tube) or
  `"around-tube"`. The boundaries are the tube-end
  sheet planes, kept `wall-gap` clear, or the z range of the `TUBE` residue
  without tube-end sheets. A `concentration` counts the region's volume.
- The `[Sheets]` summary reports each tube-end sheet (atoms, sites removed for
  the hole, edge carbons), the hole radius, the closest approach to the tube and
  the carved residues; the PBC check reports the two reservoirs. Builds write
  `<system-name>_tube_sheet_zmin.xyz` and `_zmax.xyz`, and auto system names gain
  `_TUBESHEETS`.
- Development tests for the tube-end walls: the rim planes of chiral, armchair
  and zigzag tubes, the hole and its edge, the side boxes and the ion slab, the
  new validation, on/off builds, Drude tube-end sheets, one lattice for all
  sheets, layout agreement, each `side` region, and the bonded topology of a
  tube-end sheet.
- Tutorial Section 7.4, Walls at the tube ends, and one table for the `sheets`
  keys of both kinds of wall; `side` in the tutorial's solvation keys.
- Tutorial figures in `docs/tutorial/`: the three tube families at one scale,
  functional groups and vacant sites, a solvated tube and one-sided salt,
  box-end walls, and tube-end walls at three `rim-offset` values. Each is built
  from a JSON input in `docs/tutorial/inputs/` or a shipped example, with seed
  12345, using CHARMM36m.

### Changed

- `VACANT` names an empty site in terminal and ring placement. Inputs that use
  `DUMMY`, the name before 0.8.0, are read the same way and build the same
  system; the log notes the new name. The example
  `chirality_kit_input_dummy_gaps_dry.json` is now
  `chirality_kit_input_vacant_sites_dry.json`.
- The box-end and tube-end sheets go into the main file through one shared
  step, with the same output for box-end builds.
- The residue-name checks cover the tube-end names too; the names of the sheets
  in use must all differ.

## [0.7.0] - 2026-11-17

Graphene-wall release: periodic graphene sheets that close a system at its box
ends, so the tube is the only route between the two reservoirs.

### Added

- `Sheet_Builder`: periodic graphene sheets that close a built system at its
  z-min and z-max. Turn them on with `system.sheets` and set them in an optional
  `sheets` block (orientation, wall-gap, seam-gap, residue names). The sheets
  use the tube's carbon species, including its Drude particle in polarisable
  force fields, are bonded across the x and y boundaries, and are each their own
  residue (`SHTN` and `SHTP` by default).
- A `sheets` block fits the box x and y to whole graphene cells before solvation,
  raising the xy padding only as far as needed. With the same block and seed, a
  build with `system.sheets` on matches the build with it off in every non-sheet
  atom.
- Builds with a `sheets` block write the exact box: a NAMD `.xsc`, a `CRYST1`
  line in the `.pdb`, and the box in the `.xdata`/`.data` header. The summary
  gains a `[Sheets]` section.
- `_inject_xyz_atoms_and_bonds` takes `periodic_lengths` to find bonds across
  periodic boundaries with the minimum image.
- Input validation for the sheets keywords: a periodic tube takes neither
  `system.sheets` nor a `sheets` block, and the sheet residue names also stay
  clear of the segment name `SYS`, the ion residue names (`POT`, `CLA`, ...) and
  the string `END`.
- Development tests for the walls: sheets on vs off with Drude water and salt (identical
  non-sheet atoms, one Drude per sheet carbon, unchanged lone pairs and
  anisotropy, no charge added), the `armchair-x` orientation, custom gaps and
  residue names, an explicit box, one box across `.data`, `.xsc` and `CRYST1`,
  `restraints.ref` against the `.psf`, the PBC check on a sheets build, and a
  per-sheet topology check (3N/2 bonds of degree 3, 3N angles, 6N dihedrals, no
  impropers, one core-Drude bond per carbon).
- Tutorial Section 7, Graphene walls: the keys, the box rounded up to whole
  graphene cells, and the outputs. `conventions.md` describes how the sheets
  reuse the `Nanotube-Species` entry and which residue names they must avoid.

### Changed

- The PBC box check measures the `TUBE` residue, the tube with its groups and
  caps as the `.pdb` and `.psf` name it, so sheets, solvent and ions stay out
  of the tube extent.
- Every build starts from fresh sheet state, so one IO object can run several
  builds. The `IO` class also declares that state at its defaults, so an `IO`
  made without `__init__` reads "no sheets".
- A faster bonded-term step. Once the atoms and bonds are placed,
  chirality-kit derives every angle, dihedral and improper from the bond graph
  and looks up its parameters in the force field's Coeffs. This step now reads
  the atom types once per call and keeps each Coeffs match it finds, and atom
  and bond type codes are looked up once per type. The output is byte-identical
  to 0.6.4: with seed 12345, `examples/chirality_kit_input_vacant_sites_dry.json`
  builds in 2 s (140 s in 0.6.4). Systems with many bonded terms, such as the
  graphene sheets, gain the most.
- Tutorial sections 7–11 are now 8–12.

## [0.6.4] - 2026-10-13

Bonded-topology release: complete COO- topology for Drude-polarisable tubes,
ring-carboxyl torsions for COOH, and a cleaner `.psf` dihedral section. The
same input and seed produce byte-identical `.pdb` files to 0.6.3, and
byte-identical `.psf` files for additive CHARMM36m builds without COOH and for
bare tubes in either force field.

### Added

- Drude (`charmm36m-ndp`) tubes with COO- now include the COO- planarity
  improper, `CD2O2A OD2C2A OD2C2A CD2R6A`, in the `.psf` output.
- Drude tubes with COO- now include the COO--C--C--H dihedrals to the
  neighbouring rim hydrogens (`HDR6A`), so Drude and additive tubes carry the
  same set of bonded terms.
- CHARMM36m COOH groups now include the ring-carboxyl torsions
  `CG2R61 CG2R61 CG2O2 OG2D1` and `CG2R61 CG2R61 CG2O2 OG311`, with the CGenFF
  values (1.025 kcal/mol, n = 2, 180 degrees).

### Changed

- The `.psf` writer lists each dihedral once, including multi-term (Fourier)
  dihedrals such as the four-term COOH C--C--O--H series. NAMD and CHARMM apply
  every parameter-file term for the type quartet to that single entry, so the
  `.psf` follows the CHARMM convention whatever the order of the terms. The
  `.xdata` file keeps one entry per term. NAMD energies for COOH structures are
  unchanged.

### Validation

- Every built structure was validated internally against its bond network:
  each angle and dihedral implied by the bonds, wherever the force field has a
  coefficient for it, appears exactly once in the `.psf`; each COO- and COOH
  carbon carries one centre-first improper; and Drude and additive tubes share
  one topology. Covered: 0, 2, 3 and 12 COO- per end in both force fields,
  3 COOH per end, in both the package and standalone layouts. The full
  development test suite passes, confirming no regression elsewhere.

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
