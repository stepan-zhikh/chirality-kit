from __future__ import annotations

import math
import os
from typing import TYPE_CHECKING

import numpy as np

from ._runtime import my_log_file


# Each sheet is its own residue, numbered after everything already in the box
# (the tube is residue 42; solvent and ions follow it).
SHEET_SIDES = ("z-min", "z-max")
TUBE_SHEET_SIDES = ("tube-z-min", "tube-z-max")
SHEET_ORIENTATIONS = ("zigzag-x", "armchair-x")
SHEET_DEFAULTS = {
    "lattice-xy": True,
    "orientation": "zigzag-x",
    "wall-gap": 3.0,
    "seam-gap": 15.0,
    "hole-buffer": 2.5,
    "rim-offset": 0.0,
    "resnames": {"z-min": "SHTN", "z-max": "SHTP", "tube-z-min": "MEMN", "tube-z-max": "MEMP"},
}
TUBE_RESIDUE = 42


class Sheet_Builder:

    '''
    Periodic graphene sheets that close a built system: at its two z ends, at
    the two ends of the tube, or both.

    Box-end sheets are added after the tube, its functional groups and the
    solvent are in place. One sheet sits wall-gap below the lowest atom and one
    wall-gap above the highest; nothing already in the box is moved or deleted.

    Tube-end sheets sit at the plane of each end's innermost rim carbon, moved
    along the axis by rim-offset, and have
    a round hole on the tube axis, R + hole-buffer wide, so the tube passes
    through without being bonded to them. The solvent in the annulus between
    them is removed before they go in.

    Every sheet is bonded across the x and y periodic boundaries, so the box x
    and y must be whole graphene cells: lattice_box_xy scales them up to fit
    before any solvent is placed. The sheet carbons use the nanotube's own
    carbon species (and its Drude particle in polarisable force fields), and
    each sheet is its own residue with its own residue name.
    '''

    def __init__(self,
                 force_field_manager: Force_Field_Manager,
                 file_manager: File_Manager,
                 config: Config,
                 geometry: Geometry) -> None:

        self.force_field_manager : Force_Field_Manager = force_field_manager
        self.file_manager : File_Manager = file_manager
        self.config : Config = config
        self.gt : Geometry = geometry

    # ============================================================
    # LATTICE
    # ============================================================

    @staticmethod
    def resolve_config(sheet_cfg):

        '''
        Merge a JSON "sheets" block over SHEET_DEFAULTS, including the nested
        resnames mapping.

        Parameters:
        ----------
        sheet_cfg : dict or None
            The user's "sheets" block; None gives the defaults.

        Returns:
        -------
        resolved : dict
            Every key of SHEET_DEFAULTS, user values first.

        Raises:
        ------
        None
        '''

        resolved = {key: value for key, value in SHEET_DEFAULTS.items() if key != "resnames"}
        resolved["resnames"] = dict(SHEET_DEFAULTS["resnames"])
        for key, value in (sheet_cfg or {}).items():
            if key == "resnames" and isinstance(value, dict):
                resolved["resnames"].update(value)
            else:
                resolved[key] = value
        return resolved

    @staticmethod
    def lattice_periods(bond_length, orientation):

        '''
        Repeat lengths of the rectangular 4-atom graphene cell along x and y.

        Parameters:
        ----------
        bond_length : float
            C-C bond length b in Angstrom.
        orientation : str
            'zigzag-x' puts the zigzag direction (period sqrt(3) b) along x and
            the armchair direction (period 3 b) along y; 'armchair-x' swaps them.

        Returns:
        -------
        (period_x, period_y) : tuple of float

        Raises:
        ------
        ValueError
            If the orientation is not one of SHEET_ORIENTATIONS.
        '''

        zigzag = math.sqrt(3.0) * bond_length
        armchair = 3.0 * bond_length
        if orientation == "zigzag-x":
            return zigzag, armchair
        if orientation == "armchair-x":
            return armchair, zigzag
        raise ValueError(f"Unsupported sheet orientation '{orientation}'; "
                         f"choose from {', '.join(SHEET_ORIENTATIONS)}.")

    def lattice_box_xy(self, x_range, y_range, bond_length, orientation):

        '''
        Scale an x and a y box range up to whole graphene cells about their own
        centres. A range already a whole number of cells keeps its length; no
        range is ever made smaller.

        Parameters:
        ----------
        x_range, y_range : sequence of two floats
            [lo, hi] of the requested box along x and y.
        bond_length : float
            C-C bond length b in Angstrom.
        orientation : str
            See lattice_periods.

        Returns:
        -------
        fitted : dict
            'x' and 'y' ([lo, hi]), 'nx' and 'ny' (cell counts), 'period_x' and
            'period_y' (Angstrom).

        Raises:
        ------
        ValueError
            If the orientation is unsupported.
        '''

        period_x, period_y = self.lattice_periods(bond_length, orientation)

        def fit(lo, hi, period):
            length = float(hi) - float(lo)
            # The small tolerance keeps an exact multiple at its own length.
            cells = max(1, math.ceil(length / period - 1e-9))
            centre = (float(lo) + float(hi)) / 2.0
            half = cells * period / 2.0
            return [centre - half, centre + half], cells

        x_fit, nx = fit(x_range[0], x_range[1], period_x)
        y_fit, ny = fit(y_range[0], y_range[1], period_y)
        return {
            "x": x_fit, "y": y_fit, "nx": nx, "ny": ny,
            "period_x": period_x, "period_y": period_y,
        }

    def sheet_positions(self, x_range, y_range, z, bond_length, orientation):

        '''
        Positions of one complete periodic graphene sheet filling the given x and
        y ranges at height z. The rectangular cell holds four atoms, offset by a
        quarter cell so no atom sits on a box face; with whole cells every atom
        has exactly three neighbours at bond_length under the minimum image.

        Parameters:
        ----------
        x_range, y_range : sequence of two floats
            [lo, hi] of the box along x and y; each must be a whole number of
            cells (see lattice_box_xy).
        z : float
            Height of the sheet plane.
        bond_length : float
            C-C bond length b in Angstrom.
        orientation : str
            See lattice_periods.

        Returns:
        -------
        positions : (4*nx*ny, 3) np.ndarray

        Raises:
        ------
        ValueError
            If a range is not a whole number of cells.
        '''

        b = float(bond_length)
        period_x, period_y = self.lattice_periods(b, orientation)
        lx = float(x_range[1]) - float(x_range[0])
        ly = float(y_range[1]) - float(y_range[0])
        nx = int(round(lx / period_x))
        ny = int(round(ly / period_y))
        if nx < 1 or ny < 1 or abs(nx * period_x - lx) > 1e-6 or abs(ny * period_y - ly) > 1e-6:
            raise ValueError(
                f"Sheet box {lx:.6f} x {ly:.6f} Ang is not a whole number of "
                f"{period_x:.6f} x {period_y:.6f} Ang graphene cells."
            )

        a = math.sqrt(3.0) * b
        # Zigzag along the first axis (period a), armchair along the second (period 3 b).
        basis = [(0.0, 0.0), (a / 2.0, b / 2.0), (a / 2.0, 1.5 * b), (0.0, 2.0 * b)]
        offset = (a / 4.0, b / 4.0)
        cells_1, cells_2 = (nx, ny) if orientation == "zigzag-x" else (ny, nx)

        positions = []
        for i in range(cells_1):
            for j in range(cells_2):
                for u, v in basis:
                    p1 = offset[0] + i * a + u
                    p2 = offset[1] + j * 3.0 * b + v
                    x, y = (p1, p2) if orientation == "zigzag-x" else (p2, p1)
                    positions.append([float(x_range[0]) + x, float(y_range[0]) + y, float(z)])
        return np.array(positions, dtype=float)

    # ============================================================
    # BUILD
    # ============================================================

    def _write_sheet_xyz(self, positions, element, filename, title):

        with open(filename, "w") as f:
            f.write(f"{len(positions)}\n")
            f.write(f"{title}\n")
            for x, y, z in positions:
                f.write(f"{element:<3s} {float(x):>12.5f} {float(y):>12.5f} {float(z):>12.5f}\n")
        return filename

    def _inject_sheet(self, sites, resname, xyz, title, lx, ly, expected_bonds):

        '''
        Write one sheet's sites to xyz and add them to the main file as a new
        residue: the atoms (with Drude particles in polarisable force fields),
        the C-C bonds found under the minimum image in x and y, and the angles,
        dihedrals and impropers those bonds imply.

        Parameters:
        ----------
        sites : (N, 3) np.ndarray
            Sheet carbon positions.
        resname : str
            Residue name of the sheet.
        xyz : str
            Path of the sheet's .xyz file.
        title : str
            Comment line of the .xyz file.
        lx, ly : float
            Box lengths along x and y, the periods of the sheet.
        expected_bonds : int
            C-C bonds the sheet must have; any other count means the sites and
            the box do not match.

        Returns:
        -------
        residue_idx : int
            The sheet's residue index.
        bonds_data : list
            The injected C-C bonds.

        Raises:
        ------
        ValueError
            If the bond count differs from expected_bonds.
        '''

        primary = self.config.default_species_primary
        atom_type = primary["atom_type"]
        element = primary["atom_element"]
        bond_cutoff = self.config.get_default_species_bond_param("bond_cutoff", 1.60)
        type_charges = {
            entry["atom_type"]: entry["atom_charge"]
            for entry in self.config.default_species_chain
        }

        # The next free residue index, so no solvent residue is renamed.
        residue_idx = int(self.file_manager.find_last_residue()) + 1
        self._write_sheet_xyz(sites, element, xyz, title)

        bonds_data = self.force_field_manager._inject_xyz_atoms_and_bonds(
            xyz_filename=xyz,
            element_to_atom_type={element: atom_type},
            allowed_bonds={(atom_type, atom_type): bond_cutoff},
            type_charges=type_charges,
            new_residue_idx=str(residue_idx),
            atom_index_offset=None,
            include_drude=self.config.drude_polarisable,
            periodic_lengths=(lx, ly, None),
        )
        # bonds_data is empty when no bonded section is written at all.
        bonded = bool(bonds_data) or "Bonds" in self.config.enabled_xdata_sections
        if bonded and len(bonds_data) != expected_bonds:
            raise ValueError(
                f"{resname} sheet: expected {expected_bonds} periodic C-C bonds, "
                f"found {len(bonds_data)}; the box is not lattice matched."
            )
        self.force_field_manager.inject_angles_from_bonds(bonds_data)
        self.force_field_manager.inject_dihedrals_from_bonds(bonds_data)
        self.force_field_manager.inject_impropers_from_bonds(bonds_data)
        self.config.residue_name_by_index[residue_idx] = resname
        return residue_idx, bonds_data

    def add_sheets(self, box_size, sheet_cfg):

        '''
        Add the z-min and z-max sheets to the main file and return the extended
        box. Uses no random numbers, so everything built before it is unchanged.

        Parameters:
        ----------
        box_size : list of three [lo, hi]
            Current box; its x and y must already be whole graphene cells.
        sheet_cfg : dict
            Resolved "sheets" block (see resolve_config).

        Returns:
        -------
        report : dict
            'box' (the new box, z extended to hold the sheets plus seam-gap),
            'z' (sheet heights by side), 'atoms_per_sheet', 'drudes_per_sheet',
            'nx', 'ny', 'charge' (net charge of both sheets), 'resnames',
            'residues' (residue index by side) and 'atom_type'.

        Raises:
        ------
        ValueError
            If the main file has no atoms, the box is not lattice matched, or the
            sheets carry a net charge.
        '''

        primary = self.config.default_species_primary
        atom_type = primary["atom_type"]
        bond_length = self.config.get_default_species_bond_param("bond_length", 1.44)
        orientation = sheet_cfg["orientation"]
        wall_gap = float(sheet_cfg["wall-gap"])
        seam_gap = float(sheet_cfg["seam-gap"])

        x_range = [float(box_size[0][0]), float(box_size[0][1])]
        y_range = [float(box_size[1][0]), float(box_size[1][1])]
        lx, ly = x_range[1] - x_range[0], y_range[1] - y_range[0]

        positions, _, _, _ = self.force_field_manager.get_atom_positions_all()
        if len(positions) == 0:
            raise ValueError("Sheets need a built system; the main file has no atoms.")
        z_low = float(np.min(positions[:, 2]))
        z_high = float(np.max(positions[:, 2]))
        sheet_z = {"z-min": z_low - wall_gap, "z-max": z_high + wall_gap}

        charge_before = self.force_field_manager.get_system_charge()

        counts = {}
        residues = {}
        for side in SHEET_SIDES:
            resname = sheet_cfg["resnames"][side]
            sheet = self.sheet_positions(x_range, y_range, sheet_z[side], bond_length, orientation)
            tag = side.replace("-", "")
            xyz = os.path.join(self.config.system_folder,
                               f"{self.config.system_name}_sheet_{tag}.xyz")
            residue_idx, bonds_data = self._inject_sheet(
                sheet, resname, xyz, f"chirality-kit {resname} sheet at z = {sheet_z[side]:.4f}",
                lx, ly, 3 * len(sheet) // 2,
            )
            residues[side] = residue_idx
            counts[side] = len(sheet)
            my_log_file.info(
                f"Added {resname} sheet: {len(sheet)} {atom_type} atoms at z = "
                f"{sheet_z[side]:.4f} Ang, {len(bonds_data)} periodic bonds."
            )

        sheet_charge = self.force_field_manager._round_charge(
            self.force_field_manager.get_system_charge() - charge_before)
        if abs(sheet_charge) > 1e-6:
            raise ValueError(f"The sheets carry a net charge of {sheet_charge:+.6f} e; they must be neutral.")

        new_box = [
            x_range,
            y_range,
            [sheet_z["z-min"] - seam_gap / 2.0, sheet_z["z-max"] + seam_gap / 2.0],
        ]
        period_x, period_y = self.lattice_periods(bond_length, orientation)
        return {
            "box": new_box,
            "z": sheet_z,
            "atoms_per_sheet": counts,
            "drudes_per_sheet": {side: (counts[side] if self.config.drude_polarisable else 0)
                                 for side in SHEET_SIDES},
            "nx": int(round(lx / period_x)),
            "ny": int(round(ly / period_y)),
            "charge": sheet_charge,
            "resnames": dict(sheet_cfg["resnames"]),
            "residues": residues,
            "atom_type": atom_type,
        }

    # ============================================================
    # TUBE-END SHEETS
    # ============================================================

    @staticmethod
    def _minimum_image_xy(delta, lx, ly):

        delta = np.array(delta, dtype=float, copy=True)
        delta[..., 0] -= lx * np.round(delta[..., 0] / lx)
        delta[..., 1] -= ly * np.round(delta[..., 1] / ly)
        return delta

    def _within(self, points, sites, lx, ly, cutoff):

        '''
        True for each point closer than cutoff to any site, with the minimum
        image in x and y. Works in blocks, so large systems stay in memory.
        '''

        mask = np.zeros(len(points), dtype=bool)
        for start in range(0, len(points), 512):
            delta = self._minimum_image_xy(points[start:start + 512, None, :] - sites[None, :, :], lx, ly)
            mask[start:start + 512] = (np.linalg.norm(delta, axis=-1) < cutoff).any(axis=1)
        return mask

    def _closest_distance(self, sites, atoms, lx, ly):

        closest = float("inf")
        for start in range(0, len(sites), 512):
            delta = self._minimum_image_xy(sites[start:start + 512, None, :] - atoms[None, :, :], lx, ly)
            closest = min(closest, float(np.linalg.norm(delta, axis=-1).min()))
        return closest

    def tube_frame(self, rim_offset=0.0, min_separation=0.0):

        '''
        Locate the tube in the built structure, before any solvent: its axis,
        its radius, the plane of each end's innermost rim carbon, and the z
        range of the TUBE residue (the tube with its groups and caps).

        The tube carbons are the residue-42 atoms of the primary nanotube
        species that sit on the tube wall, within 0.5 Ang of the median radius,
        so ring or group carbons of the same type are left out. A rim carbon
        has only two tube-carbon neighbours. On a chiral tube the rim is
        staggered; each end's plane is its innermost rim carbon, so the whole
        rim reaches the sheet or passes it. rim_offset moves both planes along
        the axis: positive outward, towards the box ends, negative inward.

        Parameters:
        ----------
        rim_offset : float, optional
            Shift of each sheet plane from its innermost rim carbon, in
            Angstrom; positive is outward. Default 0.0.
        min_separation : float, optional
            The two planes must stay more than this far apart. Default 0.0.

        Returns:
        -------
        frame : dict
            'axis' (x, y), 'radius', 'z_lo' and 'z_hi' (the sheet planes of the
            -z and +z ends), 'rim_lo' and 'rim_hi' (the innermost rim carbons),
            'rim_offset', 'tube_z' (z range of the TUBE residue), 'atoms'
            (positions of the TUBE residue), 'carbons' (the tube carbon
            positions) and 'rim' (rim carbons per end).

        Raises:
        ------
        ValueError
            If there are no tube carbons or an end has no rim carbon, or if
            rim_offset moves a sheet past its end's outermost rim carbon or
            leaves the two planes min_separation or less apart.
        '''

        atom_type = self.config.default_species_primary["atom_type"]
        bond_cutoff = self.config.get_default_species_bond_param("bond_cutoff", 1.60)

        positions, _, types, residues = self.force_field_manager.get_atom_positions_all()
        positions = np.asarray(positions, dtype=float)
        in_tube = np.asarray(residues) == TUBE_RESIDUE
        if not in_tube.any():
            raise ValueError("Tube-end sheets need a built tube; the TUBE residue is empty.")
        tube_z = (float(positions[in_tube, 2].min()), float(positions[in_tube, 2].max()))

        carbons = positions[in_tube & (np.asarray(types) == atom_type)]
        if len(carbons) < 3:
            raise ValueError(f"Tube-end sheets need the tube's {atom_type} carbons; found {len(carbons)}.")
        radius = np.hypot(*(carbons[:, :2] - carbons[:, :2].mean(axis=0)).T)
        carbons = carbons[np.abs(radius - np.median(radius)) < 0.5]
        axis = carbons[:, :2].mean(axis=0)
        radius = float(np.hypot(*(carbons[:, :2] - axis).T).mean())

        degree = np.zeros(len(carbons), dtype=int)
        for start in range(0, len(carbons), 512):
            dist = np.linalg.norm(carbons[start:start + 512, None, :] - carbons[None, :, :], axis=-1)
            degree[start:start + 512] = ((dist > 1e-6) & (dist <= bond_cutoff)).sum(axis=1)
        z_mid = 0.5 * (carbons[:, 2].min() + carbons[:, 2].max())
        rim = degree == 2
        top = carbons[rim & (carbons[:, 2] > z_mid), 2]
        bottom = carbons[rim & (carbons[:, 2] < z_mid), 2]
        if top.size == 0 or bottom.size == 0:
            raise ValueError("Tube-end sheets need a rim at both tube ends: no carbon with only "
                             "two tube-carbon neighbours was found at one end.")

        # Each sheet plane starts at its end's innermost rim carbon; rim_offset
        # moves both outward (positive) or inward (negative) along the axis.
        rim_lo, rim_hi = float(bottom.max()), float(top.min())
        z_lo, z_hi = rim_lo - rim_offset, rim_hi + rim_offset
        if rim_offset > 0:
            # Past the outermost rim carbon the tube no longer passes through the hole.
            reach = min(float(top.max()) - rim_hi, rim_lo - float(bottom.min()))
            if rim_offset > reach + 1e-9:
                raise ValueError(
                    f"rim-offset {rim_offset:g} Ang moves a tube-end sheet past the outermost rim "
                    f"carbon, so the tube would no longer pass through the hole; this tube allows "
                    f"at most {reach:.3f} Ang outward.")
        if (rim_offset or min_separation) and z_hi - z_lo <= min_separation:
            raise ValueError(
                f"rim-offset {rim_offset:g} Ang leaves the tube-end sheets {z_hi - z_lo:.3f} Ang apart; "
                f"they must stay more than {min_separation:.3f} Ang (2 x wall-gap) apart.")

        return {
            "axis": (float(axis[0]), float(axis[1])),
            "radius": radius,
            "z_lo": z_lo,
            "z_hi": z_hi,
            "rim_lo": rim_lo,
            "rim_hi": rim_hi,
            "rim_offset": float(rim_offset),
            "tube_z": tube_z,
            "atoms": positions[in_tube],
            "carbons": carbons,
            "rim": {"tube-z-min": int(bottom.size), "tube-z-max": int(top.size)},
        }

    def holed_sheet_positions(self, x_range, y_range, z, bond_length, orientation, axis, hole_radius):

        '''
        One periodic graphene sheet at height z with a round hole on the tube
        axis, drilled the way Qin & Zhou (2024) drill theirs: every site closer
        to the axis than hole_radius (minimum image in x and y) is removed, then
        any site left with fewer than two neighbours, until none is left, so no
        carbon hangs by a single bond.

        Parameters:
        ----------
        x_range, y_range, z, bond_length, orientation
            As for sheet_positions.
        axis : (float, float)
            x and y of the tube axis.
        hole_radius : float
            Radius of the hole in Ang: the tube radius plus hole-buffer.

        Returns:
        -------
        sheet : dict
            'sites' (the kept positions, in sheet_positions order), 'removed'
            (sites taken out, edge pruning included), 'bonds' (C-C bonds among
            the kept sites) and 'edge' (kept carbons with two neighbours).

        Raises:
        ------
        ValueError
            If the hole takes the whole sheet.
        '''

        sites = self.sheet_positions(x_range, y_range, z, bond_length, orientation)
        lx = float(x_range[1]) - float(x_range[0])
        ly = float(y_range[1]) - float(y_range[0])

        delta = self._minimum_image_xy(sites[:, :2] - np.asarray(axis, dtype=float), lx, ly)
        keep = np.hypot(delta[:, 0], delta[:, 1]) >= float(hole_radius)

        cutoff = 1.1 * float(bond_length)
        neighbours = []
        for start in range(0, len(sites), 512):
            delta = self._minimum_image_xy(sites[start:start + 512, None, :] - sites[None, :, :], lx, ly)
            dist = np.linalg.norm(delta, axis=-1)
            neighbours += [np.nonzero((row > 1e-6) & (row <= cutoff))[0] for row in dist]

        while True:
            degree = np.array([int(keep[nb].sum()) for nb in neighbours])
            loose = keep & (degree < 2)
            if not loose.any():
                break
            keep &= ~loose
        if not keep.any():
            raise ValueError(f"A hole of radius {float(hole_radius):.3f} Ang takes the whole sheet.")

        return {
            "sites": sites[keep],
            "removed": int((~keep).sum()),
            "bonds": int(degree[keep].sum()) // 2,
            "edge": int((keep & (degree == 2)).sum()),
        }

    def tube_sheet_sites(self, box_size, sheet_cfg, frame):

        '''
        The two tube-end sheets before they are added: each a holed sheet at
        its end's rim plane, with a hole R + hole-buffer wide on the tube axis.
        Pure geometry, so the solvent they would enclose or overlap can be
        removed first.

        Parameters:
        ----------
        box_size : list of three [lo, hi]
            Current box; its x and y must be whole graphene cells.
        sheet_cfg : dict
            Resolved "sheets" block.
        frame : dict
            From tube_frame.

        Returns:
        -------
        sheets : dict
            By side ('tube-z-min', 'tube-z-max'): the holed_sheet_positions
            result plus 'z'. Also 'hole_radius', 'lx' and 'ly'.

        Raises:
        ------
        ValueError
            If a hole takes the whole sheet.
        '''

        bond_length = self.config.get_default_species_bond_param("bond_length", 1.44)
        x_range = [float(box_size[0][0]), float(box_size[0][1])]
        y_range = [float(box_size[1][0]), float(box_size[1][1])]
        hole_radius = frame["radius"] + float(sheet_cfg["hole-buffer"])
        planes = {"tube-z-min": frame["z_lo"], "tube-z-max": frame["z_hi"]}

        sheets = {"hole_radius": hole_radius,
                  "lx": x_range[1] - x_range[0], "ly": y_range[1] - y_range[0]}
        for side in TUBE_SHEET_SIDES:
            sheet = self.holed_sheet_positions(x_range, y_range, planes[side], bond_length,
                                               sheet_cfg["orientation"], frame["axis"], hole_radius)
            sheet["z"] = planes[side]
            sheets[side] = sheet
        return sheets

    def annulus_residues(self, sheets, frame, wall_gap, keep_residues=()):

        '''
        The residues the tube-end sheets would enclose or overlap. A residue
        outside the TUBE residue is listed when one of its atoms lies between
        the two sheet planes and outside the tube (r > R), unless it is in
        keep_residues (the around-tube solvent groups), or when one of its
        atoms is within wall_gap of a tube-end sheet carbon.

        Parameters:
        ----------
        sheets : dict
            From tube_sheet_sites.
        frame : dict
            From tube_frame.
        wall_gap : float
            Clearance kept around every sheet carbon, in Ang.
        keep_residues : collection of int
            Residues allowed in the annulus.

        Returns:
        -------
        carved : list of dict
            One per residue, in residue order: 'residue', 'type' (its first
            atom's type), 'atoms' (atom indexes), 'ion' (True if any atom is an
            ion or an ion's Drude particle) and 'z' (mean z).

        Raises:
        ------
        None
        '''

        positions, indexes, types, residues = self.force_field_manager.get_atom_positions_all()
        positions = np.asarray(positions, dtype=float)
        residues = np.asarray(residues)
        lx, ly = sheets["lx"], sheets["ly"]
        z_lo, z_hi = sheets["tube-z-min"]["z"], sheets["tube-z-max"]["z"]
        solvent = residues != TUBE_RESIDUE

        delta = self._minimum_image_xy(positions[:, :2] - np.asarray(frame["axis"], dtype=float), lx, ly)
        r = np.hypot(delta[:, 0], delta[:, 1])
        hit = (positions[:, 2] >= z_lo) & (positions[:, 2] <= z_hi) & (r > frame["radius"])
        if keep_residues:
            hit &= ~np.isin(residues, list(keep_residues))

        for side in TUBE_SHEET_SIDES:
            sheet = sheets[side]
            near = np.nonzero(solvent & (np.abs(positions[:, 2] - sheet["z"]) < wall_gap))[0]
            if near.size:
                hit[near] |= self._within(positions[near], sheet["sites"], lx, ly, wall_gap)
        hit &= solvent

        ion_types = set(self.config.ions) | set(self.config.dions) | set(self.config.dion_particles)
        carved = []
        for residue in sorted(set(residues[hit].tolist())):
            members = np.nonzero(residues == residue)[0]
            member_types = [types[i] for i in members]
            carved.append({
                "residue": int(residue),
                "type": member_types[0],
                "atoms": [int(indexes[i]) for i in members],
                "ion": any(t in ion_types for t in member_types),
                "z": float(positions[members, 2].mean()),
            })
        return carved

    def add_tube_sheets(self, sheet_cfg, sheets, frame):

        '''
        Add the two tube-end sheets from tube_sheet_sites to the main file. The
        solvent in the annulus must already be gone. Uses no random numbers.

        Parameters:
        ----------
        sheet_cfg : dict
            Resolved "sheets" block.
        sheets : dict
            From tube_sheet_sites.
        frame : dict
            From tube_frame.

        Returns:
        -------
        report : dict
            By side: 'z', 'atoms_per_sheet', 'drudes_per_sheet', 'removed',
            'edge', 'bonds' and 'residues'. Also 'hole_radius',
            'min_tube_distance' (the closest any sheet carbon comes to an atom
            of the TUBE residue), 'charge', 'resnames' and 'atom_type'.

        Raises:
        ------
        ValueError
            If a sheet's bond count is off or the sheets carry a net charge.
        '''

        atom_type = self.config.default_species_primary["atom_type"]
        lx, ly = sheets["lx"], sheets["ly"]
        charge_before = self.force_field_manager.get_system_charge()

        report = {key: {} for key in ("z", "atoms_per_sheet", "drudes_per_sheet",
                                      "removed", "edge", "bonds", "residues")}
        closest = float("inf")
        for side in TUBE_SHEET_SIDES:
            sheet = sheets[side]
            resname = sheet_cfg["resnames"][side]
            tag = "zmin" if side == "tube-z-min" else "zmax"
            xyz = os.path.join(self.config.system_folder,
                               f"{self.config.system_name}_tube_sheet_{tag}.xyz")
            residue_idx, bonds_data = self._inject_sheet(
                sheet["sites"], resname, xyz,
                f"chirality-kit {resname} tube-end sheet at z = {sheet['z']:.4f}, "
                f"hole radius {sheets['hole_radius']:.4f}",
                lx, ly, sheet["bonds"],
            )
            count = len(sheet["sites"])
            report["z"][side] = sheet["z"]
            report["atoms_per_sheet"][side] = count
            report["drudes_per_sheet"][side] = count if self.config.drude_polarisable else 0
            report["removed"][side] = sheet["removed"]
            report["edge"][side] = sheet["edge"]
            report["bonds"][side] = sheet["bonds"]
            report["residues"][side] = residue_idx
            closest = min(closest, self._closest_distance(sheet["sites"], frame["atoms"], lx, ly))
            my_log_file.info(
                f"Added {resname} tube-end sheet: {count} {atom_type} atoms at z = {sheet['z']:.4f} Ang, "
                f"hole radius {sheets['hole_radius']:.4f} Ang, {sheet['removed']} sites removed, "
                f"{sheet['edge']} edge carbons, {len(bonds_data)} bonds."
            )

        if closest < 2.0:
            my_log_file.warning(
                f"A tube-end sheet carbon is {closest:.3f} Ang from an atom of the TUBE residue; "
                f"check rim-offset and hole-buffer against the rim and its groups.")

        charge = self.force_field_manager._round_charge(
            self.force_field_manager.get_system_charge() - charge_before)
        if abs(charge) > 1e-6:
            raise ValueError(f"The tube-end sheets carry a net charge of {charge:+.6f} e; they must be neutral.")

        report.update({
            "hole_radius": sheets["hole_radius"],
            "min_tube_distance": closest,
            "charge": charge,
            "resnames": {side: sheet_cfg["resnames"][side] for side in TUBE_SHEET_SIDES},
            "atom_type": atom_type,
        })
        return report
