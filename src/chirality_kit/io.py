from __future__ import annotations

import datetime
import json
import logging
import os
import re
import shutil
import string
import time
from collections import defaultdict
from itertools import combinations, permutations
from typing import TYPE_CHECKING

import numpy as np

from ._runtime import my_log_file, set_logger
from .functional_group_generator import (LEGACY_VACANT_FG_TYPE, TERM_GROUP_FIELDS, TERM_SIDES,
                                         VACANT_FG_TYPE, VACANT_FG_TYPES)
from .sheet_builder import SHEET_DEFAULTS, SHEET_ORIENTATIONS, SHEET_SIDES, TUBE_SHEET_SIDES, Sheet_Builder


XDATA_SECTION_SETTINGS = {
    "bonds": "Bonds",
    "angles": "Angles",
    "dihedrals": "Dihedrals",
    "impropers": "Impropers",
    "lone-pairs": "Lone Pairs",
    "anisotropy": "Anisotropy",
}


class IO: 

    # Sheet state, at the values __init__ and every build start from. Declared on
    # the class too, so an IO made without __init__ (test shells) still reads
    # "no sheets" instead of failing on a missing attribute.
    _sheet_cfg = None
    _sheet_lattice = None
    _sheet_report = None
    _tube_frame = None
    _tube_sheet_report = None
    _carved_count = 0

    def __init__(self,
                 force_field_manager: Force_Field_Manager,
                 file_manager: File_Manager,
                 structure: Structure_Generator, 
                 solvation: Solution_Generator,
                 functionalisation: Functional_Group_Generator,
                 config: Config,
                 geometry: Geometry,
                 json_input: str | dict | None = None,
                 sheets: Sheet_Builder | None = None) -> None:

        self.force_field_manager : Force_Field_Manager = force_field_manager
        self.file_manager : File_Manager = file_manager

        self.structure : Structure_Generator = structure
        self.solvation : Solution_Generator = solvation
        self.functionalisation : Functional_Group_Generator = functionalisation
        self.sheets : Sheet_Builder = (
            sheets if sheets is not None
            else Sheet_Builder(force_field_manager, file_manager, config, geometry)
        )

        self.geometry : Geometry = geometry
        self.config : Config = config

        self.json_input : str | dict | None = json_input

        self.box_size = None
        self._psf_charge_neutrality_report = None
        self._sheet_cfg = None      # resolved "sheets" block, None when absent
        self._sheet_lattice = None  # lattice fit of the box x and y
        self._sheet_report = None   # what add_sheets built, None when sheets are off
        self._tube_frame = None     # tube axis, radius and rim planes, when needed
        self._tube_sheet_report = None  # what add_tube_sheets built, None when off
        self._around_residues = set()   # residues of around-tube solvent groups
        self._carved_count = 0          # solvent residues removed for the tube-end sheets

        self.initialise_LOG(json_input)
        self.initialise_JSON(json_input)
        # Seed here so CLI and library runs use the same path.
        self.config.apply_seed(self.json_input)
        self._disabled_xdata_sections = self._get_disabled_xdata_sections()
        self.config.enabled_xdata_sections = (
            set(XDATA_SECTION_SETTINGS.values()) - self._disabled_xdata_sections
        )
        self.config.functionalisation_requested = bool(
            self.json_input["system"]["functionalisation"]
        )
        self.initialise_OUT_FOLDER(json_input)

    # ============================
    # INITIALISATION 
    # ============================

    def initialise_JSON(self, json_input): 

        '''
        Load the simulation settings into self.json_input from one of three
        sources: a dict passed directly (e.g. the module-level INPUT constant),
        a file path string pointing to a JSON file, or None which causes the
        function to look for 'chirality_kit_input.json' next to the script.

        Parameters:
        ----------
        json_input : str, dict, or None
            If a dict, used directly. If a str, treated as a file path. If None,
            the default 'chirality_kit_input.json' in the script directory is read.

        Returns:
        -------
        None

        Raises:
        ------
        FileNotFoundError
            If json_input is None and the default JSON file does not exist.
        ValueError
            If the file at the given path is not valid JSON.
        '''

        if isinstance(json_input, dict): # Dict passed directly (e.g. the INPUT constant)
            self.json_input = json_input

        elif json_input is None: # Normal JSON location
            json_file_path = os.path.join(self.config.script_location, 'chirality_kit_input.json')
            try:
                with open(json_file_path, 'r') as f:
                    self.json_input = json.load(f)  # parse JSON directly from file
            except FileNotFoundError:
                my_log_file.error("Input JSON file not found.")
                raise
            except json.JSONDecodeError:
                my_log_file.error("Input JSON file is not valid JSON.")
                raise ValueError("Input JSON file is not valid JSON.")
            
        else: # Override JSON location (file path string)
            try:
                with open(json_input, 'r') as f:
                    self.json_input = json.load(f)  # parse JSON directly from file
            except Exception as e:
                my_log_file.error(f"Failed to load JSON input: {e}")
                raise

    def initialise_LOG(self, json_input):

        '''
        Provide json simulation setting to the logger, for future support.

        Parameters:
        ----------
        json_input : str, dict, or None
            If a dict, used directly. If a str, treated as a file path. If None,
            the default 'chirality_kit_input.json' in the script directory is read.

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''
    
        set_logger(logging.getLogger(__name__))
        my_log_file.info(f"Log file for chirality_kit run on {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        my_log_file.info(f"System name: {self.config.system_name}\n")

    def initialise_OUT_FOLDER(self, json_input):

        '''
        Create the output directory, set up the Python logging handler to write
        to a log file inside that directory, and create a blank xdata file with
        all section headers. Sets config.system_folder and config.main_file.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        if self.config.out_name:
            base_dir = os.path.expanduser(self.config.out_name)
            if not os.path.isabs(base_dir):
                base_dir = os.path.join(self.config.called_from, base_dir)
            base_dir = os.path.abspath(base_dir)
            if not os.path.basename(base_dir):
                raise ValueError("--out-name must identify a folder name.")
        else:
            system_name = self.json_input['settings']['system_name']
            if system_name == "auto":
                system_name = self.create_system_name()
            base_dir = os.path.join(self.config.called_from, system_name)

        increment = 0
        while True:
            proposed_dir = base_dir if increment == 0 else f"{base_dir}_({increment})"
            try:
                os.makedirs(proposed_dir)
                break
            except FileExistsError:
                increment += 1

        self.config.system_folder = proposed_dir
        self.config.system_name = os.path.basename(proposed_dir)

        log_file_path = os.path.join(self.config.system_folder, self.config.system_name + '.log') # setup the logger
        logging.basicConfig(
            filename=log_file_path,
            filemode="a",                    
            format="%(asctime)s - %(levelname)s - %(message)s",
            level=logging.DEBUG              
        )

        with open(os.path.join(self.config.system_folder, self.config.system_name + '.xdata'), 'w') as f: # make blank xdata file
            for section in self.config.section_order: 
                if section in self._disabled_xdata_sections:
                    continue
                f.write(f"{section}\n\n\n\n")
        
        # this is the main file
        self.config.main_file = os.path.join(self.config.system_folder, self.config.system_name + '.xdata')

        # Clear old sidecars if this output folder is being reused.
        self.config.structure_record_file = os.path.join(self.config.system_folder, self.config.system_name + '.nanotube')
        self.config.salts_record_file = os.path.join(self.config.system_folder, self.config.system_name + '.salts')
        for stale_record in (self.config.structure_record_file, self.config.salts_record_file):
            if os.path.exists(stale_record):
                os.remove(stale_record)


        # copy / serialise input json to output folder
        dst = os.path.join(self.config.system_folder, 'chirality_kit_input.json')
        if isinstance(self.json_input, dict):
            # dict input: save it as JSON
            with open(dst, 'w') as f:
                json.dump(self.json_input, f, indent=4)
        else:
            # path input: copy it as-is
            shutil.copyfile(self.json_input, dst)

        # Save the seed so this run can be replayed.
        self.json_input.setdefault('metadata', {})['rng_seed_used'] = self.config.seed
        with open(dst, 'w') as f:
            json.dump(self.json_input, f, indent=4)

        # The logger starts after apply_seed(), so log the seed again here.
        my_log_file.info(f"RNG seed for this run: {self.config.seed}")
        my_log_file.info(f"Output directory set to {self.config.system_folder}")


    def run(self):

        '''
        Execute the full chirality-kit pipeline: initialise the output folder,
        build the master xdata file from JSON settings, write all requested output
        formats (.pdb, .psf, .xyz, .data), and copy the chosen force-field's
        parameter files (minus the fragments/solvents xdata) into the output
        directory.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        # ============================
        # MASTER FILE BUILDING
        # ============================

        self.write_master_file_from_JSON()

        # ============================
        # OTHER .xxx FORMATS
        # ============================

        for data_type in self.json_input['settings']["output_types"]:
            getattr(self, f"write_{data_type[1:]}_file")()

        # ============================
        # CASE SPECIFICS 
        # ============================

        if '.pdb' in self.json_input['settings']["output_types"]: self.write_pdb_restraints_file()
        if self.config.fixed_box is not None: self.write_xsc_file()

        # ============================
        # SUMMARY FILE
        # ============================

        self._write_summary_file()
        my_log_file.info("Chirality Kit run complete. Thanks for using Chirality Kit! <3")

    def _get_disabled_xdata_sections(self):

        settings = self.json_input["settings"]
        missing = [key for key in XDATA_SECTION_SETTINGS if key not in settings]
        if missing:
            raise ValueError(
                "Missing required settings: " + ", ".join(missing)
            )

        invalid = [
            key for key in XDATA_SECTION_SETTINGS
            if not isinstance(settings[key], bool)
        ]
        if invalid:
            raise ValueError(
                "Xdata section settings must be booleans: " + ", ".join(invalid)
            )

        return {
            section for key, section in XDATA_SECTION_SETTINGS.items()
            if not settings[key]
        }

    def write_master_file_from_JSON(self):

        '''
        Orchestrate the full structure-building sequence driven by the JSON input:
        (0) initialise the force field and settings, (1) generate the nanotube
        structure, (2) attach functional groups, (3) solvate the system, and
        (4) finalise non-bonded interaction coefficients. Logs the final system
        and nanotube charges when complete.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        # -------------------------------------------------------------------------
        # 0. INITIALISATION  STEP
        # -------------------------------------------------------------------------
        
        my_log_file.info("Starting Environment Build")
        self.config.reinitialise(
            force_field       = self.json_input['system']['force-field'],
            nanotube          = self.json_input['system']['nanotube'],
            functionalisation = self.json_input['system']['functionalisation'],
            solvation         = self.json_input['system']['solvation'],
            nanotube_type     = self.json_input.get('nanotube', {}).get('type'),
        )
        self.config.verbose=self.json_input['settings']['verbose']
        
        my_log_file.info(f"Using force field: {self.json_input['system']['force-field']}")
        supported_info = self.config.get_supported_info()
        my_log_file.info(f"Supported Info: {supported_info}")
        if 'low_priority_atom_cutoff' in self.json_input.get('system', {}):
            my_log_file.warning(
                "Ignoring 'low_priority_atom_cutoff': the bulk phase is derived from the "
                "solvation groups now (the highest-priority packed group). The key does "
                "nothing and can be deleted from the input."
            )

        pieces = set() # keeps track of all the pieces we are adding
        fg_types = []
        multi_species_nanotubes = {"BNNT", "MoS2NT", "MoS2", "MoSSeNT", "MoSSe"}

        if self.json_input['system']['size'] == "auto" and self.json_input['system']['nanotube'] == False:  # Box size logic.
            my_log_file.error("Only 'auto' size with 'nanotube' system other wise provide a specific size.")
        self.box_size = self.json_input['system']['size']

        my_log_file.info("Starting Build Validation")
        self._validate_input(
            force_field       = self.json_input['system']['force-field'],
            nanotube          = self.json_input['system']['nanotube'],
            functionalisation = self.json_input['system']['functionalisation'],
            solvation         = self.json_input['system']['solvation'],
        )
        my_log_file.info("Build Validation Successful")

        # A "sheets" block, or system.sheets, fits the box x and y to the graphene
        # lattice. Only system.sheets adds the sheets themselves (step 3b).
        # Every build starts with no sheet state, so a reused IO object is clean.
        self.config.fixed_box = None
        self.config.residue_name_by_index = {}
        self._sheet_lattice = None
        self._sheet_report = None
        self._tube_frame = None
        self._tube_sheet_report = None
        self._around_residues = set()
        self._carved_count = 0
        self._sheet_cfg = self._sheet_settings()

        # -------------------------------------------------------------------------
        # 1. NANOTUBE GENERATION
        # -------------------------------------------------------------------------

        my_log_file.info(f"ADDING NANOSTRUCTURE...   {self.json_input['system']['nanotube']}")
        if self.json_input['system']['nanotube'] == True:
            new_res_idx = "42" # because its the answer to everything

            # Finite tubes fold the top edge back by Lz. Periodic tubes keep one-bond edge atoms.
            periodic = self.json_input['nanotube'].get('periodic', False) is True
            if periodic:
                my_log_file.warning("nanotube.periodic is enabled: the tube is built for "
                                    "periodic simulation along z and may contain atoms with a single bond.")

            if self.json_input['nanotube']['type'] in {"CNT", "CNT-Drude"}:
                
                # -------------------------------------------------------------------------
                # CARBON NANOTUBE GENERATION
                # -------------------------------------------------------------------------

                '''
                Iijima, S. (1991). 
                Helical microtubules of graphitic carbon. Nature, 
                354(6348), 56-58. 
                https://doi.org/10.1038/354056a0
                '''

                self.force_field_manager.write_default_species()
                self.structure.write_nanotube(n=self.json_input['nanotube']['n'], 
                                                        m=self.json_input['nanotube']['m'], 
                                                        Ncells=self.json_input['nanotube']['repeats'],
                                                        PBC=periodic)
                pieces.add("Nanotube-Species")
                n_atoms = self.file_manager.find_number_of_entries_in_section('Atoms')
                self._log_cnt_geometry(self.json_input['nanotube']['n'],
                                    self.json_input['nanotube']['m'],
                                    self.json_input['nanotube']['repeats'],
                                    n_atoms=n_atoms)
            
            elif self.json_input['nanotube']['type'] == "BNNT":
                
                # -------------------------------------------------------------------------
                # BORON NITRIDE NANOTUBE GENERATION
                # -------------------------------------------------------------------------

                '''
                Rubio, A., Corkill, J. L., & Cohen, M. L. (1994). 
                Theory of graphitic boron nitride nanotubes. Physical Review B, 49(7), 5081. 
                https://doi.org/10.1103/physrevb.49.5081
                '''

                self.force_field_manager.write_default_species()
                self.structure._write_bnnt_nanotube(
                    n=self.json_input['nanotube']['n'],
                    m=self.json_input['nanotube']['m'],
                    Ncells=self.json_input['nanotube']['repeats'],
                    PBC=periodic,
                )
                pieces.add("Nanotube-Species")
                n_atoms = self.file_manager.find_number_of_entries_in_section('Atoms')
                self._log_cnt_geometry(self.json_input['nanotube']['n'],
                                    self.json_input['nanotube']['m'],
                                    self.json_input['nanotube']['repeats'],
                                    n_atoms=n_atoms)

            elif self.json_input['nanotube']['type'] in {"MoS2NT", "MoS2"}:

                # -------------------------------------------------------------------------
                # MOLYBDENUM DISULFIDE NANOTUBE GENERATION
                # -------------------------------------------------------------------------

                '''
                Luo, L. et al. (2025).
                Symmetry-broken MoS2 nanotubes through sequential sulfurization
                of MoO2 nanowires. Nature Communications, 16, 8394.
                https://doi.org/10.1038/s41467-025-63333-1
                '''

                self.force_field_manager.write_default_species()
                self.structure._write_mos2_nanotube(
                    n=self.json_input['nanotube']['n'],
                    m=self.json_input['nanotube']['m'],
                    Ncells=self.json_input['nanotube']['repeats'],
                    PBC=periodic,
                )
                pieces.add("Nanotube-Species")
                n_atoms = self.file_manager.find_number_of_entries_in_section('Atoms')
                self._log_cnt_geometry(self.json_input['nanotube']['n'],
                                    self.json_input['nanotube']['m'],
                                    self.json_input['nanotube']['repeats'],
                                    n_atoms=n_atoms)

            elif self.json_input['nanotube']['type'] in {"MoSSeNT", "MoSSe"}:

                # -------------------------------------------------------------------------
                # JANUS MOLYBDENUM SULFUR SELENIDE NANOTUBE GENERATION
                # -------------------------------------------------------------------------

                '''
                Tang, M., Pang, Y., Luo, Y. F., Song, Q., & Wang, M. (2019).
                Electronic properties of Janus MoSSe nanotubes. Computational
                Materials Science, 156, 315-320.
                https://doi.org/10.1016/j.commatsci.2018.10.012
                '''

                self.force_field_manager.write_default_species()
                self.structure._write_mosse_nanotube(
                    n=self.json_input['nanotube']['n'],
                    m=self.json_input['nanotube']['m'],
                    Ncells=self.json_input['nanotube']['repeats'],
                    outer_chalcogen=self.json_input['nanotube'].get('outer-chalcogen', 'Se'),
                    PBC=periodic,
                )
                pieces.add("Nanotube-Species")
                n_atoms = self.file_manager.find_number_of_entries_in_section('Atoms')
                self._log_cnt_geometry(self.json_input['nanotube']['n'],
                                    self.json_input['nanotube']['m'],
                                    self.json_input['nanotube']['repeats'],
                                    n_atoms=n_atoms)

            else:
                raise ValueError(f"Unsupported nanotube type: {self.json_input['nanotube']['type']}")
            

            # -------------------------------------------------------------------------
            # 2. FUNCTIONALISATION
            # -------------------------------------------------------------------------

            my_log_file.info(f"ADDING FUNCTIONAL GROUPS...   {self.json_input['system']['functionalisation']}")
            if self.json_input['system']['functionalisation'] == True:
                if self.json_input['nanotube']['type'] in multi_species_nanotubes:
                    raise NotImplementedError(
                        "Functionalisation is not implemented for multi-species nanotubes yet."
                    )

                structure_info, index_info, type_info = self.force_field_manager.get_atom_positions_residue(residue=new_res_idx)
                bond_info = (
                    self.force_field_manager._create_bonds_data(atom_offset=0)
                    if "Bonds" in self.config.enabled_xdata_sections
                    else self.force_field_manager.generated_bonds_data
                )
                
                self._validate_fg_available_carbons(fg_config=self.json_input['functionalisation'], n_atoms=n_atoms)

                fg_idx, fg_vec, fg_pos, fg_types = self.functionalisation._get_points_of_attachment(
                    structure_info, index_info, bond_info,
                    fg_config=self.json_input['functionalisation'],
                )

                if self.config.verbose:
                    my_log_file.info(f"Attaching functional groups of types:\n {fg_types}")
                    my_log_file.info(f"Attaching functional groups at indices:\n {fg_idx}")
                    my_log_file.info(f"With attachment vectors:\n {fg_vec}")

                self.functionalisation.write_functional_groups(
                    self.config.fragments_file,
                    all_types=fg_types,
                    all_indexes=fg_idx,
                    all_vectors=fg_vec,
                    all_positions=fg_pos,
                    residue_idx=new_res_idx,
                )
                self._log_fg_distribution(fg_types, fg_idx, fg_pos)

                if self.config.enabled_xdata_sections & {
                    "Angles", "Dihedrals", "Impropers"
                }:
                    bond_info = (
                        self.force_field_manager._create_bonds_data(atom_offset=0)
                        if "Bonds" in self.config.enabled_xdata_sections
                        else self.force_field_manager.generated_bonds_data
                    )
                    self.force_field_manager.inject_angles_from_bonds(bond_info, angle_idx_offset=0)
                    self.force_field_manager.inject_dihedrals_from_bonds(bond_info, dihedral_idx_offset=0)
                    self.force_field_manager.inject_impropers_from_bonds(bond_info, improper_idx_offset=0)

                pieces.update(fg_types)

        # Save the bare structure now; ion placement can skip the full solvated file later.
        self.solvation.write_structure_record()

        # Where the tube is, for the solvent sides and the tube-end sheets.
        if self.json_input['system']['nanotube'] == True and (
                self._tube_sheets_on() or self._solvent_sides_used()):
            if self._tube_sheets_on():
                self._tube_frame = self.sheets.tube_frame(
                    rim_offset=float(self._sheet_cfg['rim-offset']),
                    min_separation=2.0 * float(self._sheet_cfg['wall-gap']))
            else:
                self._tube_frame = self.sheets.tube_frame()
            frame = self._tube_frame
            my_log_file.info(
                f"Tube frame: axis ({frame['axis'][0]:.4f}, {frame['axis'][1]:.4f}), "
                f"R = {frame['radius']:.4f} Ang, innermost rim carbons at z = {frame['rim_lo']:.4f} "
                f"and {frame['rim_hi']:.4f} Ang, TUBE residue from z = {frame['tube_z'][0]:.4f} "
                f"to {frame['tube_z'][1]:.4f} Ang."
            )
            if frame['rim_offset']:
                my_log_file.info(
                    f"rim-offset {frame['rim_offset']:+.3f} Ang: tube-end sheet planes at "
                    f"z = {frame['z_lo']:.4f} and {frame['z_hi']:.4f} Ang.")


        # -------------------------------------------------------------------------
        # 3. SOLVATION
        # -------------------------------------------------------------------------

        my_log_file.info(f"ADDING SOLVENTS...   {self.json_input['system']['solvation']}")

        requested_pad_xy = None
        if self.box_size == "auto":
            
            structure_info, index_info, type_info = self.force_field_manager.get_atom_positions_residue(residue=new_res_idx)       
            PBC_pad_xy, PBC_pad_z = self._resolve_auto_box_padding()
            requested_pad_xy = PBC_pad_xy
            structure_info = np.array(structure_info)
            x_max, x_min = np.max(structure_info[:,0]), np.min(structure_info[:,0])
            y_max, y_min = np.max(structure_info[:,1]), np.min(structure_info[:,1])
            z_max, z_min = np.max(structure_info[:,2]), np.min(structure_info[:,2])
            self.box_size = [
                            [x_min - PBC_pad_xy, x_max + PBC_pad_xy],
                            [y_min - PBC_pad_xy, y_max + PBC_pad_xy],
                            [z_min - PBC_pad_z, z_max + PBC_pad_z] ]

        # The solvent must fill an x and y that are whole graphene cells, so the
        # padding is scaled up to the lattice before anything is placed.
        if self._sheet_cfg is not None:
            self._fit_box_to_sheet_lattice(requested_pad_xy)
            self.file_manager.standardise_top_counts()

        if self.json_input['system']['solvation'] == True:

            # Highest priority goes first; ties keep input order.
            solvation_groups = sorted(
                self.json_input['solvation']['groups'],
                key=lambda g: -int(g['priority']),
            )

            # Only the last ionic group neutralises, after all salt charges are present.
            last_ionic_idx = max(
                (i for i, g in enumerate(solvation_groups) if g.get('placement') == 'ionic'),
                default=None,
            )

            for group_idx, sol in enumerate(solvation_groups):

                effective_box = self.box_size
                pad = float(sol.get('pad', 0.0) or 0.0)           # effectively keeping molecules inside the box
                if pad > 0.0 and sol.get('placement') != 'ionic': # ions are points so dont need padding
                    effective_box = [
                        [self.box_size[0][0] + pad, self.box_size[0][1] - pad],
                        [self.box_size[1][0] + pad, self.box_size[1][1] - pad],
                        [self.box_size[2][0] + pad, self.box_size[2][1] - pad],
                    ]
                    my_log_file.info(
                        f"Solvent '{sol['type']}': applying symmetric pad of {pad} Angstrom "
                        f"(effective box "
                        f"{effective_box[0][1]-effective_box[0][0]:.2f} x "
                        f"{effective_box[1][1]-effective_box[1][0]:.2f} x "
                        f"{effective_box[2][1]-effective_box[2][0]:.2f} Angstrom)."
                    )

                # side: full keeps the whole box; the other sides cut it in z.
                side = sol.get('side', 'full')
                passes = self._side_passes(sol, effective_box)
                if side != 'full':
                    for _, part_box in passes:
                        my_log_file.info(
                            f"Solvent '{sol['type']}': side {side}, z from "
                            f"{part_box[2][0]:.3f} to {part_box[2][1]:.3f} Angstrom."
                        )
                residues_before = (int(self.file_manager.find_last_residue())
                                   if side == 'around-tube' else None)

                n_atom_prior = self.file_manager.find_number_of_entries_in_section('Atoms')
                for part, part_box in passes:
                    self.solvation.write_solvent(
                        part,
                        box_size = part_box,
                        bridson_params= self.json_input['solvation'].get('bridson_params', None),
                        is_last_ionic= (group_idx == last_ionic_idx),
                        no_go_slab= self._ion_slab(sol),
                    )
                if residues_before is not None:
                    self._drop_lumen_residues(residues_before)
                    last_residue = int(self.file_manager.find_last_residue())
                    self._around_residues.update(range(residues_before + 1, last_residue + 1))
                n_atoms_after = self.file_manager.find_number_of_entries_in_section('Atoms')

                if sol["placement"] == "ionic":
                    cathode_type, cathode_count, anode_type, anode_count = self.solvation._process_salt_formula(sol['type'])
                    my_log_file.info(f"Added {cathode_count} {cathode_type} and {anode_count} {anode_type} ions for solvation group {sol['type']}.")
                    pieces.add(cathode_type)
                    pieces.add(anode_type)
                else:
                    new_solvs = sum(1 for atom in self.force_field_manager.get_atom_positions_all()[2] if sol['type'] in atom)
                    self._log_solvation_stats(sol, n_atom_prior, n_atoms_after, new_solvs)
                    pieces.add(sol['type'])


        # -------------------------------------------------------------------------
        # 3a. TUBE-END SHEETS
        # -------------------------------------------------------------------------

        # The solvent between the sheet planes outside the tube, and any within
        # wall-gap of a sheet carbon, goes first. The sheets use no random numbers.
        if self._tube_sheets_on():
            my_log_file.info("ADDING TUBE-END SHEETS...")
            tube_sheets = self.sheets.tube_sheet_sites(self.box_size, self._sheet_cfg, self._tube_frame)
            self._carve_for_tube_sheets(tube_sheets)
            self._tube_sheet_report = self.sheets.add_tube_sheets(
                self._sheet_cfg, tube_sheets, self._tube_frame)
            self._tube_sheet_report["carved"] = self._carved_count
            self.file_manager.standardise_top_counts()


        # -------------------------------------------------------------------------
        # 3b. SHEETS
        # -------------------------------------------------------------------------

        # Added last and without random numbers, so a sheets-off build with the
        # same "sheets" block and seed matches this one in every other atom.
        if self.json_input['system'].get('sheets', False) is True:
            my_log_file.info("ADDING SHEETS...")
            self._sheet_report = self.sheets.add_sheets(self.box_size, self._sheet_cfg)
            self.box_size = self._sheet_report["box"]
            self.config.fixed_box = [list(axis) for axis in self.box_size]
            self.file_manager.standardise_top_counts()


        # -------------------------------------------------------------------------
        # 4. FINALISE NON-BONDED INTERACTIONS
        # -------------------------------------------------------------------------

        my_log_file.info("FINALISING NON-BONDED INTERACTIONS...")
        # Sort pieces so NBFIX output doesn't depend on Python's hash order.
        for piece in sorted(pieces):
            mol = self.file_manager.find_and_extract_subsection(self.config.fragments_file, "Fragment", piece)
            if not mol:
                mol = self.file_manager.find_and_extract_subsection(self.config.solvents_file, "Solvent", piece)
            if not mol:
                my_log_file.error(f"Could not find piece {piece} in fragments or solvents file.")
                raise ValueError(f"Could not find piece {piece} in fragments or solvents file.")
            self.force_field_manager.inject_lj_coefficients(mol)
        my_log_file.info("Environment Build Complete")

        # -------------------------------------------------------------------------
        # 5. COPY OVER FORCE FIELD PARAMETER FILES
        # -------------------------------------------------------------------------

        self.copy_other_FF_files()

        # -------------------------------------------------------------------------
        # 6. FINAL LOG MESSAGES
        # -------------------------------------------------------------------------

        self._print_charge_summary()
        self._check_pbc_box()

    def create_system_name(self): 

        '''
        Construct a descriptive system name string from the JSON settings,
        combining the force-field name, nanotube chirality, functional group
        types and counts, and solvent placement descriptors.

        Parameters:
        ----------
        None

        Returns:
        -------
        name : str
            Human-readable system identifier used for naming the output folder
            and all output files.

        Raises:
        ------
        None
        '''

        def sanitize(s):
            return (
                str(s)
                .replace("+", "")
                .replace("-", "")
                .replace(" ", "")
                .replace(".", "p")
            )

        j = self.json_input
        parts = []

        # --- Force field ---
        parts.append(sanitize(j["system"]["force-field"]))

        # --- Nanotube ---
        if j["system"].get("nanotube"):
            nt = j["nanotube"]
            parts.append(f"{sanitize(nt['type'])}_({nt['n']},{nt['m']})x{nt['repeats']}")

        # --- Functionalisation ---
        if j["system"].get("functionalisation") and j["system"].get("nanotube"):
            func = j["functionalisation"]

            summary = {}

            # --- Term groups ---
            term = func.get("term", {})
            for fg in term.get("groups", []):
                if fg["type"] in VACANT_FG_TYPES:
                    continue
                count = fg.get("count", 0)
                if count > 0:
                    summary[fg["type"]] = summary.get(fg["type"], 0) + count

            # --- Ring groups (count x ring-count) ---
            for ring in func.get("rings", []):
                ring_multiplier = ring.get("ring-count", 0)

                for fg in ring.get("groups", []):
                    if fg["type"] in VACANT_FG_TYPES:
                        continue
                    count = fg.get("count", 0)
                    if count > 0 and ring_multiplier > 0:
                        total = count * ring_multiplier
                        summary[fg["type"]] = summary.get(fg["type"], 0) + total

            # --- Loose groups ('all-*' placements have no fixed count) ---
            for fg in func.get("loose", {}).get("groups", []):
                fg_type = fg["type"]
                if str(fg.get("loose-placement", "")).startswith("all"):
                    summary[fg_type] = "all"
                    continue
                count = fg.get("count", 0)
                if count > 0 and summary.get(fg_type) != "all":
                    summary[fg_type] = summary.get(fg_type, 0) + count

            # --- Build string with x notation ---
            fg_parts = []
            for k, v in sorted(summary.items()):
                fg_parts.append(f"{sanitize(k)}x{v}")

            if fg_parts:
                parts.append("_".join(fg_parts))

            # --- Terminal hydrogenation flag ---
            if term.get("term-hydrogenate"):
                parts.append("Hterm")

        # --- Solvation ---
        if j["system"].get("solvation"):
            sol_parts = []

            for sol in j["solvation"]["groups"]:
                t = sanitize(sol.get("type"))
                count = sol.get("count", None)

                if sol.get("placement") == "ionic" and sol.get("charge-neutrality"):
                    sol_parts.append(f"{t}_neutral")
                elif count is not None and count > 0:
                    sol_parts.append(f"{t}{count}")
                else:
                    sol_parts.append(t)

            if sol_parts:
                # join separately so trailing "_" behavior is clean
                parts.append("_".join(sol_parts))

        # --- Sheets (a "sheets" block alone gives the lattice-matched box) ---
        tube_sheets = j["system"].get("tube-sheets") is True
        if j["system"].get("sheets") is True:
            parts.append("SHEETS")
        elif "sheets" in j and not tube_sheets:
            parts.append("LATTICE")
        if tube_sheets:
            parts.append("TUBESHEETS")

        return "_".join(parts)
    
    
    # -------------------------------------------------------------------------
    # REPORTING HELPERS
    # -------------------------------------------------------------------------

    def _nanotube_geometry_vectors(self, nt_type):

        if nt_type == "BNNT":
            bond_length = 1.4460
            a1 = np.array([np.sqrt(3)*bond_length, 0.0])
            a2 = np.array([np.sqrt(3)/2*bond_length, 1.5*bond_length])
            return a1, a2, bond_length

        if nt_type in {"MoS2NT", "MoS2"}:
            bond_length = 2.4100
            lattice_constant = 3.1600
            a1 = np.array([lattice_constant, 0.0])
            a2 = np.array([0.5*lattice_constant, np.sqrt(3)/2*lattice_constant])
            return a1, a2, bond_length

        if nt_type in {"MoSSeNT", "MoSSe"}:
            bond_length = 2.4100
            lattice_constant = 3.2300
            a1 = np.array([lattice_constant, 0.0])
            a2 = np.array([0.5*lattice_constant, np.sqrt(3)/2*lattice_constant])
            return a1, a2, bond_length

        bond_length = self.config.get_default_species_bond_param("bond_length", 1.44)
        a1 = np.array([np.sqrt(3)*bond_length, 0.0])
        a2 = np.array([np.sqrt(3)/2*bond_length, 1.5*bond_length])
        return a1, a2, bond_length

    def _nanochannel_summary_profile(self, nt_type):

        default_cutoff = self.config.get_default_species_bond_param("bond_cutoff", 1.60)
        default_profile = {
            "summary_elements": None,
            "expected_coordination": 3,
            "bond_cutoff": default_cutoff,
            "formula_factor": 4,
            "xyz_prefix": None,
            "xyz_outer_chalcogen": False,
        }

        profiles = {
            "BNNT": {
                "bond_cutoff": 1.65,
                "xyz_prefix": "SW-BNNT",
            },
            "MoS2NT": {
                "summary_elements": {"Mo"},
                "expected_coordination": 6,
                "bond_cutoff": 2.95,
                "formula_factor": 2,
                "xyz_prefix": "SW-MoS2NT",
            },
            "MoS2": {
                "summary_elements": {"Mo"},
                "expected_coordination": 6,
                "bond_cutoff": 2.95,
                "formula_factor": 2,
                "xyz_prefix": "SW-MoS2NT",
            },
            "MoSSeNT": {
                "summary_elements": {"Mo"},
                "expected_coordination": 6,
                "bond_cutoff": 3.05,
                "formula_factor": 2,
                "xyz_prefix": "SW-MoSSeNT",
                "xyz_outer_chalcogen": True,
            },
            "MoSSe": {
                "summary_elements": {"Mo"},
                "expected_coordination": 6,
                "bond_cutoff": 3.05,
                "formula_factor": 2,
                "xyz_prefix": "SW-MoSSeNT",
                "xyz_outer_chalcogen": True,
            },
        }

        profile = default_profile.copy()
        profile.update(profiles.get(nt_type, {}))
        return profile

    def _nanochannel_xyz_candidates(self, nt, profile):

        n, m, Ncells = nt['n'], nt['m'], nt['repeats']
        prefix = profile.get("xyz_prefix")
        if prefix is None:
            primary = self.config.default_species_primary or {}
            element = primary.get("atom_element", "C")
            prefix = f"SW-{element}NT"

        suffix = ""
        if profile.get("xyz_outer_chalcogen"):
            suffix = f"_outer-{nt.get('outer-chalcogen', 'Se')}"
        basenames = [f"{prefix}_{n}_{m}_{Ncells}{suffix}.xyz"]

        folders = [
            getattr(self.config, "system_folder", None),
            os.path.dirname(getattr(self.config, "main_file", "") or ""),
            getattr(self.config, "called_from", None),
        ]
        folders = [folder for folder in folders if folder]

        return [
            os.path.join(folder, basename)
            for folder in folders
            for basename in basenames
        ]

    def _read_xyz_records(self, filename):

        with open(filename, "r") as f:
            lines = f.readlines()

        records = []
        for line in lines[2:]:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                records.append({
                    "element": parts[0],
                    "x": float(parts[1]),
                    "y": float(parts[2]),
                    "z": float(parts[3]),
                })
            except ValueError:
                continue
        return records

    def _nanochannel_summary_indices(self, records, profile):

        summary_elements = profile.get("summary_elements")
        if summary_elements is None:
            return list(range(len(records)))
        return [
            i for i, record in enumerate(records)
            if record["element"] in summary_elements
        ]

    def _nanochannel_count_pair(self, record_i, record_j, profile):

        summary_elements = profile.get("summary_elements")
        if summary_elements is None:
            return True
        return not (
            record_i["element"] in summary_elements
            and record_j["element"] in summary_elements
        )

    def _nanochannel_counts_from_xyz(self, records, profile):

        summary_indices = self._nanochannel_summary_indices(records, profile)
        if not summary_indices:
            return None

        positions = np.array(
            [[record["x"], record["y"], record["z"]] for record in records],
            dtype=float,
        )
        bond_cutoff = float(profile["bond_cutoff"])
        expected_coordination = int(profile["expected_coordination"])

        degrees = {}
        for i in summary_indices:
            degree = 0
            for j, record_j in enumerate(records):
                if i == j:
                    continue
                if not self._nanochannel_count_pair(records[i], record_j, profile):
                    continue
                if np.linalg.norm(positions[i] - positions[j]) <= bond_cutoff:
                    degree += 1
            degrees[i] = degree

        rim_indices = [
            i for i, degree in degrees.items()
            if degree < expected_coordination
        ]
        ring_atoms = None
        if rim_indices:
            midpoint = float(np.mean(positions[summary_indices, 2]))
            entry = sum(1 for i in rim_indices if positions[i, 2] < midpoint)
            exit_ = len(rim_indices) - entry
            ring_atoms = (entry, exit_)

        return {
            "channel_atoms": len(summary_indices),
            "ring_atoms": ring_atoms,
        }

    def _nanochannel_formula_atom_count(self, nt, profile):

        n, m, Ncells = nt['n'], nt['m'], nt['repeats']
        dR = self.geometry.gcd(2*n + m, 2*m + n)
        unit_area = n*n + n*m + m*m
        return int(round(profile["formula_factor"] * unit_area / dR * Ncells))

    def _nanochannel_summary_counts(self, nt):

        profile = self._nanochannel_summary_profile(nt['type'])
        counts = None

        for filename in self._nanochannel_xyz_candidates(nt, profile):
            if not os.path.isfile(filename):
                continue
            try:
                counts = self._nanochannel_counts_from_xyz(
                    self._read_xyz_records(filename),
                    profile,
                )
            except OSError:
                counts = None
            if counts is not None:
                break

        if counts is None:
            counts = {"channel_atoms": None, "ring_atoms": None}

        if counts["channel_atoms"] is None:
            counts["channel_atoms"] = self._nanochannel_formula_atom_count(nt, profile)
        if counts["ring_atoms"] is None:
            fallback_ring_atoms = int(nt['n']) + int(nt['m'])
            counts["ring_atoms"] = (fallback_ring_atoms, fallback_ring_atoms)

        return counts

    def _format_nanochannel_ring_atoms(self, ring_atoms):

        entry, exit_ = ring_atoms
        if entry == exit_:
            return str(entry)
        return f"[{entry}, {exit_}]"

    def _log_cnt_geometry(self, n, m, Ncells, n_atoms=None):

        nt_type  = self.json_input['nanotube']['type']
        a1, a2, bond_length = self._nanotube_geometry_vectors(nt_type)

        Ch_len = np.linalg.norm(n*a1 + m*a2)
        R      = Ch_len / (2 * np.pi)

        dR     = self.geometry.gcd(2*n + m, 2*m + n)
        t1, t2 = (2*m + n) // dR, -(2*n + m) // dR
        T_len  = np.linalg.norm(t1*a1 + t2*a2)
        length = Ncells * T_len

        chiral_angle = np.degrees(np.arctan2(np.sqrt(3)*m, 2*n + m))

        if   m == 0: ctype = "zigzag"
        elif n == m: ctype = "armchair"
        else:        ctype = "chiral"

        element = self.config.default_species_primary["atom_element"]

        msg = (
            f"\n{'='*60}\n"
            f"  {nt_type} GEOMETRY  {element}({n},{m}) x {Ncells}\n"
            f"{'='*60}\n"
            f"  Type          : {ctype}\n"
            f"  Radius        : {R:.3f} Ang\n"
            f"  Length        : {length:.3f} Ang\n"
            f"  Chiral angle  : {chiral_angle:.2f} deg\n"
            f"  Bond length   : {bond_length:.3f} Ang\n"
            f"  Atoms total   : {n_atoms}\n"
            f"{'='*60}"
        )
        my_log_file.info(msg)
        print(msg)

    def _log_fg_distribution(self, fg_types, fg_idx, fg_pos):

        '''
        Log the distribution of functional groups: total count per type and
        how many landed on the entry end (z < 0) vs the exit end (z >= 0).

        Parameters:
        ----------
        fg_types : array-like of str
            Functional group type for each attachment site.
        fg_idx : np.ndarray
            LAMMPS atom indices of the attachment carbons.
        fg_pos : (M,3) np.ndarray
            Cartesian positions of the attachment sites.

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        from collections import Counter
        type_counts = Counter(fg_types)
        fg_pos = np.array(fg_pos)

        out = [f"\n{'='*60}", "  FUNCTIONAL GROUP DISTRIBUTION", f"{'='*60}"]
        for fg_type, total in sorted(type_counts.items()):
            mask    = np.array([t == fg_type for t in fg_types])
            n_entry = int(np.sum(fg_pos[mask, 2] <  0))
            n_exit  = int(np.sum(fg_pos[mask, 2] >= 0))
            out.append(f"  {fg_type:<12s}: {total:3d} total  "
                       f"(entry end: {n_entry}, exit end: {n_exit})")
        out.append(f"{'='*60}")
        msg = "\n".join(out)
        my_log_file.info(msg)
        print(msg)

    def _log_solvation_stats(self, sol, n_atoms_prior, n_atoms_after, n_atoms_added):

        '''
        Log placement results for a solvation group.

        n_atoms_added is the post-cull count of atoms belonging to this solvent
        type (computed at the call site by scanning the Atoms section). The
        system-wide delta (n_atoms_after - n_atoms_prior) is reported separately
        because the cull may also remove atoms from earlier groups or the
        structure itself, so a negative or smaller-than-expected net delta does
        NOT necessarily mean this solvent was culled.

        Parameters:
        ----------
        sol : dict
            Solvation group dict from the JSON input.
        n_atoms_prior : int
            System-wide atom count before this group was placed.
        n_atoms_after : int
            System-wide atom count after placement AND overlap removal.
        n_atoms_added : int
            Post-cull count of atoms belonging to this solvent type. Not the
            number attempted; that figure is internal to write_solvent.

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        net_system_delta = n_atoms_after - n_atoms_prior

        # Need this for both count and concentration modes.
        try:
            solvent_mol_lines  = self.file_manager.find_and_extract_subsection(
                self.config.solvents_file, 'Solvent', sol['type'])
            atoms_section      = self.file_manager.find_and_extract_lines(
                solvent_mol_lines, 'Atoms', header=False)
            atoms_per_molecule = len(atoms_section) if atoms_section else 1
        except Exception:
            atoms_per_molecule = 1

        mpu               = sol.get('molecules_per_unit', 1)
        n_kept_molecules  = n_atoms_added / atoms_per_molecule / mpu

        out = [
            f"  Solvent '{sol['type']}': "
            f"{n_atoms_added} atoms kept "
            f"({int(round(n_kept_molecules))} molecules, {atoms_per_molecule} atoms/mol).",
            f"=" * 150 + "\n"
            f"  System net atom delta: {net_system_delta:+d} "
            f"(includes any atoms culled from earlier groups or structure)."
        ]

        #  target: concentration (mol/L) 
        if sol.get('concentration') is not None and self.box_size not in (None, "auto"):
            vol_ang3 = self.solvation._box_volume_ang3(self.box_size)
            vol_L    = self.solvation._ang3_to_liters(vol_ang3)
            NA       = 6.02214076e23  # Avocado number
            achieved = n_kept_molecules / (NA * vol_L)
            out.append(
                f"  Target: {sol['concentration']:.3f} mol/L  |  "
                f"Achieved: {achieved:.3f} mol/L"
                f"\n" + f"=" * 150
            )

        #  target: count 
        if sol.get('count') is not None:
            target   = int(sol['count'])
            achieved = int(round(n_kept_molecules))
            delta    = achieved - target
            sign     = "+" if delta >= 0 else ""
            pct      = (delta / target * 100.0) if target else 0.0
            out.append(
                f"  Target: {target} molecules  |  "
                f"Achieved: {achieved} ({sign}{delta}, {sign}{pct:.1f}%)"
                f"\n" + f"=" * 150
            )

        msg = "\n".join(out)
        my_log_file.info(msg)
        print(msg)

    def _print_charge_summary(self):

        '''
        Print and log a formatted charge summary showing tube charge, solution
        charge, and total system charge. Warns if the total charge exceeds 0.01 e.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        total = self.force_field_manager._round_charge(
            self.force_field_manager.get_system_charge())

        out = [f"\n{'='*60}", "  CHARGE SUMMARY", f"{'='*60}"]

        if self.json_input['system']['nanotube']:
            tube = self.force_field_manager._round_charge(
                self.force_field_manager.get_system_charge(residue_index="42"))
            sol  = self.force_field_manager._round_charge(total - tube)
            out.append(f"  Nanotube charge  : {tube:+.6f} e")
            out.append(f"  Solution charge  : {sol:+.6f} e")

        out.append(f"  Total charge     : {total:+.6f} e")

        CHARGE_TOL = 0.01
        if abs(total) > CHARGE_TOL:
            out.append(f"  WARNING: System is NOT charge-neutral!")
            out.append(f"{'='*60}")
            msg = "\n".join(out)
            my_log_file.warning(msg)
        else:
            out.append(f"  System is charge-neutral.")
            out.append(f"{'='*60}")
            msg = "\n".join(out)
            my_log_file.info(msg)
        print(msg)

    def _write_summary_file(self):

        '''
        Write a plain-text .summary file to the output folder covering CNT
        geometry, functionalisation, solvation, box dimensions, and final charges.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        path = os.path.join(self.config.system_folder,
                            self.config.system_name + '.summary')
        out = []
        SEP = "=" * 60

        out += [SEP, "  CHIRALITY-KIT RUN SUMMARY",
                f"  {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
                SEP, ""]

        out.append(f"Force field : {self.json_input['system']['force-field']}")
        out.append(f"System name : {self.config.system_name}")
        if self.json_input['system'].get('pbc-override', False):
            out.append("PBC override: enabled (requested reduced padding retained)")
        out.append("")

        if self.json_input['system']['nanotube']:
            nt = self.json_input['nanotube']
            n, m, Ncells = nt['n'], nt['m'], nt['repeats']
            a1, a2, _ = self._nanotube_geometry_vectors(nt['type'])
            R      = np.linalg.norm(n*a1 + m*a2) / (2*np.pi)
            dR     = self.geometry.gcd(2*n + m, 2*m + n)
            t1, t2 = (2*m + n)//dR, -(2*n + m)//dR
            length = np.linalg.norm(t1*a1 + t2*a2) * Ncells
            chiral_angle = np.degrees(np.arctan2(np.sqrt(3)*m, 2*n + m))
            if   m == 0: ctype = "zigzag"
            elif n == m: ctype = "armchair"
            else:        ctype = "chiral"
            nanochannel_counts = self._nanochannel_summary_counts(nt)

            out += [f"[Nanotube]",
                    f"  Type          : {nt['type']} ({ctype})",
                    f"  Indices       : ({n},{m})",
                    f"  Repeats       : {Ncells}",
                    f"  Radius        : {R:.4f} Ang",
                    f"  Length        : {length:.4f} Ang",
                    f"  Chiral angle  : {chiral_angle:.3f} deg",
                    f"  Channel atoms : {nanochannel_counts['channel_atoms']}",
                    f"  Ring atoms    : {self._format_nanochannel_ring_atoms(nanochannel_counts['ring_atoms'])}",
                    ""]

        if self.json_input['system']['functionalisation'] and self.json_input['system']['nanotube']:
            fg_cfg = self.json_input['functionalisation']
            out.append("[Functionalisation]")

            # term
            term_cfg = fg_cfg.get('term', {}) or {}
            out.append(f"  term-hydrogenate : {term_cfg.get('term-hydrogenate', False)}")
            out.append(f"  term-start-highest-x : {term_cfg.get('term-start-highest-x', False)}")
            term_groups = [
                g for g in term_cfg.get('groups', [])
                if g.get('count', 0) > 0 and g.get('type') not in VACANT_FG_TYPES
            ]
            if term_groups:
                out.append("  term groups:")
                for fg in term_groups:
                    tag = f" (side {fg.get('side', 'both')})"
                    out.append(f"    {fg['type']:<10s}: {fg['count']}{tag}")

            # rings
            rings_cfg = fg_cfg.get('rings', []) or []
            for r_idx, ring in enumerate(rings_cfg):
                groups_here = [
                    g for g in ring.get('groups', [])
                    if g.get('count', 0) > 0 and g.get('type') not in VACANT_FG_TYPES
                ]
                if not groups_here:
                    continue
                out.append(f"  ring block {r_idx} ({ring.get('ring-placement', '?')}):")
                out.append(f"    ring-count       : {ring.get('ring-count', 0)}")
                out.append(f"    phase-increment  : {ring.get('ring-phase-increment', 0.0):.5f} rad")
                out.append(f"    phase-start      : {ring.get('ring-phase-start-offset', 0.0):.5f} rad")
                out.append(f"    padding          : {ring.get('ring-padding', True)}")
                for fg in groups_here:
                    out.append(f"    {fg['type']:<10s}: {fg['count']}")

            # loose
            loose_groups = (fg_cfg.get('loose', {}) or {}).get('groups', []) or []
            loose_groups = [g for g in loose_groups
                            if g.get('count', 0) > 0
                            or str(g.get('loose-placement', '')).startswith('all')]
            if loose_groups:
                out.append("  loose groups:")
                for fg in loose_groups:
                    mode = fg.get('loose-placement', 'random')
                    if str(mode).startswith('all'):
                        z_from = float(fg.get('z-from', 0.0))
                        z_to   = float(fg.get('z-to', 1.0))
                        out.append(f"    {fg['type']:<10s}: all free atoms, "
                                   f"z {z_from:.2f}-{z_to:.2f} ({mode})")
                    else:
                        out.append(f"    {fg['type']:<10s}: {fg['count']} ({mode})")

            out.append("")

        if self.json_input['system']['solvation']:
            out.append("[Solvation]")
            for sol in self.json_input['solvation']['groups']:
                out.append(f"  {sol['type']:<14s}: {sol['placement']}")
                if sol.get('concentration'):
                    out.append(f"    target conc : {sol['concentration']} mol/L")
                if sol.get('count'):
                    out.append(f"    count       : {sol['count']}")
            out.append("")

        if self.box_size not in (None, "auto"):
            b = self.box_size
            out += ["[Box size]",
                    f"  x : {b[0][0]:.3f} to {b[0][1]:.3f} Ang  ({b[0][1]-b[0][0]:.3f} Ang)",
                    f"  y : {b[1][0]:.3f} to {b[1][1]:.3f} Ang  ({b[1][1]-b[1][0]:.3f} Ang)",
                    f"  z : {b[2][0]:.3f} to {b[2][1]:.3f} Ang  ({b[2][1]-b[2][0]:.3f} Ang)", ""]

        if self._sheet_lattice is not None:
            lat = self._sheet_lattice
            out.append("[Sheets]")
            out.append(f"  graphene cells : {lat['nx']} x {lat['ny']} "
                       f"({self._sheet_cfg['orientation']}, b = {lat['bond_length']:.4f} Ang)")
            out.append(f"  cell periods   : {lat['period_x']:.4f} x {lat['period_y']:.4f} Ang")
            out.append(f"  box x, y       : {lat['lx']:.4f} x {lat['ly']:.4f} Ang")
            if lat['requested_pad_xy'] is not None:
                out.append(f"  xy padding     : {lat['requested_pad_xy']:.3f} Ang requested -> "
                           f"{lat['requested_pad_xy'] + lat['raised_x']:.3f} (x), "
                           f"{lat['requested_pad_xy'] + lat['raised_y']:.3f} (y) Ang")
            else:
                out.append(f"  xy scaled up   : +{2 * lat['raised_x']:.3f} (x), "
                           f"+{2 * lat['raised_y']:.3f} (y) Ang")
            report = self._sheet_report
            if report is None and self._tube_sheet_report is not None:
                out.append("  box-end sheets : off")
            elif report is None:
                out.append("  sheets         : off (lattice-matched box only)")
            else:
                for side in SHEET_SIDES:
                    drudes = report['drudes_per_sheet'][side]
                    out.append(f"  {report['resnames'][side]:<4s} ({side})   : {report['atoms_per_sheet'][side]} "
                               f"{report['atom_type']} atoms"
                               + (f" + {drudes} Drude particles" if drudes else "")
                               + f" at z = {report['z'][side]:.4f} Ang")
                out.append(f"  wall-gap       : {float(self._sheet_cfg['wall-gap']):.3f} Ang")
                out.append(f"  seam-gap       : {float(self._sheet_cfg['seam-gap']):.3f} Ang")
                out.append(f"  sheet charge   : {report['charge']:+.6f} e")
            tube = self._tube_sheet_report
            if tube is not None:
                for side in TUBE_SHEET_SIDES:
                    drudes = tube['drudes_per_sheet'][side]
                    out.append(f"  {tube['resnames'][side]:<4s} ({side}) : {tube['atoms_per_sheet'][side]} "
                               f"{tube['atom_type']} atoms"
                               + (f" + {drudes} Drude particles" if drudes else "")
                               + f" at z = {tube['z'][side]:.4f} Ang, {tube['removed'][side]} sites "
                               f"removed for the hole, {tube['edge'][side]} edge carbons")
                out.append(f"  hole radius    : {tube['hole_radius']:.4f} Ang (tube R + hole-buffer "
                           f"{float(self._sheet_cfg['hole-buffer']):.3f} Ang)")
                out.append(f"  rim-offset     : {self._tube_frame['rim_offset']:+.3f} Ang from the innermost "
                           f"rim carbons (z = {self._tube_frame['rim_lo']:.4f} and "
                           f"{self._tube_frame['rim_hi']:.4f} Ang)")
                out.append(f"  closest to tube: {tube['min_tube_distance']:.4f} Ang")
                out.append(f"  carved         : {tube['carved']} solvent residues (listed in .carved)")
                out.append(f"  tube-end charge: {tube['charge']:+.6f} e")
            out.append("")

        n_atoms = self.file_manager.find_number_of_entries_in_section("Atoms")
        total   = self.force_field_manager._round_charge(
                      self.force_field_manager.get_system_charge())
        out += ["[Final system]", f"  Total atoms   : {n_atoms}",
                f"  Total charge  : {total:+.6f} e"]

        psf_report = getattr(self, "_psf_charge_neutrality_report", None)
        if psf_report:
            if psf_report.get("applied"):
                out.append(
                    "  PSF charge fix: "
                    f"atom {psf_report['psf_atom_index']} "
                    f"{psf_report['atom_type']} "
                    f"residue {psf_report['residue_index']}; "
                    f"q {psf_report['original_charge']:+.6f} -> "
                    f"{psf_report['adjusted_charge']:+.6f}; "
                    f"PSF total {psf_report['original_total']:+.6f} -> "
                    f"{psf_report['adjusted_total']:+.6f} e"
                )
                out.append("  PSF charge fix snapshot:")
                out.extend(f"    {line}" for line in psf_report.get("snapshot", []))
            else:
                out.append(
                    "  PSF charge fix: skipped; "
                    f"PSF total {psf_report['original_total']:+.6f} e; "
                    f"{psf_report.get('reason', 'no adjustment made')}"
                )

        if self.json_input['system']['nanotube']:
            tube = self.force_field_manager._round_charge(
                self.force_field_manager.get_system_charge(residue_index="42"))
            out.append(f"  Tube charge   : {tube:+.6f} e")

        out += ["", SEP]

        with open(path, 'w') as f:
            f.write("\n".join(out) + "\n")

        my_log_file.info(f"Summary written to {path}")
        print(f"\nSummary written to: {path}")

    def _convert_to_base36(self, n):

        '''
        Convert a non-negative integer to a base-36 string using digits 0-9
        followed by uppercase letters A-Z. Used to generate compact, unique
        atom-name suffixes for nanotube carbons in PDB output (e.g. C0, CA, CZ).

        Parameters:
        ----------
        n : int
            Non-negative integer to convert.

        Returns:
        -------
        label : str
            Base-36 string representation.

        Raises:
        ------
        None
        '''

        symbols = list(string.digits + string.ascii_uppercase)
        base = len(symbols)

        label = ''
        if n == 0:
            return '0'
        while n:             
            label = symbols[n % base] + label
            n = n // base
        return label
    
    def _check_pbc_box(self):

        '''
        Check that the z-padding between the TUBE residue and the box boundary
        is sufficient to avoid periodic image interactions (minimum 10 Ang).
        The TUBE residue is the tube with its groups and caps, as the PDB and
        PSF name it; sheets, solvent and ions never count. Logs a warning if
        either side is too short.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        if self.box_size in (None, "auto") or not self.json_input['system']['nanotube']:
            return

        min_pad = 10.0

        atom_pos, _, atom_types, atom_residues = self.force_field_manager.get_atom_positions_all()

        # Same naming rules as the PDB and PSF writers.
        def _in_tube_residue(atom_type, residue_index):
            if residue_index in self.config.residue_name_by_index:
                return False
            parts = atom_type.split('-')
            if len(parts) > 1:
                return parts[1] == "DRUD" and parts[0] not in self.config.dions
            return not (atom_type in self.config.ions
                        or atom_type in self.config.dions
                        or atom_type in self.config.dion_particles)

        tube_mask = np.array([_in_tube_residue(t, r) for t, r in zip(atom_types, atom_residues)])
        if not tube_mask.any():
            return

        z_tube_min = float(np.min(atom_pos[tube_mask, 2]))
        z_tube_max = float(np.max(atom_pos[tube_mask, 2]))
        z_box_lo   = self.box_size[2][0]
        z_box_hi   = self.box_size[2][1]
        pad_lo     = z_tube_min - z_box_lo
        pad_hi     = z_box_hi  - z_tube_max

        out = [f"\n{'='*60}", "  PBC BOX CHECK", f"{'='*60}",
               f"  Tube z-extent : {z_tube_min:.2f} to {z_tube_max:.2f} Ang",
               f"  Box z-extent  : {z_box_lo:.2f} to {z_box_hi:.2f} Ang",
               f"  Padding lo    : {pad_lo:.2f} Ang",
               f"  Padding hi    : {pad_hi:.2f} Ang"]

        # With tube-end sheets each reservoir runs from its sheet to the box end
        # sheet, or round the periodic boundary to the other side without one.
        tube = self._tube_sheet_report
        if tube is not None:
            ends = self._sheet_report['z'] if self._sheet_report is not None else None
            low_end = ends['z-min'] if ends else z_box_lo
            high_end = ends['z-max'] if ends else z_box_hi
            out.append(f"  Reservoir lo  : {low_end:.2f} to {tube['z']['tube-z-min']:.2f} Ang")
            out.append(f"  Reservoir hi  : {tube['z']['tube-z-max']:.2f} to {high_end:.2f} Ang")
            if ends is None:
                out.append("  Reservoirs meet across the periodic z boundary.")

        warn = (pad_lo < min_pad or pad_hi < min_pad)
        if warn:
            out.append(f"  WARNING: z-padding < {min_pad} Ang: "
                       f"periodic image interactions likely!")
        else:
            out.append("  Box z-padding OK.")
        out.append(f"{'='*60}")
        msg = "\n".join(out)
        my_log_file.warning(msg) if warn else my_log_file.info(msg)
        print(msg)

    def _resolve_auto_box_padding(self):

        '''
        Resolve automatic box padding, enforcing the recommended padding unless
        pbc-override is explicitly enabled. Override keeps user-provided small
        boxes but still logs risk warnings.

        Returns:
        -------
        tuple[float, float]
            The xy and z padding values to use for automatic box construction.
        '''

        recommended_pad = 14.0
        system_cfg = self.json_input['system']
        PBC_pad_xy = float(system_cfg.get('size-padding-xy', recommended_pad))
        PBC_pad_z = float(system_cfg.get('size-padding-z', recommended_pad))
        pbc_override = system_cfg.get('pbc-override', False) is True

        if pbc_override:
            if PBC_pad_xy < recommended_pad:
                my_log_file.warning(
                    f"XY Padding for box size is below the recommended "
                    f"{recommended_pad} Angstroms because pbc-override is enabled; "
                    "not increasing."
                )
            if PBC_pad_z < recommended_pad:
                my_log_file.warning(
                    f"Z Padding for box size is below the recommended "
                    f"{recommended_pad} Angstroms because pbc-override is enabled; "
                    "not increasing."
                )
            return PBC_pad_xy, PBC_pad_z

        if PBC_pad_xy < recommended_pad:
            my_log_file.warning(
                f"XY Padding for box size too small for non-bonded interactions, "
                f"increasing to {recommended_pad} Angstroms."
            )
            PBC_pad_xy = recommended_pad
        if PBC_pad_z < recommended_pad:
            my_log_file.warning(
                f"Z Padding for box size too small for non-bonded interactions, "
                f"increasing to {recommended_pad} Angstroms."
            )
            PBC_pad_z = recommended_pad

        return PBC_pad_xy, PBC_pad_z

    def _sheet_settings(self):

        '''
        Resolve the "sheets" block. A build is lattice matched when the block is
        present or system.sheets or system.tube-sheets is true; the box-end
        sheets are added only when system.sheets is true, the tube-end sheets
        only when system.tube-sheets is true.

        Parameters:
        ----------
        None

        Returns:
        -------
        sheet_cfg : dict or None
            The block merged over SHEET_DEFAULTS, or None when neither is set.

        Raises:
        ------
        None
        '''

        block = self.json_input.get('sheets')
        if (block is None and self.json_input['system'].get('sheets', False) is not True
                and not self._tube_sheets_on()):
            return None
        return Sheet_Builder.resolve_config(block)

    def _fit_box_to_sheet_lattice(self, requested_pad_xy=None):

        '''
        Scale the box x and y up to whole graphene cells about their centres,
        before any solvent is placed, and record the exact box. For an auto box
        this raises the xy padding on each axis just enough; an explicit size is
        scaled up the same way. Nothing is ever made smaller.

        Parameters:
        ----------
        requested_pad_xy : float, optional
            The xy padding of an auto box after _resolve_auto_box_padding, for the
            log and summary; None for an explicit size.

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        bond_length = self.config.get_default_species_bond_param("bond_length", 1.44)
        before = [[float(lo), float(hi)] for lo, hi in self.box_size]
        fitted = self.sheets.lattice_box_xy(before[0], before[1], bond_length,
                                            self._sheet_cfg["orientation"])
        self.box_size = [fitted["x"], fitted["y"], before[2]]
        self.config.fixed_box = [list(axis) for axis in self.box_size]

        raised_x = ((fitted["x"][1] - fitted["x"][0]) - (before[0][1] - before[0][0])) / 2.0
        raised_y = ((fitted["y"][1] - fitted["y"][0]) - (before[1][1] - before[1][0])) / 2.0
        self._sheet_lattice = {
            "nx": fitted["nx"], "ny": fitted["ny"],
            "period_x": fitted["period_x"], "period_y": fitted["period_y"],
            "lx": fitted["x"][1] - fitted["x"][0], "ly": fitted["y"][1] - fitted["y"][0],
            "raised_x": raised_x, "raised_y": raised_y,
            "requested_pad_xy": requested_pad_xy,
            "bond_length": bond_length,
        }

        if requested_pad_xy is not None:
            msg = (f"XY padding raised {requested_pad_xy:.2f} -> {requested_pad_xy + raised_x:.2f} Ang (x), "
                   f"{requested_pad_xy + raised_y:.2f} Ang (y) to fit "
                   f"{fitted['nx']} x {fitted['ny']} graphene cells")
        else:
            msg = (f"XY box scaled up by {2 * raised_x:.3f} Ang (x) and {2 * raised_y:.3f} Ang (y) "
                   f"to fit {fitted['nx']} x {fitted['ny']} graphene cells")
        my_log_file.info(f"{msg}: {self._sheet_lattice['lx']:.4f} x {self._sheet_lattice['ly']:.4f} Ang.")

    def _tube_sheets_on(self):

        return self.json_input['system'].get('tube-sheets', False) is True

    def _solvent_sides_used(self):

        if self.json_input['system'].get('solvation', False) is not True:
            return False
        return any(group.get('side', 'full') != 'full'
                   for group in self.json_input.get('solvation', {}).get('groups', []))

    def _side_bounds(self):

        '''
        The z boundaries that split the box into solvent sides: the tube-end
        sheet planes kept wall-gap clear when those sheets are on, otherwise the
        z range of the TUBE residue.

        Returns:
        -------
        bounds : tuple
            (lower, upper) bounding the - and + sides, and (inner_lo, inner_hi)
            bounding around-tube.
        '''

        frame = self._tube_frame
        if self._tube_sheets_on():
            clear = float(self._sheet_cfg['wall-gap'])
            return ((frame['z_lo'] - clear, frame['z_hi'] + clear),
                    (frame['z_lo'] + clear, frame['z_hi'] - clear))
        return (tuple(frame['tube_z']), tuple(frame['tube_z']))

    def _side_box(self, box, side):

        '''
        Cut a solvation box in z to one solvent side (see _side_bounds): "+"
        runs from the upper boundary to the box end, "-" from the box end to the
        lower boundary, and "around-tube" lies between the inner boundaries.

        Parameters:
        ----------
        box : list of three [lo, hi]
            The group's box, after any pad.
        side : str
            '+', '-' or 'around-tube'.

        Returns:
        -------
        box : list of three [lo, hi]
            The same x and y, the z range of the side.

        Raises:
        ------
        ValueError
            If the side leaves no room in z.
        '''

        (z_lo, z_hi), (inner_lo, inner_hi) = self._side_bounds()
        box_lo, box_hi = float(box[2][0]), float(box[2][1])
        if side == '+':
            z_range = [max(box_lo, z_hi), box_hi]
        elif side == '-':
            z_range = [box_lo, min(box_hi, z_lo)]
        else:
            z_range = [max(box_lo, inner_lo), min(box_hi, inner_hi)]
        if z_range[1] - z_range[0] <= 0.0:
            raise ValueError(f"Solvent side '{side}' leaves no room: z from {z_range[0]:.3f} "
                             f"to {z_range[1]:.3f} Ang.")
        return [list(box[0]), list(box[1]), z_range]

    def _side_passes(self, sol, box):

        '''
        The placements for one solvation group: (group, box) pairs. full and
        every ionic group are one pass over their box (an ionic "both" group
        avoids the middle through _ion_slab). A packed "both" group is placed
        twice, on the - side and on the + side, so a concentration counts each
        side's volume; an explicit count is split between them by volume.
        '''

        side = sol.get('side', 'full')
        if side == 'full' or (side == 'both' and sol.get('placement') == 'ionic'):
            return [(sol, box)]
        if side != 'both':
            return [(sol, self._side_box(box, side))]

        lower, upper = self._side_box(box, '-'), self._side_box(box, '+')
        if sol.get('count') is None:
            return [(sol, lower), (sol, upper)]
        span_lo, span_hi = lower[2][1] - lower[2][0], upper[2][1] - upper[2][0]
        total = int(sol['count'])
        count_lo = int(round(total * span_lo / (span_lo + span_hi)))
        parts = [(dict(sol, count=count_lo), lower), (dict(sol, count=total - count_lo), upper)]
        return [(part, part_box) for part, part_box in parts if part['count'] > 0]

    def _ion_slab(self, sol):

        '''
        The slab an ionic group must avoid, or None.

        - A "both" group avoids the whole middle between the side boundaries.
        - With tube-end sheets every other ionic group except an around-tube
          one avoids the slab between the sheet planes, widened by wall-gap,
          outside the tube radius, so the annulus carve never meets an ion and
          no ion is deleted.

        Returns:
        -------
        slab : tuple or None
            (z_lo, z_hi, axis_x, axis_y, radius) for Solution_Generator; a
            radius of -1 covers the whole slab.
        '''

        side = sol.get('side', 'full')
        if sol.get('placement') != 'ionic' or side == 'around-tube':
            return None
        frame = self._tube_frame
        if side == 'both':
            (z_lo, z_hi), _ = self._side_bounds()
            return (z_lo, z_hi, frame['axis'][0], frame['axis'][1], -1.0)
        if not self._tube_sheets_on():
            return None
        clear = float(self._sheet_cfg['wall-gap'])
        return (frame['z_lo'] - clear, frame['z_hi'] + clear,
                frame['axis'][0], frame['axis'][1], frame['radius'])

    def _drop_lumen_residues(self, residues_before):

        '''
        Remove the residues of the group just placed (those numbered above
        residues_before) that have an atom inside the tube radius, so an
        around-tube group keeps out of the lumen. Ions are never removed; the
        ionic no-go cuboid already keeps them out of the tube.
        '''

        positions, indexes, types, residues = self.force_field_manager.get_atom_positions_all()
        positions = np.asarray(positions, dtype=float)
        residues = np.asarray(residues)
        frame = self._tube_frame
        ion_types = set(self.config.ions) | set(self.config.dions) | set(self.config.dion_particles)

        r = np.hypot(positions[:, 0] - frame['axis'][0], positions[:, 1] - frame['axis'][1])
        inside = set(residues[(residues > residues_before) & (r < frame['radius'])].tolist())
        drop = set()
        for residue in inside:
            members = np.nonzero(residues == residue)[0]
            if not any(types[i] in ion_types for i in members):
                drop.update(int(indexes[i]) for i in members)
        if drop:
            self.force_field_manager.obliterate_atoms(drop)
            my_log_file.info(f"around-tube: removed {len(drop)} atoms of molecules placed inside the tube.")

    def _carve_for_tube_sheets(self, sheets):

        '''
        Remove the solvent the tube-end sheets would enclose or overlap (see
        Sheet_Builder.annulus_residues) and list it in <system>.carved. Raises
        rather than delete an ion, so the charge is never changed.

        Parameters:
        ----------
        sheets : dict
            From Sheet_Builder.tube_sheet_sites.

        Returns:
        -------
        None

        Raises:
        ------
        ValueError
            If an ion lies where the carve would remove it.
        '''

        carved = self.sheets.annulus_residues(
            sheets, self._tube_frame, float(self._sheet_cfg['wall-gap']), self._around_residues)
        ions = [entry for entry in carved if entry['ion']]
        if ions:
            listed = ", ".join(f"residue {entry['residue']} ({entry['type']})" for entry in ions)
            raise ValueError(f"Tube-end sheets: ions lie where the annulus is cleared: {listed}. "
                             f"Ions are never deleted; check the ionic groups' side.")

        path = os.path.join(self.config.system_folder, f"{self.config.system_name}.carved")
        with open(path, "w") as f:
            f.write("# Solvent residues removed before the tube-end sheets went in: those with an atom\n")
            f.write("# between the sheet planes outside the tube, or within wall-gap of a sheet carbon.\n")
            f.write("# Residue numbers are the internal ones before the removal.\n")
            f.write("# residue  first-atom-type   atoms      mean-z\n")
            for entry in carved:
                f.write(f"{entry['residue']:>9d}  {entry['type']:<16s} {len(entry['atoms']):>6d}  "
                        f"{entry['z']:>10.4f}\n")

        atoms = {atom for entry in carved for atom in entry['atoms']}
        if atoms:
            self.force_field_manager.obliterate_atoms(atoms)
        self._carved_count = len(carved)
        my_log_file.info(f"Tube-end sheets: removed {len(carved)} solvent residues "
                         f"({len(atoms)} atoms); listed in {path}.")


    # -------------------------------------------------------------------------
    # VALIDATION HELPERS
    # -------------------------------------------------------------------------

    def _validate_savepath(self, og_path):

        '''
        Return the file path to use for saving an output file. Currently acts
        as a passthrough; intended as an override point for subclasses that need
        to redirect output to a different location.

        Parameters:
        ----------
        og_path : str
            The original proposed output file path.

        Returns:
        -------
        og_path : str
            The same path, unchanged.

        Raises:
        ------
        None
        '''

        return og_path 

    def _validate_input(self, force_field, nanotube, functionalisation, solvation):

        '''
        Validate all user-facing fields in self.json_input against the lists of
        supported options declared in Config.  Checks are ordered:

            1. Force-field name -> Config.supported_force_fields
            (must pass before any file-existence check, because the path
            depends on the name)
            2. Force-field folder exists on disk
            3. chirality_kit_fragments.xdata and chirality_kit_solvents.xdata
            exist inside that folder
            4. Every functional-group type -> fragments file section headers
            5. Every solvent / ion type    -> solvents file section headers

        Raises
        ------
        ValueError
            If any named option is not in the supported list.
        FileNotFoundError
            If the force-field folder or its data files are missing.
        '''

        ff = self.json_input['system']['force-field'].lower()

        # ------------------------------------------------------------------
        # 0. Optional top-level seed
        # ------------------------------------------------------------------
        seed_value = self.json_input.get('seed')
        if seed_value is not None:
            try:
                int(seed_value)
            except (TypeError, ValueError):
                msg = (
                    f"\n{'='*60}\n"
                    f"  VALIDATION ERROR: Seed\n"
                    f"{'='*60}\n"
                    f"  Attempted  : {seed_value!r}\n"
                    f"  Status     : NOT AN INTEGER\n"
                    f"  Expected   : an integer, or omit the key for a random seed\n"
                    f"{'='*60}"
                )
                my_log_file.error(msg)
                raise ValueError(msg)

        if nanotube:
            nt = self.json_input.get('nanotube', {})
            n, m, repeats = (nt.get(key) for key in ('n', 'm', 'repeats'))
            valid_indices = (
                type(n) is int
                and type(m) is int
                and n >= m >= 0
                and n + m >= 3
            )
            if not valid_indices or type(repeats) is not int or repeats < 1:
                msg = (
                    f"\n{'='*60}\n"
                    f"  VALIDATION ERROR: Nanotube Geometry\n"
                    f"{'='*60}\n"
                    f"  Attempted  : ({n}, {m}) x {repeats}\n"
                    f"  Status     : IMPOSSIBLE NANOTUBE\n"
                    f"  Expected   : integer n >= m >= 0, n + m >= 3,\n"
                    f"               and integer repeats >= 1\n"
                    f"{'='*60}"
                )
                my_log_file.error(msg)
                raise ValueError(msg)

        # ------------------------------------------------------------------
        # 1. Force-field name
        # ------------------------------------------------------------------
        available_ffs = [d for d in os.listdir(self.config.force_fields_root)
                        if os.path.isdir(os.path.join(self.config.force_fields_root, d))]

        if ff.lower() not in [f.lower() for f in available_ffs]:
            msg = (
                f"\n{'='*60}\n"
                f"  VALIDATION ERROR: Force Field\n"
                f"{'='*60}\n"
                f"  Attempted  : {ff}\n"
                f"  Status     : NOT FOUND\n"
                f"  Available  : {', '.join(sorted(available_ffs))}\n"
                f"{'='*60}"
            )
            my_log_file.error(msg)
            raise ValueError(msg)

        # ------------------------------------------------------------------
        # 2. Force-field folder
        # ------------------------------------------------------------------
        ff_folder = os.path.join(self.config.force_fields_root, ff.lower())

        if not os.path.isdir(ff_folder):
            msg = (
                f"\n{'='*60}\n"
                f"  VALIDATION ERROR: Force Field Folder\n"
                f"{'='*60}\n"
                f"  Attempted  : {ff}\n"
                f"  Status     : FOLDER NOT FOUND\n"
                f"  Expected   : {ff_folder}\n"
                f"{'='*60}"
            )
            my_log_file.error(msg)
            raise FileNotFoundError(msg)

        # ------------------------------------------------------------------
        # 3. Data files inside that folder
        # ------------------------------------------------------------------
        fragments_path = os.path.join(ff_folder, "chirality_kit_fragments.xdata")
        solvents_path  = os.path.join(ff_folder, "chirality_kit_solvents.xdata")

        if nanotube:
            if not os.path.isfile(fragments_path):
                msg = (
                    f"\n{'='*60}\n"
                    f"  VALIDATION ERROR: Fragments File\n"
                    f"{'='*60}\n"
                    f"  Force field : {ff}\n"
                    f"  Status      : FRAGMENTS FILE NOT FOUND\n"
                    f"  Expected    : {fragments_path}\n"
                    f"{'='*60}"
                )
                my_log_file.error(msg)
                raise FileNotFoundError(msg)

        if solvation:
            if not os.path.isfile(solvents_path):
                msg = (
                    f"\n{'='*60}\n"
                    f"  VALIDATION ERROR: Solvents File\n"
                    f"{'='*60}\n"
                    f"  Force field : {ff}\n"
                    f"  Status      : SOLVENTS FILE NOT FOUND\n"
                    f"  Expected    : {solvents_path}\n"
                    f"{'='*60}"
                )
                my_log_file.error(msg)
                raise FileNotFoundError(msg)

        # ------------------------------------------------------------------
        # 4. check functional-group names
        # ------------------------------------------------------------------
        nt_type = self.json_input.get('nanotube', {}).get('type')
        if nanotube and nt_type in {"MoSSeNT", "MoSSe"}:
            outer_chalcogen = self.json_input.get('nanotube', {}).get('outer-chalcogen', 'Se')
            if outer_chalcogen not in {"S", "Se"}:
                msg = (
                    f"\n{'='*60}\n"
                    f"  VALIDATION ERROR: MoSSeNT\n"
                    f"{'='*60}\n"
                    f"  outer-chalcogen : {outer_chalcogen}\n"
                    f"  Status          : INVALID\n"
                    f"  Expected        : 'S' or 'Se'\n"
                    f"{'='*60}"
                )
                my_log_file.error(msg)
                raise ValueError(msg)

        if nanotube and functionalisation and nt_type in {"BNNT", "MoS2NT", "MoS2", "MoSSeNT", "MoSSe"}:
            msg = (
                f"\n{'='*60}\n"
                f"  VALIDATION ERROR: Functionalisation\n"
                f"{'='*60}\n"
                f"  Nanotube type : {nt_type}\n"
                f"  Status        : NOT IMPLEMENTED\n"
                f"  Reason        : multi-species functionalisation requires\n"
                f"                  material-aware attachment-site logic.\n"
                f"{'='*60}"
            )
            my_log_file.error(msg)
            raise NotImplementedError(msg)

        if self.json_input['system'].get('functionalisation', False):
            available_fragments = set()
            with open(fragments_path, 'r') as f:
                for line in f:
                    stripped = line.strip()
                    if stripped.startswith("Fragment "):
                        available_fragments.add(stripped.split(None, 1)[1].strip())

            func = self.json_input.get('functionalisation', {})

            # Terminating Groups
            for group in func.get('term', {}).get('groups', []):
                fg_type = group['type']
                count = int(group['count'])
                unsupported = set(group) - TERM_GROUP_FIELDS
                if unsupported:
                    raise ValueError(
                        f"Unsupported term group field(s): "
                        f"{', '.join(sorted(unsupported))}."
                    )
                side = group.get('side', 'both')
                if side not in TERM_SIDES:
                    raise ValueError(
                        f"Invalid term side {side!r}; expected 'both', '-', or '+'."
                    )

                if (fg_type not in VACANT_FG_TYPES
                        and fg_type not in available_fragments
                        and count > 0):
                    msg = (
                        f"\n{'='*60}\n"
                        f"  VALIDATION ERROR: Functional Group Term\n"
                        f"{'='*60}\n"
                        f"  Force field : {ff}\n"
                        f"  Attempted   : {fg_type}\n"
                        f"  Status      : NOT FOUND in fragments file\n"
                        f"  Available   : {', '.join(sorted(available_fragments))}\n"
                        f"{'='*60}"
                    )
                    my_log_file.error(msg)
                    raise ValueError(msg)
                
            # "DUMMY" is the name VACANT had before 0.8.0; it is read the same way.
            legacy = sum(1 for g in func.get('term', {}).get('groups', [])
                         if g.get('type') == LEGACY_VACANT_FG_TYPE)
            legacy += sum(1 for ring in func.get('rings', []) for g in ring.get('groups', [])
                          if g.get('type') == LEGACY_VACANT_FG_TYPE)
            if legacy:
                my_log_file.warning(
                    f'"{LEGACY_VACANT_FG_TYPE}" is now called "{VACANT_FG_TYPE}"; the input is read '
                    f'the same way and builds the same system. Use "{VACANT_FG_TYPE}" in new inputs.')

            # hydrogenation needs the H-term fragment
            if func.get('term', {}).get('term-hydrogenate', False):
                if 'H-term' not in available_fragments:
                    msg = (
                        f"\n{'='*60}\n"
                        f"  VALIDATION ERROR: Term Hydrogenate\n"
                        f"{'='*60}\n"
                        f"  Force field : {ff}\n"
                        f"  Status      : 'H-term' fragment NOT FOUND\n"
                        f"  Reason      : term-hydrogenate is True but this force\n"
                        f"                field has no 'H-term' fragment defined.\n"
                        f"  Fix         : set term-hydrogenate to False, or add an\n"
                        f"                'H-term' Fragment entry to:\n"
                        f"                {fragments_path}\n"
                        f"  Available   : {', '.join(sorted(available_fragments))}\n"
                        f"{'='*60}"
                    )
                    my_log_file.error(msg)
                    raise ValueError(msg)

            # Ring Groups
            for ring in func.get('rings', []):
                for group in ring.get('groups', []):
                    fg_type = group['type']
                    count = int(group['count'])

                    if (fg_type not in VACANT_FG_TYPES
                            and fg_type not in available_fragments
                            and count > 0):
                        msg = (
                            f"\n{'='*60}\n"
                            f"  VALIDATION ERROR: Functional Group Ring\n"
                            f"{'='*60}\n"
                            f"  Force field : {ff}\n"
                            f"  Attempted   : {fg_type}\n"
                            f"  Status      : NOT FOUND in fragments file\n"
                            f"  Available   : {', '.join(sorted(available_fragments))}\n"
                            f"{'='*60}"
                        )
                        my_log_file.error(msg)
                        raise ValueError(msg)

            # Loose Groups
            for group in func.get('loose', {}).get('groups', []):
                fg_type = group['type']
                count = int(group.get('count', 0))
                mode = str(group.get('loose-placement', 'random'))

                z_from = float(group.get('z-from', 0.0))
                z_to   = float(group.get('z-to', 1.0))
                if not (0.0 <= z_from < z_to <= 1.0):
                    msg = (
                        f"\n{'='*60}\n"
                        f"  VALIDATION ERROR: Functional Group Loose\n"
                        f"{'='*60}\n"
                        f"  Attempted   : {fg_type}\n"
                        f"  Status      : INVALID z-window\n"
                        f"  Reason      : need 0 <= z-from < z-to <= 1\n"
                        f"                (got z-from={z_from}, z-to={z_to})\n"
                        f"{'='*60}"
                    )
                    my_log_file.error(msg)
                    raise ValueError(msg)

                if fg_type not in available_fragments and (count > 0 or mode.startswith('all')):
                    msg = (
                        f"\n{'='*60}\n"
                        f"  VALIDATION ERROR: Functional Group Loose\n"
                        f"{'='*60}\n"
                        f"  Force field : {ff}\n"
                        f"  Attempted   : {fg_type}\n"
                        f"  Status      : NOT FOUND in fragments file\n"
                        f"  Available   : {', '.join(sorted(available_fragments))}\n"
                        f"{'='*60}"
                    )
                    my_log_file.error(msg)
                    raise ValueError(msg)
                

        # ------------------------------------------------------------------
        # 5. check solvent and ion names
        # ------------------------------------------------------------------
        if self.json_input['system'].get('solvation', False):
            available_solvents = set()
            with open(solvents_path, 'r') as f:
                for line in f:
                    stripped = line.strip()
                    if stripped.startswith("Solvent "):
                        available_solvents.add(stripped.split(None, 1)[1].strip())

            for group in self.json_input.get('solvation', {}).get('groups', []):
                sol_type = group['type']

                if group.get('placement') == 'ionic':
                    # split the salt formula into its ions
                    cathode_type, _, anode_type, _ = self.solvation._process_salt_formula(sol_type)
                    ions_to_check = [cathode_type, anode_type]
                else:
                    ions_to_check = [sol_type]

                for ion in ions_to_check:
                    if ion not in available_solvents:
                        if group.get('placement') == 'ionic':
                            attempted_str = f"{sol_type}  ->  checking ion '{ion}'"
                        else:
                            attempted_str = ion

                        msg = (
                            f"\n{'='*60}\n"
                            f"  VALIDATION ERROR: Solvent / Ion\n"
                            f"{'='*60}\n"
                            f"  Force field : {ff}\n"
                            f"  Attempted   : {attempted_str}\n"
                            f"  Status      : NOT FOUND in solvents file\n"
                            f"  Available   : {', '.join(sorted(available_solvents))}\n"
                            f"{'='*60}"
                        )
                        my_log_file.error(msg)
                        raise ValueError(msg)

            # ------------------------------------------------------------------
            # 6. pick the bulk phase
            # ------------------------------------------------------------------
            # Only the highest-priority packed group skips self-overlap checks.
            # Ion groups still need the check, and tied bulk priorities are ambiguous.
            solvation_groups = sorted(
                self.json_input.get('solvation', {}).get('groups', []),
                key=lambda g: (int(g['priority']), str(g['type'])),
            )
            packed_groups = [g for g in solvation_groups if g.get('placement') == 'packed']

            if not packed_groups:
                self.config.bulk_priority = None
                my_log_file.info(
                    "Bulk phase: none; no packed solvation group, so every overlapping "
                    "pair is checked and resolved by priority."
                )
            else:
                bulk_priority = max(int(g['priority']) for g in packed_groups)
                holders = [g for g in solvation_groups if int(g['priority']) == bulk_priority]

                if len(holders) > 1:
                    listed = ", ".join(
                        f"'{g['type']}' ({g.get('placement', 'packed')})" for g in holders
                    )
                    msg = (
                        f"\n{'='*60}\n"
                        f"  VALIDATION ERROR: Bulk Phase\n"
                        f"{'='*60}\n"
                        f"  Priority    : {bulk_priority}  (the highest in the solvation list)\n"
                        f"  Held by     : {listed}\n"
                        f"  Status      : AMBIGUOUS; more than one group claims the bulk\n"
                        f"  Reason      : the highest-priority packed group is the bulk\n"
                        f"                phase, and its molecules are exempt from the\n"
                        f"                overlap check because Bridson placement already\n"
                        f"                sets their spacing. Only one group can hold that\n"
                        f"                exemption: shared, the two would never be checked\n"
                        f"                against each other and would interpenetrate.\n"
                        f"                An ionic group can never hold it; ionic placement\n"
                        f"                does not guarantee its own spacing.\n"
                        f"  Fix         : give all but one of them a lower priority number.\n"
                        f"                A lowered group is culled against the bulk, so keep\n"
                        f"                it dilute; a dense second solvent would need a true\n"
                        f"                co-bulk phase, which is not implemented.\n"
                        f"{'='*60}"
                    )
                    my_log_file.error(msg)
                    raise ValueError(msg)

                self.config.bulk_priority = bulk_priority
                my_log_file.info(
                    f"Bulk phase: '{holders[0]['type']}' at priority {bulk_priority} (packed). "
                    f"Its molecules are exempt from the overlap check; Bridson placement sets "
                    f"their spacing. Every other pair is checked and resolved by priority."
                )

        # ------------------------------------------------------------------
        # 7. Sheets
        # ------------------------------------------------------------------
        self._validate_sheets()
        self._validate_solvent_sides()

    def _validate_sheets(self):

        '''
        Validate system.sheets, system.tube-sheets and the optional "sheets"
        block: a carbon nanotube that is not periodic, lattice-xy on, a
        supported orientation, positive gaps, a hole-buffer of zero or more, a
        numeric rim-offset, and
        four-character residue names, distinct for the sheets in use, that clash
        with no other residue name, segment name or ion name and do not contain
        END.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        ValueError
            On the first group of problems found, listed together.
        '''

        system = self.json_input['system']
        enabled = system.get('sheets', False)
        tube = system.get('tube-sheets', False)
        block = self.json_input.get('sheets')
        problems = []

        if not isinstance(enabled, bool):
            problems.append(f"system.sheets must be true or false, not {enabled!r}")
        if not isinstance(tube, bool):
            problems.append(f"system.tube-sheets must be true or false, not {tube!r}")
        if block is None and enabled is not True and tube is not True:
            if problems:
                self._raise_sheet_error(problems)
            return
        if block is not None and not isinstance(block, dict):
            self._raise_sheet_error([f'"sheets" must be an object, not {type(block).__name__}'])

        unknown = sorted(set(block or {}) - set(SHEET_DEFAULTS))
        if unknown:
            problems.append(f"unknown key(s) {', '.join(unknown)}; "
                            f"allowed: {', '.join(sorted(SHEET_DEFAULTS))}")
        cfg = Sheet_Builder.resolve_config(block)

        nanotube_type = self.json_input.get('nanotube', {}).get('type')
        if system.get('nanotube') is not True or nanotube_type not in {"CNT", "CNT-Drude"}:
            problems.append("sheets need a carbon nanotube (nanotube type CNT or CNT-Drude); "
                            "the sheet carbons use the tube's carbon species")
        if self.json_input.get('nanotube', {}).get('periodic', False) is True:
            problems.append('a periodic nanotube cannot take sheets: remove system.sheets, '
                            'system.tube-sheets and the "sheets" block, or set nanotube.periodic '
                            'to false')
        if cfg['lattice-xy'] is not True:
            problems.append("lattice-xy must be true: a periodic sheet needs a whole number "
                            "of graphene cells in x and y")
        if cfg['orientation'] not in SHEET_ORIENTATIONS:
            problems.append(f"orientation {cfg['orientation']!r} is not one of "
                            f"{', '.join(SHEET_ORIENTATIONS)}")
        for key in ('wall-gap', 'seam-gap'):
            value = cfg[key]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
                problems.append(f"{key} must be a positive number of Angstrom, not {value!r}")
        value = cfg['hole-buffer']
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0:
            problems.append(f"hole-buffer must be zero or a positive number of Angstrom, not {value!r}")
        value = cfg['rim-offset']
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            problems.append(f"rim-offset must be a number of Angstrom (positive moves the tube-end "
                            f"sheets outward), not {value!r}")

        resnames = cfg['resnames']
        all_sides = SHEET_SIDES + TUBE_SHEET_SIDES
        if not isinstance(resnames, dict) or set(resnames) != set(all_sides):
            problems.append(f"resnames must name exactly {', '.join(all_sides)}")
        else:
            # The PDB and PSF also use SYS as a segment name and the ion names
            # as residue names.
            reserved = {"TUBE", "IONS", "SYS"}
            reserved.update(str(name).upper() for name in self.config.ions)
            reserved.update(str(name).upper() for name in self.config.dion_pdb_names.values())
            for sol in self.json_input.get('solvation', {}).get('groups', []):
                label = re.sub(r"[^A-Za-z0-9]", "", str(sol.get('type', ''))).upper()
                if label:
                    reserved.add(label[:4])
            for side in all_sides:
                name = resnames[side]
                if not isinstance(name, str) or not re.fullmatch(r"[A-Z0-9]{1,4}", name):
                    problems.append(f"resnames.{side} {name!r} must be 1-4 uppercase letters or digits")
                elif name in reserved:
                    problems.append(f"resnames.{side} {name!r} clashes with an existing residue name")
                elif "END" in name:
                    problems.append(f"resnames.{side} {name!r} contains END, which ends the "
                                    f"atom records when restraints.ref is written")
            in_use = SHEET_SIDES + (TUBE_SHEET_SIDES if tube is True else ())
            if len({resnames[side] for side in in_use}) != len(in_use):
                problems.append("the sheets need different resnames")

        if problems:
            self._raise_sheet_error(problems)

    def _validate_solvent_sides(self):

        '''
        Validate each solvation group's side, in the same terms as the terminal
        groups' side: full (the default, the whole box), "-" or "+" (one side
        of the tube), "both" (the two sides, nothing around or inside the tube)
        or around-tube. Every side but full needs a nanotube. A one-sided group in a box with no walls at all is allowed,
        with a warning, since its two sides meet across the periodic boundary.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        ValueError
            On the first group of problems found, listed together.
        '''

        system = self.json_input['system']
        if system.get('solvation', False) is not True:
            return
        allowed = ("full", "-", "+", "both", "around-tube")
        walled = system.get('sheets', False) is True or system.get('tube-sheets', False) is True
        problems = []
        for group in self.json_input.get('solvation', {}).get('groups', []):
            side = group.get('side', 'full')
            label = f"solvation group '{group.get('type')}'"
            if side not in allowed:
                problems.append(f"{label}: side {side!r} is not one of {', '.join(allowed)}")
            elif side != 'full' and system.get('nanotube') is not True:
                problems.append(f"{label}: side {side!r} needs a nanotube to measure from; use full")
            elif side in ('-', '+') and not walled:
                my_log_file.warning(
                    f"Solvation group '{group.get('type')}' is placed on the {side} side, but the box "
                    f"has no walls, so the two sides meet across the periodic boundary."
                )

        if problems:
            lines = "\n".join(f"  - {problem}" for problem in problems)
            msg = (
                f"\n{'='*60}\n"
                f"  VALIDATION ERROR: Solvent side\n"
                f"{'='*60}\n"
                f"{lines}\n"
                f"{'='*60}"
            )
            my_log_file.error(msg)
            raise ValueError(msg)

    def _raise_sheet_error(self, problems):

        lines = "\n".join(f"  - {problem}" for problem in problems)
        msg = (
            f"\n{'='*60}\n"
            f"  VALIDATION ERROR: Sheets\n"
            f"{'='*60}\n"
            f"{lines}\n"
            f"{'='*60}"
        )
        my_log_file.error(msg)
        raise ValueError(msg)

    def _validate_fg_available_carbons(self, fg_config, n_atoms):

        '''
        Validate that there are enough available carbons on the nanotube to
        '''

        n_fgs = 0 

        # Terminating Groups
        for group in fg_config.get('term', {}).get('groups', []):
            fg_type = group['type']
            count = int(group['count'])
            n_fgs += count

        # Ring Groups
        for ring in fg_config.get('rings', []):
            for group in ring.get('groups', []):
                fg_type = group['type']
                count = int(group['count'])
            n_fgs += count

        # all-* placement is bounded by atom count, so it doesn't use this budget
        for group in fg_config.get('loose', {}).get('groups', []):
            fg_type = group['type']
            if str(group.get('loose-placement', 'random')).startswith('all'):
                continue
            count = int(group.get('count', 0))
            n_fgs += count

        if n_fgs > n_atoms:
            msg = (
                f"\n{'='*60}\n"
                f"  VALIDATION ERROR: Functional Group Count\n"
                f"{'='*60}\n"
                f"  Attempted   : {n_fgs} functional groups\n"
                f"  Available   : {n_atoms} carbon atoms in the nanotube\n"
                f"  Status      : NOT ENOUGH AVAILABLE CARBONS\n"
                f"  Fix         : reduce the total count of functional groups, or\n"
                f"                increase the nanotube size (n, m, repeats)\n"
                f"{'='*60}"
            )
            my_log_file.error(msg)
            raise ValueError(msg)

        return True 
   
   
    # -------------------------------------------------------------------------
    # MD DATA FILES
    # -------------------------------------------------------------------------

    def write_data_file(self):

        '''
        Write a LAMMPS .data file from the main xdata file by stripping all
        comment fields and non-ASCII characters. The result is saved alongside
        the other output files in the system folder.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        src = self.config.main_file
        lammps_file = src[:-6] + ".data"
        lammps_file = self._validate_savepath(lammps_file)

        shutil.copy(src, lammps_file)

        with open(lammps_file, "r", encoding="utf-8", errors="ignore") as f:
            text = f.read()

        text = re.sub(r"#.*$", "", text, flags=re.MULTILINE)

        text = "".join(
            ch for ch in text 
            if ch in ("\t", "\n", "\r") or 0x20 <= ord(ch) <= 0x7E
        )

        if text and not text.endswith("\n"):
            text += "\n"

        with open(lammps_file, "w", encoding="utf-8") as f:
            f.write(text)

        my_log_file.info(f"File type .data file written to {lammps_file}")

    def write_xsc_file(self):

        '''
        Write a NAMD extended-system (.xsc) file holding the exact box of a sheet
        or lattice-matched build: step 0, the three cell vectors and the cell
        centre as origin. NAMD reads the periodic cell from it, which a periodic
        sheet needs exactly; the PDB atom extents would be too small.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        src = self.config.main_file
        xsc_file = self._validate_savepath(src[:-6] + ".xsc")
        (x_lo, x_hi), (y_lo, y_hi), (z_lo, z_hi) = (
            (float(lo), float(hi)) for lo, hi in self.config.fixed_box
        )
        lx, ly, lz = x_hi - x_lo, y_hi - y_lo, z_hi - z_lo
        ox, oy, oz = (x_lo + x_hi) / 2.0, (y_lo + y_hi) / 2.0, (z_lo + z_hi) / 2.0

        with open(xsc_file, "w") as f:
            f.write("# NAMD extended system configuration written by chirality-kit\n")
            f.write("#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z\n")
            f.write(f"0 {lx:.6f} 0 0 0 {ly:.6f} 0 0 0 {lz:.6f} {ox:.6f} {oy:.6f} {oz:.6f}\n")

        my_log_file.info(f"File type .xsc file written to {xsc_file}")

    def write_xyz_file(self): 

        '''
        Write an XYZ-format coordinate file from the Atoms section of the main
        file. Each atom's element symbol is resolved from its type name and the
        output includes the standard two-line XYZ header.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        src = self.config.main_file
        filename = src[:-6] + ".xyz"
        filename = self._validate_savepath(filename)

        position_lines = self.file_manager.find_and_extract_lines(self.config.main_file, "Atoms", header=False)

        with open(filename, 'w') as f:
            f.write(f"{len(position_lines)}\n")
            f.write("XYZ generated with chirality_kit\n")
            for line in position_lines:
                parts = self.file_manager.line_to_array(line)
                x, y, z = float(parts[4]), float(parts[5]), float(parts[6])

                atom_type = parts[-1]
                element = self.force_field_manager.get_atom_element(atom_type)  
                f.write(f"{element:<3s} {x:>12.5f} {y:>12.5f} {z:>12.5f}\n")

        my_log_file.info(f"File type .xyz file written to {filename}")

    def write_pdb_file(self):

        '''
        Write a PDB-format coordinate file from the main xdata file. Residue
        names, segment IDs, and B-factor fields are assigned based on atom type:
        nanotube carbons receive residue 'TUBE' and B-factor 1.00; ions keep
        their CHARMM atom/residue names and use their net charge as the B-factor;
        Drude particles and lone pairs use B-factor 0.00.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        src = self.config.main_file
        pdb_file = src[:-6] + ".pdb"
        pdb_file = self._validate_savepath(pdb_file)

        position_lines = self.file_manager.find_and_extract_lines(self.config.main_file, "Atoms", header=False)

        with open(pdb_file, 'w') as f:
            f.write(f"REMARK  GENERATED WITH CHIRALITY-KIT\n")
            f.write(f"REMARK  SETUP PERIODIC BOUNDARY CONDITION\n")
            f.write(f"REMARK  DATE {time.strftime('%Y-%m-%d')}\n")
            if self.config.fixed_box is not None:
                lx, ly, lz = (float(hi) - float(lo) for lo, hi in self.config.fixed_box)
                f.write(f"CRYST1{lx:9.3f}{ly:9.3f}{lz:9.3f}{90.0:7.2f}{90.0:7.2f}{90.0:7.2f} P 1           1\n")
            
            tube_species_counts = defaultdict(int)
            atom_count = 1
            previous_core_props = None

            for line_number, line in enumerate(position_lines):
                parts = self.file_manager.line_to_array(line)
                atom_index, residue_index, charge, atom_type = int(parts[0]), int(parts[1]), float(parts[3]), parts[-1]
                x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
                original_atom_type = atom_type

                atom_type_parts = atom_type.split('-')
                is_drude_particle = (
                    (len(atom_type_parts) > 1 and atom_type_parts[1] == "DRUD")
                    or atom_type in self.config.dion_particles
                )
                is_ion_core = atom_type in self.config.ions or atom_type in self.config.dions
                is_ion_particle = atom_type in self.config.dion_particles

                if len(atom_type_parts) > 1 and atom_type_parts[1] != "DRUD":
                    element = self.force_field_manager.get_atom_element(atom_type)
                    res_name = atom_type_parts[1]
                    atom_type = atom_type_parts[0]

                elif len(atom_type_parts) > 1 and atom_type_parts[1] == "DRUD":
                    element = self.force_field_manager.get_atom_element(atom_type)
                    res_name = "IONS" if atom_type_parts[0] in self.config.dions else "TUBE"
                    atom_type = "DRUD"

                else:
                    element = self.force_field_manager.get_atom_element(atom_type)
                    if (atom_type in self.config.ions
                            or atom_type in self.config.dions
                            or atom_type in self.config.dion_particles):
                        res_name = "IONS"
                    else:
                        res_name = "TUBE"

                # Sheet residues carry their own name; their Drudes follow the core.
                sheet_name = self.config.residue_name_by_index.get(residue_index)
                if sheet_name is not None:
                    res_name = sheet_name

                if is_drude_particle and previous_core_props is not None:
                    residue_index, res_name = previous_core_props

                atom_name = element # for other atoms

                species_entry = self.config.default_species_by_type.get(original_atom_type)
                if (species_entry is not None
                        and not self.config._default_species_entry_is_drude_particle(species_entry)):
                    count = tube_species_counts[original_atom_type]
                    atom_name = species_entry["atom_element"] + str(self._convert_to_base36(count))
                    tube_species_counts[original_atom_type] += 1

                if is_ion_core:
                    atom_name = self.config.dion_pdb_names.get(original_atom_type, original_atom_type)
                    res_name = atom_name
                elif is_ion_particle:
                    atom_name = original_atom_type

                if sheet_name is not None:
                    seg_id = sheet_name
                else:
                    seg_id = "TUBE" if res_name == "TUBE" else "SYS"

                occupancy = 1.00
                temp = 0.00 # default temp factor

                if is_ion_core:
                    temp = charge
                    if line_number + 1 < len(position_lines):
                        next_parts = self.file_manager.line_to_array(position_lines[line_number + 1])
                        if (int(next_parts[1]) == residue_index
                                and next_parts[-1] in self.config.dion_particles):
                            temp += float(next_parts[3])
                elif res_name == "TUBE" or sheet_name is not None:
                    temp = 1.00    # set temp factor to 1.00 for tube and sheet carbons

                if element in ['DP', 'LP']: # lone pairs or drude particles
                    temp = 0.00

                f.write(
                    "ATOM  {atom_index:5d} {atom_name:^4s} {res_name:>4s} {residue_index:4d}    "
                    "{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{temp:6.2f}      {seg_id:>5s}\n"
                    .format(
                        atom_index=atom_count,
                        atom_name=atom_name,
                        element=element,
                        res_name=res_name,
                        residue_index=residue_index,
                        x=x, y=y, z=z,
                        occupancy=occupancy,
                        temp=temp,
                        seg_id=seg_id
                    )
                )

                atom_count += 1
                if not is_drude_particle:
                    previous_core_props = (residue_index, res_name)

            f.write("END\n")

        my_log_file.info(f"File type .pdb file written to {pdb_file}")

    def write_psf_file(self): 

        '''
        Write a CHARMM PSF file from the main xdata file. Handles both standard
        and Drude-polarisable topologies (EXT CMAP vs EXT CMAP DRUDE header).
        Sections written: NATOM, NBOND, NTHETA, NPHI, NIMPHI, NDON, NACC, NNB,
        NUMLP/NUMLPH (lone pairs), NUMANISO (anisotropy), and NCRTERM. Raises
        ValueError if a Drude bond coefficient cannot be found for an anisotropic atom.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        ValueError
            If the Drude-Drude bond coefficient for an anisotropic atom type
            cannot be found in the Bond Coeffs section.
        '''

        src = self.config.main_file
        psf_file = src[:-6] + ".psf"
        psf_file = self._validate_savepath(psf_file)
        atom_index_corrector = []
        self._psf_charge_neutrality_report = None

        def _base_atom_type(atom_type):
            atom_type_parts = str(atom_type).split('-')
            return atom_type_parts[0], atom_type_parts

        def _is_drude_particle(atom_type):
            base_atom_type, atom_type_parts = _base_atom_type(atom_type)
            return (
                base_atom_type == "DRUD"
                or base_atom_type in self.config.dion_particles
                or (len(atom_type_parts) > 1 and atom_type_parts[1] == "DRUD")
            )

        def _is_ion_core(atom_type):
            base_atom_type, _ = _base_atom_type(atom_type)
            return (
                not _is_drude_particle(atom_type)
                and (
                    base_atom_type in self.config.ions
                    or base_atom_type in self.config.dions
                )
            )

        def _same_sign(a, b):
            return (a > 0 and b > 0) or (a < 0 and b < 0)

        def _make_psf_charge_report(record, original_charge,
                                    adjusted_charge, original_total,
                                    adjusted_total):
            return {
                "applied": True,
                "record_index": record["record_index"],
                "psf_atom_index": record["psf_atom_index"],
                "xdata_atom_index": record["xdata_atom_index"],
                "residue_index": record["residue_index"],
                "atom_type": record["atom_type"],
                "original_charge": float(original_charge),
                "adjusted_charge": float(adjusted_charge),
                "original_total": float(original_total),
                "adjusted_total": float(adjusted_total),
                "snapshot": [],
            }

        def _make_skipped_psf_charge_report(original_total, reason):
            return {
                "applied": False,
                "original_total": float(original_total),
                "reason": reason,
                "snapshot": [],
            }

        def _parse_psf_natom_line(line, record_index, atom_index_corrector):
            parts = line.split()
            psf_atom_index = int(parts[0])
            return {
                "record_index": record_index,
                "psf_atom_index": psf_atom_index,
                "xdata_atom_index": atom_index_corrector[psf_atom_index - 1],
                "residue_index": int(parts[2]),
                "res_name": parts[3],
                "atom_name": parts[4],
                "atom_type": parts[5],
                "charge": float(parts[6]),
            }

        def _replace_psf_charge_field(line, old_charge, new_charge):
            token = f"{old_charge:.6f}"
            token_start = line.find(token, 45)
            if token_start < 0:
                token = str(old_charge)
                token_start = line.find(token, 45)
            if token_start < 0:
                raise ValueError(
                    f"Could not locate PSF charge field for atom line: {line.rstrip()}"
                )
            charge_start = token_start - (14 - len(token))
            charge_end = token_start + len(token)
            return line[:charge_start] + f"{new_charge:14.6f}" + line[charge_end:]

        def _apply_psf_charge_neutrality_patch(psf_file, atom_index_corrector):

            if not self.json_input.get("settings", {}).get("adjust-ion-partial-charge", False):
                return

            with open(psf_file, "r") as f:
                psf_lines = f.readlines()

            natom_header_index = next(
                (i for i, line in enumerate(psf_lines) if "!NATOM" in line),
                None,
            )
            if natom_header_index is None:
                return

            natom_count = int(psf_lines[natom_header_index].split()[0])
            natom_start = natom_header_index + 1
            natom_end = natom_start + natom_count
            records = [
                _parse_psf_natom_line(line, record_index, atom_index_corrector)
                for record_index, line in enumerate(psf_lines[natom_start:natom_end])
            ]

            raw_total = sum(record["charge"] for record in records)
            original_total = self.force_field_manager._round_charge(raw_total)
            if abs(original_total) <= 0.01:
                return

            target_record = None
            target_new_charge = None
            target_adjusted_total = None

            if self.config.drude_polarisable:
                ion_units = []
                for record_index, record in enumerate(records[:-1]):
                    if record["res_name"] != "IONS":
                        continue
                    if not _is_ion_core(record["atom_type"]):
                        continue

                    particle = records[record_index + 1]
                    if particle["res_name"] != "IONS":
                        continue
                    if not _is_drude_particle(particle["atom_type"]):
                        continue

                    ion_units.append({
                        "core_index": record_index,
                        "particle_index": record_index + 1,
                        "net_charge": record["charge"] + particle["charge"],
                    })

                for ion_unit in reversed(ion_units):
                    ion_charge = self.force_field_manager._round_charge(
                        ion_unit["net_charge"]
                    )
                    if not _same_sign(ion_charge, original_total):
                        continue

                    target_record = records[ion_unit["core_index"]]
                    old_charge = target_record["charge"]
                    target_new_charge = old_charge - original_total
                    target_adjusted_total = self.force_field_manager._round_charge(
                        raw_total - old_charge + target_new_charge
                    )
                    break

            else:
                for record_index in range(len(records) - 1, -1, -1):
                    record = records[record_index]
                    atom_type = record["atom_type"]
                    if record["res_name"] != "IONS":
                        continue
                    if not _is_ion_core(atom_type):
                        continue

                    old_charge = record["charge"]
                    if not _same_sign(old_charge, original_total):
                        continue

                    new_charge = old_charge - original_total
                    adjusted_total = self.force_field_manager._round_charge(
                        raw_total - old_charge + new_charge
                    )
                    if abs(new_charge) > abs(old_charge) + 1e-9:
                        continue
                    if abs(adjusted_total) >= abs(original_total):
                        continue

                    target_record = record
                    target_new_charge = new_charge
                    target_adjusted_total = adjusted_total
                    break

            if target_record is not None:
                old_charge = target_record["charge"]
                psf_line_index = natom_start + target_record["record_index"]
                psf_lines[psf_line_index] = _replace_psf_charge_field(
                    psf_lines[psf_line_index],
                    old_charge,
                    target_new_charge,
                )
                self._psf_charge_neutrality_report = _make_psf_charge_report(
                    target_record, old_charge, target_new_charge,
                    original_total, target_adjusted_total
                )
                snapshot_start = max(natom_start, psf_line_index - 3)
                snapshot_end = min(natom_end, psf_line_index + 4)
                self._psf_charge_neutrality_report["snapshot"] = [
                    line.rstrip("\n") for line in psf_lines[snapshot_start:snapshot_end]
                ]

                with open(psf_file, "w") as f:
                    f.writelines(psf_lines)

                system_type = "Drude ion core" if self.config.drude_polarisable else "ion core"
                my_log_file.info(
                    f"PSF charge neutrality adjusted {system_type} "
                    f"{target_record['atom_type']} at PSF atom "
                    f"{target_record['psf_atom_index']}: "
                    f"{old_charge:+.6f} -> {target_new_charge:+.6f}."
                )
                return

            self._psf_charge_neutrality_report = _make_skipped_psf_charge_report(
                original_total,
                "no eligible ion core could be adjusted without violating the charge rules",
            )
            my_log_file.warning(
                "PSF charge neutrality was requested, but no eligible ion core "
                "could be adjusted without violating the charge rules."
            )

        with open(psf_file, 'w') as f:
            
            if self.config.drude_polarisable == False:
                f.write("PSF EXT CMAP\n\n")
            else:
                f.write("PSF EXT CMAP DRUDE\n\n")
            f.write("     0 !NTITLE\n\n")
            
            #!NATOM lines
            lines = self.file_manager.find_and_extract_lines(self.config.main_file, "Atoms", header=False)
            f.write(f"{len(lines):8d} !NATOM\n")
            
            tube_species_counts = defaultdict(int)
            atom_count = 1
            previous_core_props = None

            for line in lines:
                parts = self.file_manager.line_to_array(line)
                atom_index, residue_index, atom_number_type, charge, atom_type = int(parts[0]), int(parts[1]), int(parts[2]), float(parts[3]), str(parts[-1])
                x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
                mass = self.force_field_manager.get_atom_mass(atom_number_type)
                original_atom_type = atom_type

                atom_index_corrector.append(atom_index)
                atom_type_parts = atom_type.split('-')
                is_drude_particle = (
                    (len(atom_type_parts) > 1 and atom_type_parts[1] == "DRUD")
                    or atom_type in self.config.dion_particles
                )

                if len(atom_type_parts) > 1 and atom_type_parts[1] != "DRUD":
                    element = self.force_field_manager.get_atom_element(atom_type)
                    alpha, thole = self.force_field_manager._get_drude_coeffs(
                        drude_type=atom_type, atom_number_type=atom_number_type
                    )
                    res_name = atom_type_parts[1]
                    atom_type = atom_type_parts[0]

                elif len(atom_type_parts) > 1 and atom_type_parts[1] == "DRUD":
                    element = self.force_field_manager.get_atom_element(atom_type)
                    alpha, thole = self.force_field_manager._get_drude_coeffs(
                        drude_type=atom_type, atom_number_type=atom_number_type
                    )
                    res_name = "IONS" if atom_type_parts[0] in self.config.dions else "TUBE"

                else:
                    element = self.force_field_manager.get_atom_element(atom_type)
                    alpha, thole = self.force_field_manager._get_drude_coeffs(
                        drude_type=atom_type, atom_number_type=atom_number_type
                    )
                    if (atom_type in self.config.ions
                            or atom_type in self.config.dions
                            or atom_type in self.config.dion_particles):
                        res_name = "IONS"
                    else:
                        res_name = "TUBE"

                # Sheet residues carry their own name; their Drudes follow the core.
                sheet_name = self.config.residue_name_by_index.get(residue_index)
                if sheet_name is not None:
                    res_name = sheet_name

                if is_drude_particle and previous_core_props is not None:
                    residue_index, res_name = previous_core_props

                atom_name = element
                species_entry = self.config.default_species_by_type.get(original_atom_type)
                if (species_entry is not None
                        and not self.config._default_species_entry_is_drude_particle(species_entry)):
                    count = tube_species_counts[original_atom_type]
                    atom_name = species_entry["atom_element"] + str(self._convert_to_base36(count))
                    tube_species_counts[original_atom_type] += 1

                if sheet_name is not None:
                    seg_id = sheet_name
                else:
                    seg_id = "TUBE" if res_name == "TUBE" else "SYS"

                if self.config.drude_polarisable == False:
                    line_write = (
                        "{:10d} {:<8s} {:8d} {:<8s} {:<8s} {:<6s}{:14.6f}{:14.4f}           0\n".format(
                            atom_count, seg_id, residue_index, res_name, atom_name,
                            atom_type, charge, mass
                        )
                    )
                else:
                    if element == "LP":
                        val = -1
                    else:
                        val = 0

                    if is_drude_particle:
                        atom_type = "DRUD"
                        atom_name = "DP"
                    # Drude ion cores keep their real atom type.

                    line_write = (
                        "{:10d} {:<8s} {:8d} {:<8s} {:<8s} {:<6s}{:14.6f}{:14.4f}          {:2d}  {:12.3f} {:12.3f}\n".format(
                            atom_count, seg_id, residue_index, res_name, atom_name,
                            atom_type, charge, mass, val, alpha, thole
                        )
                    )

                f.write(line_write)
                atom_count += 1
                if not is_drude_particle:
                    previous_core_props = (residue_index, res_name)

            # Map atom ids once instead of scanning the correction list every time.
            atom_index_map = {orig: new + 1 for new, orig in enumerate(atom_index_corrector)}

            #!NBONDS lines
            BOND_INT_WIDTH = 10
            lines = (
                self.file_manager.find_and_extract_lines(self.config.main_file, "Bonds", header=False)
                if self.json_input["settings"]["bonds"] else []
            )
            f.write(f"\n{len(lines):10d} !NBOND: bonds\n")
            for i in range(len(lines))[::4]:
                parts1 = self.file_manager.line_to_array(lines[i]) if i < len(lines) else None
                parts2 = self.file_manager.line_to_array(lines[i+1]) if i+1 < len(lines) else None
                parts3 = self.file_manager.line_to_array(lines[i+2]) if i+2 < len(lines) else None
                parts4 = self.file_manager.line_to_array(lines[i+3]) if i+3 < len(lines) else None
                parts = [parts1, parts2, parts3, parts4]
                line = []
                for part in parts:
                    if part is not None:
                        atom1 = int(part[2])
                        atom2 = int(part[3])

                        atom1_corrected = atom_index_map[atom1]
                        atom2_corrected = atom_index_map[atom2]
                        
                        line.append(atom1_corrected)
                        line.append(atom2_corrected)

                line_write = "".join(f"{v:{BOND_INT_WIDTH}d}" for v in line) + "\n"
                f.write(line_write)

            #!NANGLES lines 
            ANGLE_INT_WIDTH = 10
            lines = (
                self.file_manager.find_and_extract_lines(self.config.main_file, "Angles", header=False)
                if self.json_input["settings"]["angles"] else []
            )
            f.write(f"\n{len(lines):10d} !NTHETA: angles\n")
            for i in range(len(lines))[::3]:
                parts1 = self.file_manager.line_to_array(lines[i]) if i < len(lines) else None
                parts2 = self.file_manager.line_to_array(lines[i+1]) if i+1 < len(lines) else None
                parts3 = self.file_manager.line_to_array(lines[i+2]) if i+2 < len(lines) else None
                parts = [parts1, parts2, parts3]
                line = []
                for part in parts:
                    if part is not None:
                        atom1 = int(part[2])
                        atom2 = int(part[3])
                        atom3 = int(part[4])
                        atom1_corrected = atom_index_map[atom1]
                        atom2_corrected = atom_index_map[atom2]
                        atom3_corrected = atom_index_map[atom3]
                        line.append(atom1_corrected)
                        line.append(atom2_corrected)
                        line.append(atom3_corrected)

                line_write = "".join(f"{v:{ANGLE_INT_WIDTH}d}" for v in line) + "\n"
                f.write(line_write)

            #!NDIHEDRALS lines
            DIHED_INT_WIDTH = 10
            lines = (
                self.file_manager.find_and_extract_lines(self.config.main_file, "Dihedrals", header=False)
                if self.json_input["settings"]["dihedrals"] else []
            )
            # A multi-term dihedral (a Fourier series) is one Dihedrals entry per term, as
            # LAMMPS needs. A PSF holds only atom indices and NAMD/CHARMM apply every
            # parameter line of the type quartet to each entry, so each quartet goes in once.
            quartets = []
            seen_quartets = set()
            for raw_line in lines:
                part = self.file_manager.line_to_array(raw_line)
                quartet = tuple(atom_index_map[int(part[k])] for k in (2, 3, 4, 5))
                key = min(quartet, quartet[::-1])
                if key in seen_quartets:
                    continue
                seen_quartets.add(key)
                quartets.append(quartet)

            f.write(f"\n{len(quartets):10d} !NPHI: dihedrals\n")
            for i in range(0, len(quartets), 2):
                line = [atom for quartet in quartets[i:i + 2] for atom in quartet]
                line_write = "".join(f"{v:{DIHED_INT_WIDTH}d}" for v in line) + "\n"
                f.write(line_write)

            #!NIMPROPERS lines
            IMPROPER_INT_WIDTH = 10 
            lines = (
                self.file_manager.find_and_extract_lines(self.config.main_file, "Impropers", header=False)
                if self.json_input["settings"]["impropers"] else []
            )
            f.write(f"\n{len(lines):10d} !NIMPHI: impropers\n")
            for i in range(len(lines))[::2]:
                parts1 = self.file_manager.line_to_array(lines[i]) if i < len(lines) else None
                parts2 = self.file_manager.line_to_array(lines[i+1]) if i+1 < len(lines) else None
                parts = [parts1, parts2]
                line = []
                for part in parts:
                    if part is not None:
                        atom1 = int(part[2])
                        atom2 = int(part[3])
                        atom3 = int(part[4])
                        atom4 = int(part[5])
                        atom1_corrected = atom_index_map[atom1]
                        atom2_corrected = atom_index_map[atom2]
                        atom3_corrected = atom_index_map[atom3]
                        atom4_corrected = atom_index_map[atom4]
                        line.append(atom1_corrected)
                        line.append(atom2_corrected)
                        line.append(atom3_corrected)
                        line.append(atom4_corrected)
                line_write = "".join(f"{v:{IMPROPER_INT_WIDTH}d}" for v in line) + "\n"
                f.write(line_write)

            _ = 0

            f.write(f"\n{_:10d} !NDON: donors\n")

            f.write(f"\n{_:10d} !NACC: acceptors\n")

            f.write(f"\n{_:10d} !NNB\n")

            #!NUMLP NUMLPH lines
            lines = (
                self.file_manager.find_and_extract_lines(self.config.main_file, "Lone Pairs", header=False)
                if self.json_input["settings"]["lone-pairs"] else []
            )
            n_lp    = len(lines)
            n_lph   = 4 * n_lp
            f.write(f"\n{n_lp:10d}{n_lph:10d} !NUMLP NUMLPH\n")
            for i in range(len(lines)):
                parts = self.file_manager.line_to_array(lines[i])
                lp_idx = i*4 + 1 # each lone pair has 4 anisotropic atoms associated with it
                line_write = f"{'3':>10s}{lp_idx:10d}{parts[2]:>3s}{parts[3]:>15s}{parts[4]:>15s}{parts[5]:>15s}" + "\n"
                f.write(line_write) 

            LONEPAIR_INT_WIDTH = 10
            for i in range(len(lines))[::2]:
                parts1 = self.file_manager.line_to_array(lines[i]) if i < len(lines) else None
                parts2 = self.file_manager.line_to_array(lines[i+1]) if i+1 < len(lines) else None
                parts = [parts1, parts2]
                line = []
                for part in parts:
                    if part is not None:
                        dependecies = part[-1].strip("()").split(',')
                        line.append(atom_index_map[int(part[1])]) # lone pair index
                        for dep in dependecies:
                            atom_dep = int(dep)
                            atom_dep_corrected = atom_index_map[atom_dep]
                            line.append(atom_dep_corrected)
                line_write = "".join(f"{v:{LONEPAIR_INT_WIDTH}d}" for v in line) + "\n"
                f.write(line_write)

            #!NUMANISO lines
            if self.config.drude_polarisable: # bricks for non-polarisable sims, anisotropy for polarisable sims
                if self.json_input["settings"]["anisotropy"]:
                    lines = self.file_manager.find_and_extract_lines(
                        self.config.main_file, "Anisotropy", header=False
                    )
                    bond_coefs = self.file_manager.find_and_extract_lines(
                        self.config.main_file, "Bond Coeffs", header=False
                    )
                else:
                    lines = []
                    bond_coefs = []
                f.write(f"\n{len(lines):10d} !NUMANISO\n")
                for i in range(len(lines)):
                    parts = self.file_manager.line_to_array(lines[i])
                    anisotropy_type = parts[-2]
                    drude_type = f"{anisotropy_type.split('-', 1)[0]}-DRUD"
                    Kdrude = None
                    for bond in bond_coefs:
                        bond_parts = self.file_manager.line_to_array(bond)
                        if {anisotropy_type, drude_type} == set(bond_parts[-2:]):
                            Kdrude = float(bond_parts[1])
                            break
                    if Kdrude is None:
                        # parts[-1] is the dependency tuple string; the atom type is parts[-2].
                        my_log_file.error(f"Could not find Drude bond coefficient for anisotropic atom type {parts[-2]}.")
                        raise ValueError(f"Could not find Drude bond coefficient for anisotropic atom type {parts[-2]}.")
                    
                    K11, K22, K33 = self.force_field_manager._get_anisotropic_coeffs_CHARMM36_NDP(A11=parts[2], A22=parts[3], Kdrude=Kdrude)
                    f.write(f"{K11:>20f}{K22:>15f}{K33:>15f}\n")
                
                ANISO_INT_WIDTH = 10
                for i in range(len(lines))[::2]:
                    parts1 = self.file_manager.line_to_array(lines[i]) if i < len(lines) else None
                    parts2 = self.file_manager.line_to_array(lines[i+1]) if i+1 < len(lines) else None
                    parts = [parts1, parts2]
                    line = []
                    for part in parts:
                        if part is not None:
                            dependecies = part[-1].strip("()").split(',')
                            line.append(atom_index_map[int(part[1])])  # atom index
                            for dep in dependecies:
                                atom_dep = int(dep)
                                atom_dep_corrected = atom_index_map[atom_dep]
                                line.append(atom_dep_corrected)
                    line_write = "".join(f"{v:{ANISO_INT_WIDTH}d}" for v in line) + "\n"
                    f.write(line_write)

            f.write(f"\n{_:10d} !NCRTERM: cross-terms\n")

            f.write(f"\n{_:10d} !AUTOGEN\n") # useless, tell CHARMM which angles etc were guessed

        _apply_psf_charge_neutrality_patch(psf_file, atom_index_corrector)

        my_log_file.info(f"File type .psf file written to {psf_file}")

    def write_pdb_restraints_file(self):

        '''
        Generate a positional-restraint reference PDB for the Drude simulation.
        Nanotube carbons receive occupancy and B-factor 1.00 to activate the
        harmonic restraint; all others receive 0.00. Skips if no PDB is found.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        spring_constant = 1.00 # how strongly the harmonic potential acts
        B_factor = 1.00 # 1.00 for yes 0.00 for no

        src = self.config.main_file
        pdb_file = src[:-6] + ".pdb"

        if os.path.exists(pdb_file) == False:
            my_log_file.error("PDB file not found for generating restraints.")
            return
        
        restraints_file = os.path.join(self.config.system_folder, f"restraints.ref")

        # Tube-end sheets are held with the tube; box-end sheets stay free.
        held_sheets = set()
        if self._tube_sheet_report is not None:
            held_sheets = set(self._tube_sheet_report['resnames'].values())
        
        with open(pdb_file, 'r') as f:
            with open(restraints_file, 'w') as rf:
                lines = f.readlines()
                previous_core_line = None
                for line in lines:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        atom_name = line[12:16].strip()
                        if atom_name == "DP" and previous_core_line is not None:
                            atom_index = int(previous_core_line[6:11]) + 1
                            line = (
                                line[:6] + f"{atom_index:5d}" + line[11:17]
                                + previous_core_line[17:26] + line[26:66]
                                + previous_core_line[66:]
                            )
                        res_name = line[17:21].strip()  # residue name (columns 18-20)

                        # it is a carbon in the nanotube, but not the functional group.
                        if atom_name != "DP" and res_name == "TUBE" and "C" in line[12:16].strip() and line[12:16].strip() != "C":

                            # Format occupancy (cols 55-150) and B-factor (cols 61-66)
                            new_occ = f"{B_factor:6.2f}"
                            new_b = f"{spring_constant:6.2f}"
                            line = line[:54] + new_occ + new_b + line[66:]
                        elif atom_name != "DP" and res_name in held_sheets:
                            new_occ = f"{B_factor:6.2f}"
                            new_b = f"{spring_constant:6.2f}"
                            line = line[:54] + new_occ + new_b + line[66:]
                        else:
                            new_occ = f"{0.00:6.2f}" # No restraint
                            new_b = f"{0.00:6.2f}"
                            line = line[:54] + new_occ + new_b + line[66:]
                        if atom_name != "DP":
                            previous_core_line = line
                    if "END" in line or "REMARK" in line:
                        continue
                    rf.write(line)

    # -------------------------------------------------------------------------
    # MD FORCE FIELD PARAMETER FILES
    # -------------------------------------------------------------------------

    def copy_other_FF_files(self):

        '''
        Copy all other files apart from the _solvents and _fragments files.

        Parameters:
        ----------
        None

        Returns:
        -------
        None

        Raises:
        ------
        None
        '''

        ff_low = self.config.force_field.lower()
        base_folder = os.path.join(self.config.force_fields_root, ff_low)
        dest = os.path.join(self.config.system_folder, os.path.basename(base_folder.rstrip("/")))

        ignore = shutil.ignore_patterns(
            "chirality_kit_solvents.data",
            "chirality_kit_fragments.data"
        )
        shutil.copytree(base_folder, dest, ignore=ignore, dirs_exist_ok=True)

        my_log_file.info(f"{self.config.force_field} other force field files copied to system folder.")
