#!/bin/bash

#Define paths to Rosetta and the database
#---------------------------------------------------------------------------------------------------------
rosetta=/scratch/ameger/software/Rosetta/main/source/bin
database=/scratch/ameger/software/Rosetta/main/database
cd ../setup/
#---------------------------------------------------------------------------------------------------------

#Make 4A0T_cleaned.pdb by manually removing everything in 4A0T cyrstal structure except residues 463-553 of chains A/B/C
#Idealize the crystal structure to remove any unusual bonds, angels and dihedrals
#---------------------------------------------------------------------------------------------------------
${rosetta}/idealize_jd2.linuxgccrelease -database ${database} -in::file::fullatom -s 4A0T_cleaned.pdb -no_optH false -flip_HNQ > idealized.log
mv 4A0T_cleaned_0001.pdb 4A0T_idealized.pdb
#---------------------------------------------------------------------------------------------------------

#Renumber the pdb
#---------------------------------------------------------------------------------------------------------
${rosetta}/score_jd2.default.linuxgccrelease -renumber_pdb -ignore_unrecognized_res -s 4A0T_idealized.pdb -out:pdb > renumbered.log
mv 4A0T_idealized_0001.pdb 4A0T_renumbered.pdb
#---------------------------------------------------------------------------------------------------------

#Use Relax to optimize the starting structure
#---------------------------------------------------------------------------------------------------------
${rosetta}/relax.linuxgccrelease -database ${database} -relax::sequence_file always_constrained_relax_script -constrain_relax_to_start_coords -relax::coord_cst_width 0.25 -relax::coord_cst_stdev 0.25 -s 4A0T_renumbered.pdb -in::file::fullatom -no_optH false -flip_HNQ > relaxed.log
mv 4A0T_renumbered_0001.pdb 4A0T_relaxed.pdb
#---------------------------------------------------------------------------------------------------------

#Perform an energy minimization and generate a coordinate constraint file
#---------------------------------------------------------------------------------------------------------
${rosetta}/minimize_with_cst.linuxgccrelease -database ${database} -in:file:s 4A0T_relaxed.pdb -in:file:fullatom -ignore_unrecognized_res -fa_max_dis 9.0 -ddg::harmonic_ca_tether 0.5 -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix min_cst_0.5 -ddg::sc_min_only false -score:patch ${database}/scoring/weights/score12.wts_patch > mincst.log
./convert_to_cst_file.sh mincst.log > input.cst
mv min_cst_0.5.4A0T_relaxed_0001.pdb 4A0T.pdb
#---------------------------------------------------------------------------------------------------------


