#!/bin/bash

tar -xzf static_database.tar.gz
chmod a+x ddg_monomer.static.linuxgccrelease

./ddg_monomer.static.linuxgccrelease @ddg_flags -ddg::mut_file pos_number_identity.mutfile > RUN_4A0T_pos_number_identity.log

tar -czf pos_number_identity_pdbs.tar.gz *.pdb

rm -rf database/
rm static_database.tar.gz ddg_monomer.static.linuxgccrelease *.pdb 

