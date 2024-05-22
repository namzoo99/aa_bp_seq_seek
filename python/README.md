### About ```breakpoints_to_bed.py``` used in ```aa_bp_seq_seek```
#### I only remove ```regionstring``` in ```AmpliconSuite-pipeline/scripts/breakpoints_to_bed.py``` line 104.
#### Original ```breakpoints_to_bed.py``` used ```regionstring = aa interval region``` as name of output file, and if aa interval region is too long for the file name, error occurs.
