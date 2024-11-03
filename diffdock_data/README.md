# diffdock_data

Run diffdock with
```bash
sh run_inference.sh
```
which runs diffdock on the examples given in `diffdock_input.csv`.

To generate `diffdock_input.csv`, run 
```bash
python generate_diffdock_input_csv.py
```
which, by default, creates an input csv using the belka ids specified in `ids.txt`.
If you run
```bash
python generate_diffdock_input_csv.py {K}
```
where `K` is an integer, the script will create an input csv using K randomly selected belka ids.

The output will be in the `results` directory, with each output file (`rank1_confidence{CONFIDENCE SCORE HERE}.sdf`) being saved in a subdirectory named after the belka id.


## Folder structure required
This requires belka's test set saved in `belka/test.csv` and also requires DiffDock's repository to be cloned in `../DiffDock/` (relative to this repository).