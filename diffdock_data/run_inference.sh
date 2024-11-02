cd ../DiffDock

inference_args="../diffdock_data/inference_args.yaml"
protein_ligand_csv="../diffdock_data/diffdock_input.csv"
out_dir="../diffdock_data/results/"

python -m inference --config $inference_args  --protein_ligand_csv $protein_ligand_csv --out_dir $out_dir

cd ../diffdock_data

python -m postprocess.py 