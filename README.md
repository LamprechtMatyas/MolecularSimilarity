# Modeling of fragment-based molecular similarity
## Installation
### Requirements
* Python 3.5
* RDkit
* scikit-learn
* matplotlib

### Download
git clone https://github.com/LamprechtMatyas/MolecularSimilarity.git

## Usage
Program can be used in 2 ways - script by script or run_all_scripts:
1) script by script: If you want to see intermediate results you can run scripts in order:
    extract_fragments.py -> compute_descriptors.py -> create_model.py -> 
    -> score_molecules.py -> add_activity.py -> compute_evaluation.py
    At the beginning of each script is a description of the input files and other requirements.
2) run_all_scripts: If you want to run whole process, then use run_all_scripts.py.
    You will get same results as in first case after running compute_evaluation.py.
    Here you can also see some intermediate results, if you specify that you want them
    in input parameters.
    As default input files there are used the files that are stored in data section,
    but you can use your own ones.
      
In model folder there are models that can be used for molecular similarity computing.
If you want to use one, read the beginning of the script and by instruction configure
your own configuration.json file or you can use one of these connfiguration files
that are stored in data folder.

For comparing the outputs of more models you can use scripts that have got "graph"
or "graphs" in the name (e.g. print_graph.py or print_evaluation_graphs.py), that
will show you outputs of compute_evaluation.py or run_all_scripts in graphical way.
