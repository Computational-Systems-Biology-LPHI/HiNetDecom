# Nested hierarchical decomposition for two real-life biomodels b2AR and EGFR.

This GitHub repository contains a Python implementation of the nested hierarchical decomposition for the models of B2AR and EGFR. This package corresponds to the journal article "Hierarchical Optimization of Biochemical Networks", by Nisha Ann Viswan, Alexandre Tribut, Manvel Gasparyan, Ovidiu Radulescu, and Upinder S. Bhalla,
which is publicly available at . For more details, please refer to that paper.

The package includes six SBML files corresponding to two biomodels: b2AR_PKA_v5 and D4_model_EGFR_v13b. Each model has three associated files:

1. b2AR_PKA_v5_original.xml - This SBML file represents the original model, referred to as Model1.
2. b2AR_PKA_v5_removed_reaction.xml - This SBML file represents the original Model1 with certain reversible reactions modified
   to be irreversible (refer to the cited journal article for details). We further refer to this model as Model2.
4. b2AR_PKA_v5_reduced.xml - This SBML file represents the Michaelis-Menten reduced version of Model2.
   
Similarly, the files for EGFR are:

1. D4_model_EGFR_v13b_original.xml
2. D4_model_EGFR_v13b_removed_reaction.xml
3. D4_model_EGFR_v13b_reduced.xml

The Python code accepts an SBML file (either 'b2AR_PKA_v5_reduced.xml' or 'D4_model_EGFR_v13b_reduced.xml') as input and produces hierarchical levels of species subsets. Each subset is nested within the next, collectively forming the complete decomposition.

To run the code:
1. Ensure Python3 is Installed:
Make sure Python 3.x is installed on your system. You can check this by running:
python3 --version
If Python 3 is not installed, download and install it from the official Python website.
2. Install Required Packages:
The script uses several Python packages. Install them using pip:
pip3 install libsbml numpy networkx matplotlib
3. Download the Script:
Save the decomposition_procedure.py script to your working directory.
4. Prepare the Input File:
Ensure the D4_model_EGFR_v13b_removed_reaction.xml file (or another specified SBML file) is in the same directory as the script.
5. Run the Script:
Execute the script from the command line. Navigate to the directory containing decomposition_procedure.py and run:
python3 decomposition_procedure.py


