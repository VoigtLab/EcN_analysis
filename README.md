# EcN_analysis
Files accompanying Triassi et al., Redesign of an _Escherichia coli_ Nissle treatment for phenylketonuria using insulated genomic landing pads and genetic circuits to reduce burden, Cell Systems (2023), https://doi.org/10.1016/j.cels.2023.05.004

There are four folders:
* dna-maps_genbank
	- genbank files for plasmids and genome integrations in this paper
* rna-seq
	- All codes utilized in producing RNA-seq-related data and plots can be found in: EcN_analysis.py
	- Metadata including experimental information can be found in: 2021_01_07_EcN_hdoost_h_dataframe.csv
	- Analyzed output differential expression and p-values are shown in: 2021_04_09_DE_all_23h.csv
	- Other .csv files include QC, metrics, and raw & FPKM values for gene expression.

* plotted_data_only
	- Excel files containing values used to generate plots, no data anaylsis

* raw-data_code_and_notebooks
	- a bunch of Python notebooks used to analyze raw data and generate plots
		* notebooks are nested within subfolders
		* notebooks often contain superfluous code from prototyping or copied from previous analysis, some thinking will be required to run these
		* some raw data files were too large to upload to GitHub, but files can be shared upon request
	- Python scripts that are utilized by these notebooks
		* cytoflow_modules
			- cf_helpers.py and cf_imports.py
		* transferfunction.py
			- defines a class for creating transfer functions
		* transferfunction_beta.py
			- updated version used in later notebooks
