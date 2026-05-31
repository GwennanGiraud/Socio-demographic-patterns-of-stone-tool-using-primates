### Citation

Please cite our paper if you use these data, code, or results in your own work.


### Repository Contents

## Getting Started

The main analytical workflow is provided in `R_Script_Analysis_Publication.R`, which contains all statistical analyses presented in the study. Any additional scripts required during the workflow are referenced directly within this file at the appropriate stage of the analysis. The descriptions below summarize the content of each file and its dependencies.


## Dataset construction (N=49 individuals)

# Attributes_Table.csv

Contains the characteristics of each individual, including:

* Individual ID code
* Individual name
* Age class
* Age group
* Sex
* Dominance rank
* Hair pattern (phenotype)
* Tool-use status (tool user = Yes; non-tool user = No)
* Number of focal observations collected
* Focal sampling effort (minutes)

Used with: `R_Script_Metrics_Calcul_Publication.R`

---

# Grooming_Dataset.csv

Contains the raw grooming interaction data for all focal individuals (infants excluded). Variables include:

* Date
* Time
* Focal ID
* Groomer ID
* Receiver ID
* Start time
* End time
* Grooming duration (seconds)

Used with: `R_Script_Metrics_Calcul_Publication.R`

---

# R_Script_Metrics_Calcul_Publication.R

R script used to:

* Create the grooming matrix (`mGroom`)
* Calculate social network analysis (SNA) grooming metrics corrected for sampling effort

Requires:

* `Attributes_Table.csv`
* `Grooming_Dataset.csv`



## Dataset analysis (N=42 individuals)

# R_Script_Analysis_Publication.R

R script containing all statistical analyses conducted in the study.

Requires:

* `Environment_Analysis_Published.RData`
* `diagnostic_fcns.R`

This script also calls several complementary scripts:

* `R_Script_Grooming_Sociogram_Publication.R`
* `R_Script_to_prepare_for_Python.R`
* `Python_Script_for_homophily_correction.py`

---

# Environment_Analysis_Published.RData

R environment containing the datasets used for the analyses.

Used with:

* `R_Script_Analysis_Publication.R`
* `R_Script_Grooming_Sociogram_Publication.R`

---

# diagnostic_fcns.R

Functions used for model diagnostics.

Used with:

* `R_Script_Analysis_Publication.R`

---

# R_Script_Grooming_Sociogram_Publication.R

R script used to generate the grooming network sociogram illustrating the Social Position Model presented in `R_Script_Analysis_Publication.R`.

Requires:

* `Environment_Analysis_Published.RData`

---

# R_Script_to_prepare_for_Python.R

R script used to configure the R environment for running the Python homophily correction analysis.

Requires:

* `src` folder

Please, run this script before `Python_Script_for_homophily_correction.py`

---

# Python_Script_for_homophily_correction.py

Python script used to calculate assortativity values corrected using the method proposed by Karimi & Oliveira (2023).

Requires:

* `src` folder

Please, run this script after `R_Script_to_prepare_for_Python.R`

The corrected assortativity results are also reported in:

* `R_Script_Analysis_Publication.R`

---

# src

Folder containing files and functions required for the assortativity correction procedure adapted from Karimi & Oliveira (2023).

Contents include:

* `df` (Dataset from `Environment_Analysis_Published.RData`)
* `mat` (`mGroom` matrix from `Environment_Analysis_Published.RData`)
* Additional helper functions

Used with:

* `R_Script_to_prepare_for_Python.R`
* `Python_Script_for_homophily_correction.py`
