
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Data

Below is information on the data in this subfolder, as well as
instructions for obtaining from the [Triticeae Toolbox
(T3)](https://triticeaetoolbox.org/barley) the experimental data used in
this study. These instructions valid as of 27 March, 2019.

## Data in this subfolder

1.  `project_entries.csv` - line name, program name, population group,
    predigree, family, and notes for all barley lines used in this
    study.
2.  `trial_metadata.csv` - information on the experimental trials used
    in this study.
3.  `PVV_BLUE.RData` - a binary file containing the adjusted phenotypic
    data on the training population and validation families.

## Data from T3

### Training population and parent candidates

#### Genomewide marker data

Currently, the genomewide marker data is not publicly available on T3.
Please [email me](neyha001@umn.edu) to obtain this data.

#### Phenotype data

1.  Go to <https://triticeaetoolbox.org/barley>.
2.  Under the “Select” tab, go to “Wizard (Lines, Traits, Trials)”.
3.  Fill in the following information in the selection prompts:
      - Breeding Program: University of Minnesota (MN)
      - Year: 2014 and 2015
      - Trials: S2TP\_2014\_Crookston, S2TP\_2014\_FHB\_Crookston,
        S2TP\_2014\_FHB\_StPaul, S2TP\_2015\_Crookston,
        S2TP\_2015\_Crookston\_FHB, S2TP\_2015\_StPaul, and
        S2TP\_2015\_StPaul\_FHB
      - Traits: plant height, heading date, and FHB Severity
4.  Click the “Save current selection” button.
5.  Under the “Download” tab, go to “Genotype and Phenotype Data”.
6.  Make sure only the “Phenotype” box is checked.
7.  Click the “Create file” button with instructions for one column for
    each trait (not used by TASSEL).
8.  Click the “Download Zip file of results” button.

### Validation population

#### Phenotype data

1.  Go to <https://triticeaetoolbox.org/barley>.
2.  Under the “Select” tab, go to “Wizard (Lines, Traits, Trials)”.
3.  Fill in the following information in the selection prompts:
      - Breeding Program: University of Minnesota (MN)
      - Year: 2017 and 2018
      - Trials: PVV\_2017\_Crookston, PVV\_2017\_StPaul,
        PVV\_FHB\_2017\_Crookston, PVV\_FHB\_2017\_StPaul,
        PVV\_FHB\_2018\_Crookston, and PVV\_FHB\_2018\_StPaul
      - Traits: plant height, heading date, and FHB Severity
4.  Click the “Save current selection” button.
5.  Under the “Download” tab, go to “Genotype and Phenotype Data”.
6.  Make sure only the “Phenotype” box is checked.
7.  Click the “Create file” button with instructions for one column for
    each trait (not used by TASSEL).
8.  Click the “Download Zip file of results” button.
