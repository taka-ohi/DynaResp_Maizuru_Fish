# R scripts for Ohigashi et al. (*in prep*)

This repository contains the codes to reproduce the results and figures presented in Ohigashi et al. (*in prep*).

The preprint is available here: <a href="https://doi.org/10.64898/2026.02.16.706069"><img src="https://img.shields.io/badge/Preprint-bioRxiv-red.svg" valign="middle"/></a>

------------------------------------------------------------------------

## Data

The analysis is built upon three primary data files located in the `01_dataformatting_out/` directory:

-   **`fishcount.txt`**: Time-series records of fish abundance (underwater visual census) from Maizuru Bay. Each row represents a specific observation date and species count.
-   **`metadata.txt`**: Time-series of environmental variables, including water temperature, corresponding to the census dates.
-   **`fishinfo_w_ecology.csv`**: Information for each fish species (e.g., latitudinal center).

## Computational environment

-   **R version**: 4.5.0
-   **Key dependencies**: `rEDM` (v0.9.12), `rUIC` (v0.7.5)
