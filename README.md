# The Impact of Air conditioning on Residential Electricity Demand across World Countries.
Replication package for Enrica De Cian, Giacomo Falchetta, Filippo Pavanello, Yasmin Romitti, Ian Sue Wing, "**The Impact of Air conditioning on Residential Electricity Demand across World Countries**" \[Journal of Environmental Economics and Management\]

[Old Working Paper](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4604871)

For questions: [pavanello@ifo.de](mailto:pavanello@ifo.de) (econometrics) and [giacomo.falchetta@cmcc.it](mailto:giacomo.falchetta@cmcc.it) (projections)

# Description
This repository provides the codes required to reproduce the tables, figures, and in-text summary statistics in De Cian et al. (2025, JEEM). 

## Folders

 - `rscript/` - contains all R codes that replicates tables and figures in the paper
 - `interm/` - where intermediate inputs are stored
 - `output/` - where all results are stored

## Data for replications
The data used in the paper are available for download in the following repository: [Link](https://doi.org/10.5281/zenodo.14990955)

## Steps

In the `rscript/` folder the analysis is conducted following **nine** steps:

1. **Descriptives** - code to produce the descriptive analysis in the paper (Table 1, Table 2, Figure 1, Figure A1)
2. **Main** - code to produce the main specification of the paper (Table 3, Figure 2, Table A1, Table A2)
3. **Heterogeneity** - code to produce the heterogeneity analysis in the paper, based on income quintiles and country-level (Figure 3, Figure 4, Table A5)
4. **Standardised** - code to produce the country-specific standardised regressions and the descriptive meta-analysis of the coefficients (Figure 5, Figure A2)
6. **Additional analysis** - code to produce auxiliary analysis with: (a) solar PV generation (Table 5, Table A11-A14), (b) other appliances (Table A6-A9)
5. **Robustness** - code to produce the robustness checks (Table A3, Table A4)
7. **Projections** - code to produce the projections and the relative analysis (Figure A4, Figure A5, Figure A6)
8. **Implications** - code to produce the 'implications' results: future AC ownership and consumption, budget share, distribution of AC electricity consumption, and CO2 emissions (Table 4, Figure 6, Figure 7, Figure A3, Table A15)
9. **Supplementary** - code to produce Supplementary Material results (Table S1-S18, Figure S1-S6)
