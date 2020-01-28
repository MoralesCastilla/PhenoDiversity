# Forecasting Phenology and Climatic Suitability of Crop Cultivar Diversity 

Code and data sources utilized in Morales-Castilla et al. (2020)[https://www.pnas.org/content/early/2020/01/21/1906731117](url) to forecast phenology and climatic suitability for winegrape cultivars over the 1950-2100 period.


## Workflow

The modelling workflow is as follows:

1. Parameterization of phenological models (done externally through PMP V5.0 (Chuine et al. 2013).
2. Climate data extraction, bias-correction and storing.
3. Phenological forecasting (makes predictions of dates of budbreak, flowering and veraison, for a given year and variety). 
4. Determination of the climatic envelope for maturation during a baseline period (1950-1980).
5. Forecasting of the climatic envelope for the maturation process, for a given variety and year.
