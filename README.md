# K-AGB
Stochastic simulation of aboveground biomass driven by land cover changes.

## Summary
Tropical Montane Forests (TMFs) are known for their freshwater resources, which feed rivers passing through major cities in the Andes. TMFs are highly biodiverse and essential for capturing carbon. However, measuring carbon pools in these topographically complex ecosystems is challenging. K-AGB implements a novel method for quantifying the carbon stored in the Aboveground Biomass (AGB) that combines several data sources, quantifies the predictive uncertainty, and uses land cover as the predictive variable.

**Reference:** Álvarez-Villa et al. (2023) Spatiotemporal dynamics of above-ground biomass in a high tropical montane basin. Submitted to Environmental Modeling and Software.

# Features

K-AGB is a Python framework that enables the use of stochastic simulation of aboveground biomass (stored by the vegetation) and the equivalente of stored carbon. Some features that you can use from K-AGB are:

- Stochastic geostatistical simulation supported by geostatspy (https://github.com/GeostatsGuy/GeostatsPy) to fill uncertainty models of the secondary biomass data (e.g. GEDI database, from: https://daac.ornl.gov/GEDI/guides/GEDI_L4B_Gridded_Biomass.html).
- Statistic downscaling using multivariate linear regression with NDVI, EVI and LAI (use MODIS database as a reference: https://modis.gsfc.nasa.gov/).
- Estimation of biomass probability distrbutions and quantile transformation using prumery data.
- MonteCarlo simulation of biomass using land cover category as the main driver.
