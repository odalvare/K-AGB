![sponsors](/auxiliar_git/logos/logosKAGB.jpg)

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

## Conceptualization

![K-AGB conceptualization](/auxiliar_git/images/ConceptualFrameWork.jpg)

Figure extracted from Álvarez-Villa et al. (2023), rights protected.

# Required packages

K-AGB requires the installation of the following packages for its correct performance. It is strongly encoraged the use of anaconda to create an individual environmente (https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

### Basic packages

- python v. 3.7.12.
- matplotlib v. 3.5.7.
- pandas v. 1.3.4.

### Computing packages

- numpy v. 1.21.5.
- numba v. 0.55.1.
- scipy v. 1.7.3.
- scikit-learn v. 1.0.2.
- statsmodels v. 0.13.1.
- fsspec v. 2022.11.10.

### GIS packages

- gdal v. 3.0.2.
- fiona v. 1.8.13.
- rasterio v. 1.1.0.
- geopandas v. 0.10.2.
- shapely v. 1.7.1.
- pyproj v. 6.2.1.

### Geostatistics packages

- geostatspy v. 0.0.26.

# Data availability

The information required to perform a test biomass simulation via K-AGB can be downloaded from: https://nextcloud.emergente.com.co/nextcloud/index.php/s/n64iyAoLfnYkqAf

# Support

Long-term support of the K-AGB utilities will be provided by Emergente (https://www.emergente.com.co/). Dr Oscar D. Álvarez-Villa (email: oscar.alvarez@emergente.com.com) is the principal developer of the project.
