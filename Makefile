UTILS_PARTIAL = partials/utils.R
export UTILS_PARTIAL
DISPLAY_PARTIAL = partials/display.R
export DISPLAY_PARTIAL
UTILS_CPP_PARTIAL = partials/utils.cpp
export UTILS_CPP_PARTIAL
HMC_EXACT_CPP_PARTIAL = partials/hmc-exact.cpp
export HMC_EXACT_CPP_PARTIAL

INVERSION_START_DATE = 2014-09-01
INVERSION_END_DATE = 2021-04-01

GEOS_CHEM_SOURCE = external/GEOS_Chem/Code.v12.3.2
GEOS_CHEM_EXTDATA = data/GEOS_Chem/ExtData

AREA_1X1 = data/area-1x1.nc
OCO2_SIF_DIRECTORY = data/OCO2_L2_Lite_SIF_10r
OCO2_OBSERVATIONS = data/OCO2_b10c_10sec_GOOD_r5.nc4
OBSPACK_DIRECTORY = data/obspack_co2_1_OCO2MIP_v3.2_2021-05-20/data/daily
TCCON_SOUNDING_DIRECTORY = data/downloaded_20211217
LAUDER_DATA = data/lauder_co2_2014_2021.50_ooofti_lhr.csv
SIB4_HOURLY_DIRECTORY = data/sib4-hourly
LANDSCHUTZER_INVENTORY_INPUT = data/spco2_MPI-SOM_FFN_v2020.nc
TRANSCOM_MASK_ORIGINAL_1X1 = data/TRANSCOM_mask_original_1x1.nc
GEOS_CHEM_RESTART_FILE = data/GEOSChem.Restart.20140901_0000z.nc4
GEOS_2X25_GRID = data/geos.2x25.grid
FOSSIL_INVENTORY_PATH = data/fossil-mipv10
FIRES_INVENTORY_PATH = data/GFED4.1s_Aleya_hemco_080321
XBASE_DIRECTORY = data/FLUXCOM_XBase
XBASE_05X05_GRID = data/xbase.05x05.grid
SAMPLES_WOMBAT_V2 = data/samples-LNLGIS.rds

all:
	echo "Please refer to the README for instructions on how to run this project"

include 1_transport/module.mk
include 2_matching/module.mk
include 3_sif/module.mk
include 4_inversion/module.mk
include 5_results/module.mk

1_transport_targets: $(1_TRANSPORT_TARGETS)
2_matching_targets: $(2_MATCHING_TARGETS)
3_sif_targets: $(3_SIF_TARGETS)
4_inversion_targets: $(4_INVERSION_TARGETS)
5_results_targets: $(5_RESULTS_TARGETS)

.SECONDARY: $(SECONDARY_TARGETS)
