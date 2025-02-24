5_RESULTS_SRC_DIR = 5_results/src
5_RESULTS_INTERMEDIATES_DIR = 5_results/intermediates
5_RESULTS_FIGURES_DIR = 5_results/figures
5_RESULTS_OSSE_DIR = 5_results/figures/osse_flux_decompositions
5_RESULTS_PRODUCTS_DIR = 5_results/products

$(shell mkdir -p $(5_RESULTS_INTERMEDIATES_DIR))
$(shell mkdir -p $(5_RESULTS_FIGURES_DIR))
$(shell mkdir -p $(5_RESULTS_OSSE_DIR))
$(shell mkdir -p $(5_RESULTS_PRODUCTS_DIR))


# Intermediates

PERTURBATIONS_AUGMENTED = $(5_RESULTS_INTERMEDIATES_DIR)/perturbations-augmented.fst
PERTURBATIONS_AUGMENTED_ZONAL = $(5_RESULTS_INTERMEDIATES_DIR)/perturbations-augmented-zonal.fst
OSSE_FLUX_AGGREGATES_SAMPLES_BASE = $(5_RESULTS_INTERMEDIATES_DIR)/osse-flux-aggregates-samples
OSSE_FLUX_AGGREGATES_SAMPLES_CASES = $(foreach OSSE_CASE,$(OSSE_CASES),$(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-$(OSSE_CASE).rds)

REGION_GRID = $(5_RESULTS_INTERMEDIATES_DIR)/region-grid.rds
REGION_SF = $(5_RESULTS_INTERMEDIATES_DIR)/region-sf.rds
SIX_YEAR_AVERAGE = $(5_RESULTS_INTERMEDIATES_DIR)/six-year-average.fst
LANDSCHUTZER_CLIMATOLOGY_2X25 = $(5_RESULTS_INTERMEDIATES_DIR)/spco2_MPI-SOM_FFN_v2020_climatology-2x25.nc

XBASE_FLUXES = GPP TER NEE
XBASE_YEARS = 2015 2016 2017 2018 2019 2020
XBASE_TER_BASE = $(5_RESULTS_INTERMEDIATES_DIR)/xbase_TER_unweighted
XBASE_TER_FILES = $(foreach YEAR,$(XBASE_YEARS),$(XBASE_TER_BASE)_$(YEAR).nc)
XBASE_MONTHLY_BASE = $(5_RESULTS_INTERMEDIATES_DIR)/xbase_monthly
XBASE_MONTHLY_FILES = \
	$(foreach FLUX,$(XBASE_FLUXES),\
	$(foreach YEAR,$(XBASE_YEARS),\
	$(XBASE_MONTHLY_BASE)_$(FLUX)_$(YEAR).nc))
XBASE_MONTHLY_2x25_BASE = $(5_RESULTS_INTERMEDIATES_DIR)/xbase-monthly-2x25
XBASE_MONTHLY_2x25_FILES = $(foreach FLUX,$(XBASE_FLUXES),$(XBASE_MONTHLY_2x25_BASE)-$(FLUX).nc)
XBASE_MONTHLY_2x25 = $(5_RESULTS_INTERMEDIATES_DIR)/xbase-monthly-2x25.fst
XBASE_MONTHLY_2x25_ZONAL = $(5_RESULTS_INTERMEDIATES_DIR)/xbase-monthly-2x25-zonal.fst

5_RESULTS_TARGETS += \
	$(5_RESULTS_PRODUCTS_DIR)/WOMBAT_v2S_CO2_gridded_flux_samples.nc4 \
	$(5_RESULTS_PRODUCTS_DIR)/WOMBAT_v2S_CO2_gridded_climatology_samples.nc4 \
	$(5_RESULTS_FIGURES_DIR)/region-map.pdf \
	$(5_RESULTS_FIGURES_DIR)/observation-count.pdf \
	$(5_RESULTS_FIGURES_DIR)/sif-gpp-average-slope.pdf \
	$(5_RESULTS_FIGURES_DIR)/sif-gpp-map-slope.pdf \
	$(5_RESULTS_FIGURES_DIR)/sif-gpp-map-intercept.pdf \
	$(5_RESULTS_FIGURES_DIR)/error-params-table.tex \
	$(5_RESULTS_FIGURES_DIR)/prior-flux-decomposition.pdf \
	$(5_RESULTS_FIGURES_DIR)/osse-true-fluxes.pdf \
	$(5_RESULTS_FIGURES_DIR)/osse-metrics-table-rmse.tex \
	$(5_RESULTS_FIGURES_DIR)/osse-metrics-table-crps.tex \
	$(5_RESULTS_FIGURES_DIR)/osse-flux-decomposition-main.pdf \
	$(5_RESULTS_FIGURES_DIR)/flux-annual-average-table.csv \
	$(5_RESULTS_FIGURES_DIR)/flux-global.pdf \
	$(5_RESULTS_FIGURES_DIR)/flux-decomposition.pdf \
	$(5_RESULTS_FIGURES_DIR)/flux-net-zonal.pdf \
	$(5_RESULTS_FIGURES_DIR)/seasonal-cycle-zonal.pdf \
	$(5_RESULTS_FIGURES_DIR)/seasonal-latitude-profile.pdf \
	$(5_RESULTS_FIGURES_DIR)/average-maps-main-gpp.pdf \
	$(5_RESULTS_FIGURES_DIR)/average-maps-supp-gpp.pdf \
	$(5_RESULTS_FIGURES_DIR)/average-maps-resp.pdf \
	$(5_RESULTS_FIGURES_DIR)/average-maps-nee.pdf \
	$(5_RESULTS_FIGURES_DIR)/effective-sample-size.txt \
	$(5_RESULTS_FIGURES_DIR)/traceplots.pdf

REGIONS = global Region01 Region02 Region03 Region04 Region05 Region06 Region07 Region08 Region09 Region10 Region11
RESP_TYPES = FIXRESP FREERESP
OSSE_FLUX_DECOMPOSITIONS = \
	$(foreach OSSE_CASE,$(OSSE_BASE_CASES),\
	$(foreach REGION,$(REGIONS),\
	$(foreach RESP_TYPE,$(RESP_TYPES),\
	$(5_RESULTS_OSSE_DIR)/$(OSSE_CASE)-$(REGION)-$(RESP_TYPE).pdf)))
OSSE_FLUX_DECOMPOSITIONS += \
	$(foreach OSSE_CASE,$(OSSE_BASE_CASES),\
	$(foreach REGION,$(REGIONS),\
	$(5_RESULTS_OSSE_DIR)/$(OSSE_CASE)-$(REGION)-cross.pdf))
SECONDARY_TARGETS += $(OSSE_FLUX_DECOMPOSITIONS)


## Products

$(5_RESULTS_PRODUCTS_DIR)/WOMBAT_v2S_CO2_gridded_flux_samples.nc4: \
	$(5_RESULTS_SRC_DIR)/gridded-flux-samples.R \
	$(BASIS_VECTORS) \
	$(CONTROL_EMISSIONS) \
	$(PERTURBATIONS_AUGMENTED) \
	$(SAMPLES_LNLGISSIF)
	Rscript $< \
		--basis-vectors $(BASIS_VECTORS) \
		--control-emissions $(CONTROL_EMISSIONS) \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--samples $(SAMPLES_LNLGISSIF) \
		--output $@

$(5_RESULTS_PRODUCTS_DIR)/WOMBAT_v2S_CO2_gridded_climatology_samples.nc4: \
	$(5_RESULTS_SRC_DIR)/gridded-climatology-samples.R \
	$(CONTROL_EMISSIONS) \
	$(REGION_GRID) \
	$(SIB4_CLIMATOLOGY_ASSIM_2X25) \
	$(SIB4_CLIMATOLOGY_RESP_TOT_2X25) \
	$(LANDSCHUTZER_CLIMATOLOGY_2X25) \
	$(SAMPLES_LNLGISSIF)
	Rscript $< \
		--control-emissions $(CONTROL_EMISSIONS) \
		--region-grid $(REGION_GRID) \
		--sib4-climatology-assim $(SIB4_CLIMATOLOGY_ASSIM_2X25) \
		--sib4-climatology-resp-tot $(SIB4_CLIMATOLOGY_RESP_TOT_2X25) \
		--landschutzer-climatology $(LANDSCHUTZER_CLIMATOLOGY_2X25) \
		--samples $(SAMPLES_LNLGISSIF) \
		--output $@


## Figures

$(5_RESULTS_FIGURES_DIR)/osse-true-fluxes.pdf: \
	$(5_RESULTS_SRC_DIR)/osse-true-fluxes.R \
	$(PERTURBATIONS_AUGMENTED) \
	$(OSSE_ADJUSTED_ALPHAS) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--alpha-v2 $(ALPHA_WOMBAT_V2) \
		--alpha-positive $(ALPHA_ADJUSTMENT_BASE)-positive.fst \
		--alpha-negative $(ALPHA_ADJUSTMENT_BASE)-negative.fst \
		--output $@

$(5_RESULTS_FIGURES_DIR)/osse-metrics-table-%.tex: \
	$(5_RESULTS_SRC_DIR)/osse-metrics-table.R \
	$(OSSE_FLUX_AGGREGATES_SAMPLES_CASES)
	Rscript $< \
		--metric $* \
		--flux-samples-alpha0-fixresp-wsif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHA0-FIXRESP-WSIF.rds \
		--flux-samples-alpha0-fixresp-wosif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHA0-FIXRESP-WOSIF.rds \
		--flux-samples-alpha0-freeresp-wsif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHA0-FREERESP-WSIF.rds \
		--flux-samples-alpha0-freeresp-wosif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHA0-FREERESP-WOSIF.rds \
		--flux-samples-alphav2-fixresp-wsif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHAV2-FIXRESP-WSIF.rds \
		--flux-samples-alphav2-fixresp-wosif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHAV2-FIXRESP-WOSIF.rds \
		--flux-samples-alphav2-freeresp-wsif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHAV2-FREERESP-WSIF.rds \
		--flux-samples-alphav2-freeresp-wosif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHAV2-FREERESP-WOSIF.rds \
		--flux-samples-alphap-fixresp-wsif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHAP-FIXRESP-WSIF.rds \
		--flux-samples-alphap-fixresp-wosif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHAP-FIXRESP-WOSIF.rds \
		--flux-samples-alphap-freeresp-wsif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHAP-FREERESP-WSIF.rds \
		--flux-samples-alphap-freeresp-wosif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHAP-FREERESP-WOSIF.rds \
		--flux-samples-alphan-fixresp-wsif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHAN-FIXRESP-WSIF.rds \
		--flux-samples-alphan-fixresp-wosif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHAN-FIXRESP-WOSIF.rds \
		--flux-samples-alphan-freeresp-wsif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHAN-FREERESP-WSIF.rds \
		--flux-samples-alphan-freeresp-wosif $(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-ALPHAN-FREERESP-WOSIF.rds \
		--output $@

$(5_RESULTS_FIGURES_DIR)/osse-flux-decomposition-main.pdf: \
	$(5_RESULTS_SRC_DIR)/osse-flux-decomposition.R \
	$(PERTURBATIONS_AUGMENTED) \
	$(OSSE_SAMPLES_BASE)-ALPHAN-FREERESP-WSIF.rds \
	$(OSSE_SAMPLES_BASE)-ALPHAN-FREERESP-WOSIF.rds \
	$(ALPHA_ADJUSTMENT_BASE)-negative.fst \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--samples-freeresp-wsif $(OSSE_SAMPLES_BASE)-ALPHAN-FREERESP-WSIF.rds \
		--samples-freeresp-wosif $(OSSE_SAMPLES_BASE)-ALPHAN-FREERESP-WOSIF.rds \
		--true-alpha $(ALPHA_ADJUSTMENT_BASE)-negative.fst \
		--region global \
		--trim-labels TRUE \
		--paper TRUE \
		--output $@

$(5_RESULTS_OSSE_DIR)/%-FIXRESP.pdf: \
	$(5_RESULTS_SRC_DIR)/osse-flux-decomposition.R \
	$(PERTURBATIONS_AUGMENTED) \
	$(OSSE_SAMPLES_CASES) \
	$(OSSE_ALPHAS) \
	$(DISPLAY_PARTIAL)
	Rscript $< $(OSSE_FLAGS_ALPHA) \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--samples-fixresp-wsif $(OSSE_SAMPLES_BASE)-$(firstword $(subst -, ,$*))-FIXRESP-WSIF.rds \
		--samples-fixresp-wosif $(OSSE_SAMPLES_BASE)-$(firstword $(subst -, ,$*))-FIXRESP-WOSIF.rds \
		--region $(lastword $(subst -, ,$*)) \
		--output $@

$(5_RESULTS_OSSE_DIR)/%-FREERESP.pdf: \
	$(5_RESULTS_SRC_DIR)/osse-flux-decomposition.R \
	$(PERTURBATIONS_AUGMENTED) \
	$(OSSE_SAMPLES_CASES) \
	$(OSSE_ALPHAS) \
	$(DISPLAY_PARTIAL)
	Rscript $< $(OSSE_FLAGS_ALPHA) \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--samples-freeresp-wsif $(OSSE_SAMPLES_BASE)-$(firstword $(subst -, ,$*))-FREERESP-WSIF.rds \
		--samples-freeresp-wosif $(OSSE_SAMPLES_BASE)-$(firstword $(subst -, ,$*))-FREERESP-WOSIF.rds \
		--region $(lastword $(subst -, ,$*)) \
		--output $@

$(5_RESULTS_OSSE_DIR)/%-cross.pdf: \
	$(5_RESULTS_SRC_DIR)/osse-flux-decomposition.R \
	$(PERTURBATIONS_AUGMENTED) \
	$(OSSE_SAMPLES_CASES) \
	$(OSSE_ALPHAS) \
	$(DISPLAY_PARTIAL)
	Rscript $< $(OSSE_FLAGS_ALPHA) \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--samples-fixresp-wosif $(OSSE_SAMPLES_BASE)-$(firstword $(subst -, ,$*))-FIXRESP-WOSIF.rds \
		--samples-freeresp-wsif $(OSSE_SAMPLES_BASE)-$(firstword $(subst -, ,$*))-FREERESP-WSIF.rds \
		--region $(lastword $(subst -, ,$*)) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/prior-flux-decomposition.pdf: \
	$(5_RESULTS_SRC_DIR)/prior-flux-decomposition.R \
	$(PRIOR) \
	$(CONSTRAINTS) \
	$(BASIS_VECTORS) \
	$(PERTURBATIONS_AUGMENTED)
	Rscript $< \
		--prior $(PRIOR) \
		--constraints $(CONSTRAINTS) \
		--basis-vectors $(BASIS_VECTORS) \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--region global \
		--output $@

$(5_RESULTS_FIGURES_DIR)/region-map.pdf: \
	$(5_RESULTS_SRC_DIR)/region-map.R \
	$(REGION_SF) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--region-sf $(REGION_SF) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/observation-count.pdf: \
	$(5_RESULTS_SRC_DIR)/observation-count.R \
	$(OBSERVATIONS) \
	2_matching/intermediates/runs/base/oco2-hourly.fst \
	2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
	3_sif/intermediates/oco2-hourly-sif.fst \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--observations $(OBSERVATIONS) \
		--control \
			2_matching/intermediates/runs/base/oco2-hourly.fst \
			2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
			3_sif/intermediates/oco2-hourly-sif.fst \
		--region-sf $(REGION_SF) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/sif-gpp-average-slope.pdf: \
	$(5_RESULTS_SRC_DIR)/sif-gpp-average-slope.R \
	$(MODEL_SIF_ASSIM) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--model-sif-assim $(MODEL_SIF_ASSIM) \
		--region-sf $(REGION_SF) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/sif-gpp-map-%.pdf: \
	$(5_RESULTS_SRC_DIR)/sif-gpp-map-model.R \
	$(MODEL_SIF_ASSIM) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--model-sif-assim $(MODEL_SIF_ASSIM) \
		--region-sf $(REGION_SF) \
		--term $* \
		--output $@

$(5_RESULTS_FIGURES_DIR)/flux-annual-average-table.csv: \
	$(5_RESULTS_SRC_DIR)/flux-annual-average.R \
	$(PERTURBATIONS_AUGMENTED) \
	$(SAMPLES_WOMBAT_V2) \
	$(SAMPLES_LNLGISSIF) \
	$(XBASE_MONTHLY_2x25)
	Rscript $< \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--samples-LNLGIS $(SAMPLES_WOMBAT_V2) \
		--samples-LNLGISSIF $(SAMPLES_LNLGISSIF) \
		--xbase-monthly-2x25 $(XBASE_MONTHLY_2x25) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/flux-global.pdf: \
	$(5_RESULTS_SRC_DIR)/flux-net-global-seasonal.R \
	$(PERTURBATIONS_AUGMENTED) \
	$(SAMPLES_WOMBAT_V2) \
	$(SAMPLES_LNLGISSIF) \
	$(XBASE_MONTHLY_2x25) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--samples-LNLGIS $(SAMPLES_WOMBAT_V2) \
		--samples-LNLGISSIF $(SAMPLES_LNLGISSIF) \
		--xbase-monthly-2x25 $(XBASE_MONTHLY_2x25) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/flux-decomposition.pdf: \
	$(5_RESULTS_SRC_DIR)/flux-decomposition.R \
	$(PERTURBATIONS_AUGMENTED) \
	$(SAMPLES_WOMBAT_V2) \
	$(SAMPLES_LNLGISSIF) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--samples-LNLGIS $(SAMPLES_WOMBAT_V2) \
		--samples-LNLGISSIF $(SAMPLES_LNLGISSIF) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/flux-net-zonal.pdf: \
	$(5_RESULTS_SRC_DIR)/flux-net-zonal.R \
	$(PERTURBATIONS_AUGMENTED_ZONAL) \
	$(SAMPLES_WOMBAT_V2) \
	$(SAMPLES_LNLGISSIF) \
	$(XBASE_MONTHLY_2x25_ZONAL) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--perturbations-augmented-zonal $(PERTURBATIONS_AUGMENTED_ZONAL) \
		--samples-LNLGIS $(SAMPLES_WOMBAT_V2) \
		--samples-LNLGISSIF $(SAMPLES_LNLGISSIF) \
		--xbase-monthly-2x25-zonal $(XBASE_MONTHLY_2x25_ZONAL) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/seasonal-cycle-zonal.pdf: \
	$(5_RESULTS_SRC_DIR)/seasonal-cycle-zonal.R \
	$(PERTURBATIONS_AUGMENTED_ZONAL) \
	$(SAMPLES_WOMBAT_V2) \
	$(SAMPLES_LNLGISSIF) \
	$(XBASE_MONTHLY_2x25_ZONAL) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--perturbations-augmented-zonal $(PERTURBATIONS_AUGMENTED_ZONAL) \
		--samples-LNLGIS $(SAMPLES_WOMBAT_V2) \
		--samples-LNLGISSIF $(SAMPLES_LNLGISSIF) \
		--xbase-monthly-2x25-zonal $(XBASE_MONTHLY_2x25_ZONAL) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/seasonal-latitude-profile.pdf: \
	$(5_RESULTS_SRC_DIR)/seasonal-latitude-profile.R \
	$(PERTURBATIONS_AUGMENTED) \
	$(SAMPLES_WOMBAT_V2) \
	$(SAMPLES_LNLGISSIF) \
	$(XBASE_MONTHLY_2x25) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--samples-LNLGIS $(SAMPLES_WOMBAT_V2) \
		--samples-LNLGISSIF $(SAMPLES_LNLGISSIF) \
		--xbase-monthly-2x25 $(XBASE_MONTHLY_2x25) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/average-maps-%.pdf: \
	$(5_RESULTS_SRC_DIR)/average-maps.R \
	$(SIX_YEAR_AVERAGE) \
	$(REGION_SF) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--six-year-average $(SIX_YEAR_AVERAGE) \
		--flux-component $* \
		--region-sf $(REGION_SF) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/average-maps-main-gpp.pdf: \
	$(5_RESULTS_SRC_DIR)/average-maps-gpp-fluxcom.R \
	$(SIX_YEAR_AVERAGE) \
	$(REGION_SF) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--six-year-average $(SIX_YEAR_AVERAGE) \
		--region-sf $(REGION_SF) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/average-maps-supp-gpp.pdf: \
	$(5_RESULTS_SRC_DIR)/average-maps-gpp-wombat.R \
	$(SIX_YEAR_AVERAGE) \
	$(REGION_SF) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--six-year-average $(SIX_YEAR_AVERAGE) \
		--region-sf $(REGION_SF) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/error-params-table.tex: \
	$(5_RESULTS_SRC_DIR)/error-params-table.R \
	$(HYPERPARAMETER_ESTIMATES)
	Rscript $< \
		--hyperparameter-estimates $(HYPERPARAMETER_ESTIMATES) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/traceplots.pdf: \
	$(5_RESULTS_SRC_DIR)/traceplots.R \
	$(SAMPLES_LNLGISSIF) \
	$(DISPLAY_PARTIAL)
	Rscript $< \
		--samples $(SAMPLES_LNLGISSIF) \
		--output $@

$(5_RESULTS_FIGURES_DIR)/effective-sample-size.txt: \
	$(5_RESULTS_SRC_DIR)/effective-sample-size.R \
	$(SAMPLES_LNLGISSIF)
	Rscript $< \
		--samples $(SAMPLES_LNLGISSIF) \
		--output $@


## Intermediates

$(SIX_YEAR_AVERAGE): \
	$(5_RESULTS_SRC_DIR)/six-year-average.R \
	$(PERTURBATIONS_AUGMENTED) \
	$(SAMPLES_WOMBAT_V2) \
	$(SAMPLES_LNLGISSIF) \
	$(XBASE_MONTHLY_2x25) \
	$(UTILS_PARTIAL)
	Rscript $< \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--samples-LNLGIS $(SAMPLES_WOMBAT_V2) \
		--samples-LNLGISSIF $(SAMPLES_LNLGISSIF) \
		--xbase-monthly-2x25 $(XBASE_MONTHLY_2x25) \
		--output $@

$(XBASE_MONTHLY_2x25_ZONAL): \
	$(5_RESULTS_SRC_DIR)/xbase-monthly-zonal.R \
	$(XBASE_MONTHLY_2x25) \
	$(AREA_1X1)
	Rscript $< \
		--xbase-monthly-2x25 $(XBASE_MONTHLY_2x25) \
		--area-1x1 $(AREA_1X1) \
		--output $@

$(XBASE_MONTHLY_2x25): \
	$(5_RESULTS_SRC_DIR)/xbase-monthly.R \
	$(XBASE_MONTHLY_2x25_FILES) \
	$(CONTROL_EMISSIONS)
	Rscript $< \
		--input-files $(XBASE_MONTHLY_2x25_FILES) \
		--control-emissions $(CONTROL_EMISSIONS) \
		--output $@

$(XBASE_MONTHLY_2x25_BASE)-%.nc: \
	$(GEOS_2X25_GRID) \
	$(XBASE_05X05_GRID) \
	$(XBASE_MONTHLY_FILES)
	cdo -v -z zip_6 \
		-setctomiss,0 \
		-remapcon,$(GEOS_2X25_GRID) \
		-setgrid,$(XBASE_05X05_GRID) \
		-select,name=$* \
		$(XBASE_MONTHLY_BASE)_$*_{2015,2016,2017,2018,2019,2020}.nc \
		$@

$(XBASE_MONTHLY_BASE)_NEE_%.nc:
	cdo -v -z zip_6 \
		-mul \
		-selvar,NEE \
		$(XBASE_DIRECTORY)/NEE_$*_050_monthly.nc \
		-selvar,land_fraction \
		$(XBASE_DIRECTORY)/NEE_$*_050_monthly.nc \
		$@

$(XBASE_MONTHLY_BASE)_GPP_%.nc:
	cdo -v -z zip_6 \
		-mul \
		-selvar,GPP \
		$(XBASE_DIRECTORY)/GPP_$*_050_monthly.nc \
		-selvar,land_fraction \
		$(XBASE_DIRECTORY)/GPP_$*_050_monthly.nc \
		$@

$(XBASE_MONTHLY_BASE)_TER_%.nc: \
	$(XBASE_TER_FILES)
	cdo -v -z zip_6 \
		-mul \
		-selvar,TER \
		$(XBASE_TER_BASE)_$*.nc \
		-selvar,land_fraction \
		$(XBASE_DIRECTORY)/GPP_$*_050_monthly.nc \
		$@

$(XBASE_TER_BASE)_%.nc:
	cdo -v -z zip_6 \
		-chname,NEE,TER \
		-add \
		-selvar,NEE \
		$(XBASE_DIRECTORY)/NEE_$*_050_monthly.nc \
		-selvar,GPP \
		$(XBASE_DIRECTORY)/GPP_$*_050_monthly.nc \
		$@

$(OSSE_FLUX_AGGREGATES_SAMPLES_BASE)-%.rds: \
	$(5_RESULTS_SRC_DIR)/flux-aggregates-samples.R \
	$(OSSE_SAMPLES_BASE)-%.rds \
	$(PERTURBATIONS_AUGMENTED)
	Rscript $< $(OSSE_FLAGS_ALPHA) \
		--samples $(OSSE_SAMPLES_BASE)-$*.rds \
		--perturbations-augmented $(PERTURBATIONS_AUGMENTED) \
		--output $@

$(PERTURBATIONS_AUGMENTED_ZONAL): \
	$(5_RESULTS_SRC_DIR)/perturbations-augmented-zonal.R \
	$(BASIS_VECTORS) \
	$(CONTROL_EMISSIONS) \
	$(PERTURBATIONS) \
	$(AREA_1X1)
	Rscript $< \
		--basis-vectors $(BASIS_VECTORS) \
		--control-emissions $(CONTROL_EMISSIONS) \
		--perturbations $(PERTURBATIONS) \
		--area-1x1 $(AREA_1X1) \
		--output $@

$(PERTURBATIONS_AUGMENTED): \
	$(5_RESULTS_SRC_DIR)/perturbations-augmented.R \
	$(BASIS_VECTORS) \
	$(CONTROL_EMISSIONS) \
	$(PERTURBATIONS)
	Rscript $< \
		--basis-vectors $(BASIS_VECTORS) \
		--control-emissions $(CONTROL_EMISSIONS) \
		--perturbations $(PERTURBATIONS) \
		--output $@

$(REGION_SF): \
	$(5_RESULTS_SRC_DIR)/region-sf.R \
	$(REGION_GRID)
	Rscript $< \
		--region-grid $(REGION_GRID) \
		--output $@

$(REGION_GRID): \
	$(5_RESULTS_SRC_DIR)/region-grid.R \
	$(TRANSCOM_MASK_GEOS_2X25)
	Rscript $< \
		--transcom-grid $(TRANSCOM_MASK_GEOS_2X25) \
		--output $@

$(LANDSCHUTZER_CLIMATOLOGY_2X25): \
	$(GEOS_2X25_GRID) \
	$(LANDSCHUTZER_CLIMATOLOGY)
	cdo -f nc2 remapcon,$(GEOS_2X25_GRID) $(LANDSCHUTZER_CLIMATOLOGY) $@
	ncks -A -v variable $(LANDSCHUTZER_CLIMATOLOGY) $@
