$(shell mkdir -p 3_sif/intermediates)
$(shell mkdir -p 3_sif/figures)

SIB4_SIF_HOURLY = $(foreach SIB4_YEAR,$(SIB4_INPUT_YEARS),3_sif/intermediates/sib4-hourly-sif-$(SIB4_YEAR).nc)
SIB4_SIF_ASSIM_HOURLY_2X25 = $(foreach SIB4_YEAR,$(SIB4_INPUT_YEARS),3_sif/intermediates/sib4-hourly-sif-assim-2x25-$(SIB4_YEAR).nc)
MODEL_SIF_ASSIM = 3_sif/intermediates/model-sif-assim.fst

OCO2_OBSERVATIONS_SIF = 3_sif/intermediates/observations-sif.fst
CONTROL_SIF = 3_sif/intermediates/oco2-hourly-sif.fst

SIB4_CLIMATOLOGY_INVENTORY_ASSIM_HOURLY_2X25 = \
	$(foreach SIB4_YEAR,$(INVENTORY_OUTPUT_YEARS),\
	3_sif/intermediates/sib4-hourly-climatology-inventory-assim-2x25-$(SIB4_YEAR).nc)
SIB4_RESIDUAL_ASSIM_HOURLY_2X25 = \
	$(foreach SIB4_YEAR,$(SIB4_INPUT_YEARS),3_sif/intermediates/sib4-hourly-residual-assim-2x25-$(SIB4_YEAR).nc)
SENSITIVITIES_SIF = 3_sif/intermediates/sensitivities-oco2-hourly-sif.fst

3_SIF_TARGETS = \
	3_sif/figures/map-fit-statistics.pdf \
	3_sif/figures/validate-matches.pdf \
	$(SENSITIVITIES_SIF)

# Sensitivities

$(SENSITIVITIES_SIF) &: \
	3_sif/src/sensitivities-sif.R \
	$(REGION_MASK) \
	$(CONTROL_SIF) \
	$(SIB4_CLIMATOLOGY_INVENTORY_ASSIM_HOURLY_2X25) \
	$(SIB4_RESIDUAL_ASSIM_HOURLY_2X25)
	Rscript $< \
		--region-mask $(REGION_MASK) \
		--control-sif $(CONTROL_SIF) \
		--basis-climatology $(SIB4_CLIMATOLOGY_INVENTORY_ASSIM_HOURLY_2X25) \
		--basis-residual $(SIB4_RESIDUAL_ASSIM_HOURLY_2X25) \
		--output $@

## Regrid time-series components to basis grid

$(SIB4_CLIMATOLOGY_INVENTORY_ASSIM_HOURLY_2X25) &: \
	3_sif/src/regrid-decomposition.sh \
	$(GEOS_2X25_GRID) \
	$(SIB4_CLIMATOLOGY_INVENTORY_ASSIM_HOURLY)
	bash $< \
		$(GEOS_2X25_GRID) \
		"$(SIB4_CLIMATOLOGY_INVENTORY_ASSIM_HOURLY)" \
		3_sif/intermediates

$(SIB4_RESIDUAL_ASSIM_HOURLY_2X25) &: \
	3_sif/src/regrid-decomposition.sh \
	$(GEOS_2X25_GRID) \
	$(SIB4_RESIDUAL_ASSIM_HOURLY)
	bash $< \
		$(GEOS_2X25_GRID) \
		"$(SIB4_RESIDUAL_ASSIM_HOURLY)" \
		3_sif/intermediates

# Control SIF

3_sif/figures/validate-matches.pdf: \
	3_sif/src/validate-matches.Rmd \
	$(OCO2_OBSERVATIONS_SIF) \
	$(CONTROL_SIF) \
	$(SIB4_SIF_ASSIM_HOURLY_2X25)
	Rscript -e "rmarkdown::render('$<', output_file = '$(shell basename $@)', output_dir = '$(shell dirname $@)', knit_root_dir = '../../', params = list(oco2_observations_sif = '$(OCO2_OBSERVATIONS_SIF)', control_sif = '$(CONTROL_SIF)', inventory = '$(SIB4_SIF_ASSIM_HOURLY_2X25)'))"

$(CONTROL_SIF): \
	3_sif/src/match-sif.R \
	$(OCO2_OBSERVATIONS_SIF) \
	$(SIB4_SIF_ASSIM_HOURLY_2X25) \
	$(MODEL_SIF_ASSIM)
	Rscript $< \
		--oco2-observations-sif $(OCO2_OBSERVATIONS_SIF) \
		--inventory $(SIB4_SIF_ASSIM_HOURLY_2X25) \
		--linear-models $(MODEL_SIF_ASSIM) \
		--output $@

# Fitted slope coefficients

3_sif/figures/map-fit-statistics.pdf: \
	3_sif/src/map-fit-statistics.Rmd \
	$(MODEL_SIF_ASSIM)
	Rscript -e "rmarkdown::render('$<', output_file = '$(shell basename $@)', output_dir = '$(shell dirname $@)', knit_root_dir = '../../', params = list(model_sif_assim = '$(MODEL_SIF_ASSIM)'))"

$(MODEL_SIF_ASSIM): \
	3_sif/src/linear-model.R \
	$(SIB4_SIF_ASSIM_HOURLY_2X25)
	Rscript $< \
		--input $(SIB4_SIF_ASSIM_HOURLY_2X25) \
		--output $@

# SIF and ASSIM inventories

$(SIB4_SIF_ASSIM_HOURLY_2X25) &: \
	3_sif/src/regrid-sif-assim.sh \
	$(GEOS_2X25_GRID) \
	$(SIB4_HOURLY) \
	$(SIB4_SIF_HOURLY)
	bash $< \
		$(GEOS_2X25_GRID) \
		"$(SIB4_HOURLY)" \
		"$(SIB4_SIF_HOURLY)" \
		3_sif/intermediates

# SiB4 SIF units are listed as W/m2/nm/sr but are actually W/m^2/µm/sr
# See SIF section in SiB4 technical description: https://hdl.handle.net/10217/200691
# Also see the related paper on the canopy scaling method: https://doi.org/10.1111/gcb.12948

$(SIB4_SIF_HOURLY) &: \
	$(SIB4_HOURLY_DIRECTORY)
	cdo -w -z zip_6 \
		-splityear \
		-shifttime,-30minutes \
		-settunits,minutes \
		-setattribute,sif@units=W/m^2/µm/sr \
		-vertsum \
		-select,name=sif \
		$(SIB4_HOURLY_DIRECTORY)/sib4-hourly-{2014,2015,2016,2017,2018,2019,2020}-* \
		3_sif/intermediates/sib4-hourly-sif-

# Observations (units: W/m^2/µm/sr)

$(OCO2_OBSERVATIONS_SIF): \
	3_sif/src/observations-sif.R \
	$(OCO2_SIF_DIRECTORY)
	Rscript $< \
		--oco2-sif-directory $(OCO2_SIF_DIRECTORY) \
		--start-date $(INVERSION_START_DATE) \
		--end-date $(INVERSION_END_DATE) \
		--output $@
