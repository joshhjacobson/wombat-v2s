$(shell mkdir -p 4_inversion/intermediates)

ALPHA_PRECISION_PARTIAL = 4_inversion/src/partials/alpha-precision.R
export ALPHA_PRECISION_PARTIAL

SIB4_CLIMATOLOGY_ASSIM_2X25 = 4_inversion/intermediates/sib4-climatology-assim-2x25.nc
SIB4_CLIMATOLOGY_RESP_TOT_2X25 = 4_inversion/intermediates/sib4-climatology-resp-tot-2x25.nc

CONTROL_EMISSIONS = 4_inversion/intermediates/control-emissions.fst
PERTURBATIONS = 4_inversion/intermediates/perturbations.fst
OBSERVATIONS = 4_inversion/intermediates/observations.fst
BASIS_VECTORS = 4_inversion/intermediates/basis-vectors.fst
CONSTRAINTS = 4_inversion/intermediates/constraints.rds
PRIOR_BASE = 4_inversion/intermediates/prior-base.rds
PRIOR = 4_inversion/intermediates/prior.rds
SENSITIVITIES_BASE_PART_1 = 4_inversion/intermediates/sensitivities-base-oco2-hourly-part-1.fst \
	4_inversion/intermediates/sensitivities-base-obspack-hourly-assim-0-part-1.fst \
	4_inversion/intermediates/sensitivities-base-obspack-hourly-assim-1-part-1.fst \
	4_inversion/intermediates/sensitivities-base-obspack-hourly-assim-2-part-1.fst \
	4_inversion/intermediates/sensitivities-base-tccon-hourly-part-1.fst \
	4_inversion/intermediates/sensitivities-base-oco2-daily-part-1.fst \
	4_inversion/intermediates/sensitivities-base-obspack-daily-assim-0-part-1.fst \
	4_inversion/intermediates/sensitivities-base-obspack-daily-assim-1-part-1.fst \
	4_inversion/intermediates/sensitivities-base-obspack-daily-assim-2-part-1.fst \
	4_inversion/intermediates/sensitivities-base-tccon-daily-part-1.fst
SENSITIVITIES_R10_R15_RNZ_PART_1 = 4_inversion/intermediates/sensitivities-r10-r15-rNZ-oco2-hourly-part-1.fst \
	4_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-hourly-assim-0-part-1.fst \
	4_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-hourly-assim-1-part-1.fst \
	4_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-hourly-assim-2-part-1.fst \
	4_inversion/intermediates/sensitivities-r10-r15-rNZ-tccon-hourly-part-1.fst \
	4_inversion/intermediates/sensitivities-r10-r15-rNZ-oco2-daily-part-1.fst \
	4_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-daily-assim-0-part-1.fst \
	4_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-daily-assim-1-part-1.fst \
	4_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-daily-assim-2-part-1.fst \
	4_inversion/intermediates/sensitivities-r10-r15-rNZ-tccon-daily-part-1.fst

H_IS = 4_inversion/intermediates/H-IS.mat.lz4
H_LNLG = 4_inversion/intermediates/H-LNLG.mat.lz4
H_SIF = 4_inversion/intermediates/H-SIF.mat.lz4

RESIDUAL_1ST_STAGE = 4_inversion/intermediates/residual-1st-stage.fst
HYPERPARAMETER_ESTIMATES = 4_inversion/intermediates/hyperparameter-estimates.fst

SAMPLES_LNLGIS = 4_inversion/intermediates/samples-LNLGIS.rds
SAMPLES_LNLGISSIF = 4_inversion/intermediates/samples-LNLGISSIF.rds

ALPHA_WOMBAT_V2 = 4_inversion/intermediates/alpha-wombat-v2.fst
ALPHA_ADJUSTMENT_CASES = positive negative
ALPHA_ADJUSTMENT_BASE = 4_inversion/intermediates/osse-alpha
OSSE_ADJUSTED_ALPHAS = $(foreach ADJUSTMENT,$(ALPHA_ADJUSTMENT_CASES),$(ALPHA_ADJUSTMENT_BASE)-$(ADJUSTMENT).fst)
OSSE_ALPHAS = $(ALPHA_WOMBAT_V2) $(OSSE_ADJUSTED_ALPHAS)

OSSE_BASE_CASES = ALPHA0 ALPHAV2 ALPHAP ALPHAN
OSSE_CASES = ALPHA0-FIXRESP-WSIF \
	ALPHA0-FIXRESP-WOSIF \
	ALPHA0-FREERESP-WSIF \
	ALPHA0-FREERESP-WOSIF \
	ALPHAV2-FIXRESP-WSIF \
	ALPHAV2-FIXRESP-WOSIF \
	ALPHAV2-FREERESP-WSIF \
	ALPHAV2-FREERESP-WOSIF \
	ALPHAP-FIXRESP-WSIF \
	ALPHAP-FIXRESP-WOSIF \
	ALPHAP-FREERESP-WSIF \
	ALPHAP-FREERESP-WOSIF \
	ALPHAN-FIXRESP-WSIF \
	ALPHAN-FIXRESP-WOSIF \
	ALPHAN-FREERESP-WSIF \
	ALPHAN-FREERESP-WOSIF
OSSE_FLAGS_ALPHA0 = --seed 0 --bio-season-slice-w 1
OSSE_FLAGS_ALPHAV2 = --seed 1 --true-alpha $(ALPHA_WOMBAT_V2)
OSSE_FLAGS_ALPHAP = --seed 2 --true-alpha $(ALPHA_ADJUSTMENT_BASE)-positive.fst
OSSE_FLAGS_ALPHAN = --seed 3 --true-alpha $(ALPHA_ADJUSTMENT_BASE)-negative.fst
OSSE_FLAGS_ALPHA = $(foreach OSSE_BASE_CASE,$(OSSE_BASE_CASES),$(OSSE_FLAGS_$(findstring $(OSSE_BASE_CASE), $*)))
OSSE_FLAGS_FREERESP = --fix-resp-linear Region03
OSSE_FLAGS_RESP = $(OSSE_FLAGS_$(findstring FREERESP, $*))
OSSE_FLAGS = $(OSSE_FLAGS_ALPHA) $(OSSE_FLAGS_RESP)
OSSE_OBSERVATIONS_BASE = 4_inversion/intermediates/osse-observations
OSSE_OBSERVATIONS_CASES = $(foreach OSSE_BASE_CASE,$(OSSE_BASE_CASES),$(OSSE_OBSERVATIONS_BASE)-$(OSSE_BASE_CASE).fst)
OSSE_SAMPLES_BASE = 4_inversion/intermediates/osse-samples
OSSE_SAMPLES_CASES = $(foreach OSSE_CASE,$(OSSE_CASES),$(OSSE_SAMPLES_BASE)-$(OSSE_CASE).rds)


# OSSE inversions

$(OSSE_SAMPLES_BASE)-%-WOSIF.rds: \
	4_inversion/src/samples.R \
	$(ALPHA_PRECISION_PARTIAL) \
	$(OSSE_OBSERVATIONS_CASES) \
	$(BASIS_VECTORS) \
	$(HYPERPARAMETER_ESTIMATES) \
	$(CONSTRAINTS) \
	$(PRIOR) \
	2_matching/intermediates/runs/base/oco2-hourly.fst \
	2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
	$(H_LNLG) \
	$(H_IS)
	Rscript $< $(OSSE_FLAGS) \
		--n-samples 1000 \
		--n-warm-up 200 \
		--observations $(OSSE_OBSERVATIONS_BASE)-$(firstword $(subst -, ,$*)).fst \
		--basis-vectors $(BASIS_VECTORS) \
		--prior $(PRIOR) \
		--constraints $(CONSTRAINTS) \
		--hyperparameter-estimates $(HYPERPARAMETER_ESTIMATES) \
		--overall-observation-mode LN LG IS \
		--control \
			2_matching/intermediates/runs/base/oco2-hourly.fst \
			2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
		--component-name LNLG IS \
		--component-parts "LN|LG" IS \
		--component-transport-matrix \
			$(H_LNLG) \
			$(H_IS) \
		--output $@

$(OSSE_SAMPLES_BASE)-%-WSIF.rds: \
	4_inversion/src/samples.R \
	$(ALPHA_PRECISION_PARTIAL) \
	$(OSSE_OBSERVATIONS_CASES) \
	$(BASIS_VECTORS) \
	$(HYPERPARAMETER_ESTIMATES) \
	$(CONSTRAINTS) \
	$(PRIOR) \
	2_matching/intermediates/runs/base/oco2-hourly.fst \
	2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
	3_sif/intermediates/oco2-hourly-sif.fst \
	$(H_LNLG) \
	$(H_IS) \
	$(H_SIF)
	Rscript $< $(OSSE_FLAGS) \
		--n-samples 1000 \
		--n-warm-up 200 \
		--observations $(OSSE_OBSERVATIONS_BASE)-$(firstword $(subst -, ,$*)).fst \
		--basis-vectors $(BASIS_VECTORS) \
		--prior $(PRIOR) \
		--constraints $(CONSTRAINTS) \
		--hyperparameter-estimates $(HYPERPARAMETER_ESTIMATES) \
		--overall-observation-mode LN LG IS LN_SIF LG_SIF \
		--control \
			2_matching/intermediates/runs/base/oco2-hourly.fst \
			2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
			3_sif/intermediates/oco2-hourly-sif.fst \
		--component-name LNLG IS SIF \
		--component-parts "LN|LG" IS "LN_SIF|LG_SIF" \
		--component-transport-matrix \
			$(H_LNLG) \
			$(H_IS) \
			$(H_SIF) \
		--output $@

$(OSSE_OBSERVATIONS_BASE)-%.fst: \
	4_inversion/src/osse-observations.R \
	$(OBSERVATIONS) \
	$(BASIS_VECTORS) \
	$(HYPERPARAMETER_ESTIMATES) \
	$(OSSE_ALPHAS) \
	2_matching/intermediates/runs/base/oco2-hourly.fst \
	2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
	3_sif/intermediates/oco2-hourly-sif.fst \
	$(H_LNLG) \
	$(H_IS) \
	$(H_SIF)
	Rscript $< $(OSSE_FLAGS_$*) \
		--basis-vectors $(BASIS_VECTORS) \
		--hyperparameter-estimates $(HYPERPARAMETER_ESTIMATES) \
		--observations $(OBSERVATIONS) \
		--overall-observation-mode LN LG IS LN_SIF LG_SIF \
		--control \
			2_matching/intermediates/runs/base/oco2-hourly.fst \
			2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
			3_sif/intermediates/oco2-hourly-sif.fst \
		--component-name LNLG IS SIF \
		--component-parts "LN|LG" IS "LN_SIF|LG_SIF" \
		--component-transport-matrix \
			$(H_LNLG) \
			$(H_IS) \
			$(H_SIF) \
		--output $@

ADJUSTED_ALPHA_DEPS = 4_inversion/src/osse-alpha-adjusted.R \
	$(REGION_MASK) \
	$(SIB4_CLIMATOLOGY_ASSIM_2X25) \
	$(SIB4_CLIMATOLOGY_RESP_TOT_2X25) \
	$(BASIS_VECTORS) \
	$(CONTROL_EMISSIONS)
ADJUSTED_ALPHA_CALL = Rscript 4_inversion/src/osse-alpha-adjusted.R \
	--region-mask $(REGION_MASK) \
	--sib4-climatology-assim $(SIB4_CLIMATOLOGY_ASSIM_2X25) \
	--sib4-climatology-resp-tot $(SIB4_CLIMATOLOGY_RESP_TOT_2X25) \
	--basis-vectors $(BASIS_VECTORS) \
	--control-emissions $(CONTROL_EMISSIONS) \
	--alpha-wombat-v2 $(ALPHA_WOMBAT_V2) \
	--fix-resp-linear Region03

$(ALPHA_ADJUSTMENT_BASE)-positive.fst: \
	$(ADJUSTED_ALPHA_DEPS)
	$(ADJUSTED_ALPHA_CALL) \
		--delta 0.1 \
		--output $@

$(ALPHA_ADJUSTMENT_BASE)-negative.fst: \
	$(ADJUSTED_ALPHA_DEPS)
	$(ADJUSTED_ALPHA_CALL) \
		--delta -0.1 \
		--output $@

$(ALPHA_WOMBAT_V2): \
	4_inversion/src/osse-alpha-v2.R \
	$(SAMPLES_WOMBAT_V2)
	Rscript $< \
		--samples-wombat-v2 $(SAMPLES_WOMBAT_V2) \
		--output $@


# Real-data inversions

$(SAMPLES_LNLGIS): \
	4_inversion/src/samples.R \
	$(ALPHA_PRECISION_PARTIAL) \
	$(OBSERVATIONS) \
	$(BASIS_VECTORS) \
	$(HYPERPARAMETER_ESTIMATES) \
	$(CONSTRAINTS) \
	$(PRIOR) \
	2_matching/intermediates/runs/base/oco2-hourly.fst \
	2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
	$(H_LNLG) \
	$(H_IS)
	Rscript $< \
		--observations $(OBSERVATIONS) \
		--basis-vectors $(BASIS_VECTORS) \
		--prior $(PRIOR) \
		--constraints $(CONSTRAINTS) \
		--hyperparameter-estimates $(HYPERPARAMETER_ESTIMATES) \
		--overall-observation-mode LN LG IS \
		--control \
		    2_matching/intermediates/runs/base/oco2-hourly.fst \
				2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
		--component-name LNLG IS \
		--component-parts "LN|LG" IS \
		--component-transport-matrix \
			$(H_LNLG) \
			$(H_IS) \
		--output $@

## NOTE(jhj): In the inversion with SIF data, we fix the respiration linear component for Region03 only;
## in the other inversions, we default to fixing it for all land regions.
$(SAMPLES_LNLGISSIF): \
	4_inversion/src/samples.R \
	$(ALPHA_PRECISION_PARTIAL) \
	$(OBSERVATIONS) \
	$(BASIS_VECTORS) \
	$(HYPERPARAMETER_ESTIMATES) \
	$(CONSTRAINTS) \
	$(PRIOR) \
	2_matching/intermediates/runs/base/oco2-hourly.fst \
	2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
	3_sif/intermediates/oco2-hourly-sif.fst \
	$(H_LNLG) \
	$(H_IS) \
	$(H_SIF)
	Rscript $< \
		--fix-resp-linear Region03 \
		--observations $(OBSERVATIONS) \
		--basis-vectors $(BASIS_VECTORS) \
		--prior $(PRIOR) \
		--constraints $(CONSTRAINTS) \
		--hyperparameter-estimates $(HYPERPARAMETER_ESTIMATES) \
		--overall-observation-mode LN LG IS LN_SIF LG_SIF \
		--control \
			2_matching/intermediates/runs/base/oco2-hourly.fst \
			2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
			3_sif/intermediates/oco2-hourly-sif.fst \
		--component-name LNLG IS SIF \
		--component-parts "LN|LG" IS "LN_SIF|LG_SIF" \
		--component-transport-matrix \
			$(H_LNLG) \
			$(H_IS) \
			$(H_SIF) \
		--output $@


# Hyperparameter estimates

$(HYPERPARAMETER_ESTIMATES): \
	4_inversion/src/hyperparameter-estimates.R \
	$(OBSERVATIONS) \
	$(RESIDUAL_1ST_STAGE)
	Rscript $< \
		--observations $(OBSERVATIONS) \
		--residuals $(RESIDUAL_1ST_STAGE) \
		--output $@

## Inversions (1st stage) to get a residual for hyperparameter estimates

$(RESIDUAL_1ST_STAGE): \
	4_inversion/src/residual.R \
	2_matching/intermediates/runs/base/oco2-hourly.fst \
	2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
	3_sif/intermediates/oco2-hourly-sif.fst \
	$(OBSERVATIONS) \
	$(BASIS_VECTORS) \
	$(PRIOR) \
	$(CONSTRAINTS) \
	$(H_LNLG) \
	$(H_IS) \
	$(H_SIF)
	Rscript $< \
		--fix-resp-linear Region03 \
		--observations $(OBSERVATIONS) \
		--basis-vectors $(BASIS_VECTORS) \
		--prior $(PRIOR) \
		--constraints $(CONSTRAINTS) \
		--overall-observation-mode LN LG IS LN_SIF LG_SIF \
		--control \
			2_matching/intermediates/runs/base/oco2-hourly.fst \
			2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
			3_sif/intermediates/oco2-hourly-sif.fst \
		--component-name LNLG IS SIF \
		--component-parts "LN|LG" IS "LN_SIF|LG_SIF" \
		--component-transport-matrix $(H_LNLG) $(H_IS) $(H_SIF) \
		--output $@

# Transport matrix (H)

TRANSPORT_MATRIX_DEPS = 4_inversion/src/transport-matrix.R \
	4_inversion/intermediates/sensitivities-base-oco2-hourly-part-1.fst \
	4_inversion/intermediates/sensitivities-base-obspack-hourly-assim-1-part-1.fst \
	4_inversion/intermediates/sensitivities-base-oco2-daily-part-1.fst \
	4_inversion/intermediates/sensitivities-base-obspack-daily-assim-1-part-1.fst \
	4_inversion/intermediates/sensitivities-r10-r15-rNZ-oco2-hourly-part-1.fst \
	4_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-hourly-assim-1-part-1.fst \
	4_inversion/intermediates/sensitivities-r10-r15-rNZ-oco2-daily-part-1.fst \
	4_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-daily-assim-1-part-1.fst \
	$(BASIS_VECTORS) \
	$(OBSERVATIONS)
TRANSPORT_MATRIX_CALL = Rscript 4_inversion/src/transport-matrix.R \
	--basis-vectors $(BASIS_VECTORS) \
	--observations $(OBSERVATIONS) \
	--sensitivities-base \
		4_inversion/intermediates/sensitivities-base \
		4_inversion/intermediates/sensitivities-r10-r15-rNZ

$(H_IS): \
	$(TRANSPORT_MATRIX_DEPS)
	$(TRANSPORT_MATRIX_CALL) \
		--overall-observation-mode IS \
		--control 2_matching/intermediates/runs/base/obspack-hourly-assim-1.fst \
		--output $@

$(H_LNLG): \
	$(TRANSPORT_MATRIX_DEPS)
	$(TRANSPORT_MATRIX_CALL) \
		--overall-observation-mode LN LG \
		--control 2_matching/intermediates/runs/base/oco2-hourly.fst \
		--output $@

$(H_SIF): \
	4_inversion/src/transport-matrix.R \
	$(SENSITIVITIES_SIF) \
	$(BASIS_VECTORS) \
	$(OBSERVATIONS)
	Rscript 4_inversion/src/transport-matrix.R \
		--basis-vectors $(BASIS_VECTORS) \
		--observations $(OBSERVATIONS) \
		--sensitivities-base 3_sif/intermediates/sensitivities \
		--overall-observation-mode LN_SIF LG_SIF \
		--control 3_sif/intermediates/oco2-hourly-sif.fst \
		--output $@

# Sensitivities

$(SENSITIVITIES_R10_R15_RNZ_PART_1) &: \
	4_inversion/src/sensitivities.R
	Rscript $< \
		--input \
			oco2-hourly.fst \
			obspack-hourly-assim-1.fst \
			obspack-hourly-assim-2.fst \
			tccon-hourly.fst \
			oco2-daily.fst \
			obspack-daily-assim-1.fst \
			obspack-daily-assim-2.fst \
			tccon-daily.fst \
			obspack-hourly-assim-0.fst \
			obspack-daily-assim-0.fst \
		--resolution hourly hourly hourly hourly daily daily daily daily \
			hourly daily \
		--runs 1_transport/intermediates/runs-r10-r15-rNZ \
		--matched-runs 2_matching/intermediates/runs-r10-r15-rNZ \
		--output-base 4_inversion/intermediates/sensitivities-r10-r15-rNZ-oco2-hourly \
			4_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-hourly-assim-1 \
			4_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-hourly-assim-2 \
			4_inversion/intermediates/sensitivities-r10-r15-rNZ-tccon-hourly \
			4_inversion/intermediates/sensitivities-r10-r15-rNZ-oco2-daily \
			4_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-daily-assim-1 \
			4_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-daily-assim-2 \
			4_inversion/intermediates/sensitivities-r10-r15-rNZ-tccon-daily \
			4_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-hourly-assim-0 \
			4_inversion/intermediates/sensitivities-r10-r15-rNZ-obspack-daily-assim-0

$(SENSITIVITIES_BASE_PART_1) &: \
	4_inversion/src/sensitivities.R
	Rscript $< \
		--input \
			oco2-hourly.fst \
			obspack-hourly-assim-1.fst \
			obspack-hourly-assim-2.fst \
			tccon-hourly.fst \
			oco2-daily.fst \
			obspack-daily-assim-1.fst \
			obspack-daily-assim-2.fst \
			tccon-daily.fst \
			obspack-hourly-assim-0.fst \
			obspack-daily-assim-0.fst \
		--resolution hourly hourly hourly hourly daily daily daily daily \
			hourly daily \
		--runs 1_transport/intermediates/runs \
		--matched-runs 2_matching/intermediates/runs \
		--output-base 4_inversion/intermediates/sensitivities-base-oco2-hourly \
			4_inversion/intermediates/sensitivities-base-obspack-hourly-assim-1 \
			4_inversion/intermediates/sensitivities-base-obspack-hourly-assim-2 \
			4_inversion/intermediates/sensitivities-base-tccon-hourly \
			4_inversion/intermediates/sensitivities-base-oco2-daily \
			4_inversion/intermediates/sensitivities-base-obspack-daily-assim-1 \
			4_inversion/intermediates/sensitivities-base-obspack-daily-assim-2 \
			4_inversion/intermediates/sensitivities-base-tccon-daily \
			4_inversion/intermediates/sensitivities-base-obspack-hourly-assim-0 \
			4_inversion/intermediates/sensitivities-base-obspack-daily-assim-0

# Preliminaries

$(PRIOR): \
	4_inversion/src/prior.R \
	$(ALPHA_PRECISION_PARTIAL) \
	$(PRIOR_BASE) \
	$(BASIS_VECTORS) \
	$(CONTROL_EMISSIONS) \
	$(PERTURBATIONS)
	Rscript $< \
		--prior-base $(PRIOR_BASE) \
		--basis-vectors $(BASIS_VECTORS) \
		--control-emissions $(CONTROL_EMISSIONS) \
		--perturbations $(PERTURBATIONS) \
		--output $@

$(PRIOR_BASE): \
	4_inversion/src/prior-base.R \
	$(PERTURBATIONS) \
	$(BASIS_VECTORS)
	# NOTE(mgnb/jhj): the GpGp package may complain if the number of threads is greater
	# than one; try setting OMP_NUM_THREADS=1 as needed
	Rscript $< \
		--basis-vectors $(BASIS_VECTORS) \
		--perturbations $(PERTURBATIONS) \
		--output $@

$(CONSTRAINTS): \
	4_inversion/src/constraints.R \
	$(BASIS_VECTORS) \
	$(CONTROL_EMISSIONS) \
	$(PERTURBATIONS)
	Rscript $< \
		--basis-vectors $(BASIS_VECTORS) \
		--control-emissions $(CONTROL_EMISSIONS) \
		--perturbations $(PERTURBATIONS) \
		--output $@

$(BASIS_VECTORS): \
	4_inversion/src/basis-vectors.R \
	$(PERTURBATIONS)
	Rscript $< \
		--perturbations $(PERTURBATIONS) \
		--output $@

$(PERTURBATIONS): \
	4_inversion/src/perturbations.R
	Rscript $< \
		--runs 1_transport/intermediates/runs 1_transport/intermediates/runs-r10-r15-rNZ \
		--matched-runs 2_matching/intermediates/runs 2_matching/intermediates/runs-r10-r15-rNZ \
		--output $@

$(CONTROL_EMISSIONS): \
	4_inversion/src/control-emissions.R
	Rscript $< \
		--matched-runs 2_matching/intermediates/runs \
		--output $@

$(OBSERVATIONS): \
	4_inversion/src/observations.R \
	$(OCO2_OBSERVATIONS) \
	$(OCO2_OBSERVATIONS_SIF) \
	$(CONTROL_SIF)
	Rscript $< \
		--oco2-observations $(OCO2_OBSERVATIONS) \
		--obspack-directory $(OBSPACK_DIRECTORY) \
		--tccon-sounding-directory $(TCCON_SOUNDING_DIRECTORY) \
		--oco2-observations-sif $(OCO2_OBSERVATIONS_SIF) \
		--control-sif $(CONTROL_SIF) \
		--start-date $(INVERSION_START_DATE) \
		--end-date $(INVERSION_END_DATE) \
		--output $@

$(SIB4_CLIMATOLOGY_ASSIM_2X25): \
	$(GEOS_2X25_GRID) \
	$(SIB4_CLIMATOLOGY_ASSIM)
	cdo -f nc2 remapcon,$(GEOS_2X25_GRID) $(SIB4_CLIMATOLOGY_ASSIM) $@
	ncks -A -v variable $(SIB4_CLIMATOLOGY_ASSIM) $@

$(SIB4_CLIMATOLOGY_RESP_TOT_2X25): \
	$(GEOS_2X25_GRID) \
	$(SIB4_CLIMATOLOGY_RESP_TOT)
	cdo -f nc2 remapcon,$(GEOS_2X25_GRID) $(SIB4_CLIMATOLOGY_RESP_TOT) $@
	ncks -A -v variable $(SIB4_CLIMATOLOGY_RESP_TOT) $@
