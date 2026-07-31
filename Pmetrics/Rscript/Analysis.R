# =============================================================================
# FDA-FEP-ELF :: Analysis.R #####
# Canonical analysis script for cefepime CRRT / ELF Pmetrics workflow
#
# Rules
# - Open from FDA-FEP-ELF.Rproj
# - Do not use setwd()
# - Use repository-root-relative paths only
# - Keep source and derived analysis data in Pmetrics/src
# - Keep simulation inputs in Pmetrics/Sim
# - Keep run folders in Pmetrics/Runs
# - Keep manuscript-ready outputs in Manuscript/Figures and Manuscript/Tables
# - Preserve scientific behavior unless explicitly approved otherwise
# =============================================================================

# ---- project paths -----------------------------------------------------------
PATHS <- list(
  root = ".",
  src  = "src",
  sim  = "Sim",
  runs = "Runs",
  r    = "Rscript",
  fig  = "Manuscript/Figures",
  tab  = "Manuscript/Tables",
  arch = "archive"
)

dir.create(PATHS$fig, recursive = TRUE, showWarnings = FALSE)
dir.create(PATHS$tab, recursive = TRUE, showWarnings = FALSE)

# ---- optional support sourcing ----------------------------------------------
# source(file.path(PATHS$r, "utils_io_plotting.R"), chdir = FALSE)
# source(file.path(PATHS$r, "utils_tables.R"), chdir = FALSE)

reticulate::use_condaenv("r-reticulate")
reticulate::py_run_string("import sys")

#devtools::install_github("LAPKB/Pmetrics_rust")

library(Pmetrics)
library(tidyverse)

# utils
get_run <- function(run_no, env = parent.frame()) {
  run_name <- paste0("run", run_no)
  
  if (exists(run_name, envir = env, inherits = TRUE)) {
    get(run_name, envir = env, inherits = TRUE)
  } else {
    obj <- PM_load(run_no, path = PATHS$runs)
    assign(run_name, obj, envir = .GlobalEnv)
    obj
  }
}

# =============================================================================
# Dataset and run roadmap ####
# =============================================================================
#
# Dataset objects:
#   dat  = raw Excel import
#   dat1 = rich CRRT plasma/post-filter/effluent dataset before PM_data creation
#   dat3 = PM_data object from dat1/dat_ML.csv; used for rich CRRT structural screening
#   dat4 = combined CRRT + ELF dataset before train/test split
#   dat5 = training PM_data object from train2.csv
#   dat6 = validation PM_data object from test2.csv
#
# Run registry:
#   run4  = mod_cef1, rich CRRT, base linear CrCL
#   run5  = mod_cef2, rich CRRT, BSA scalar
#   run6  = mod_cef3, rich CRRT, CrCL power
#   run7  = mod_cef4, rich CRRT, WT allometry on CL
#   run8  = mod_cef5, rich CRRT, CrCL power + WT allometry on CL
#   run9  = mod_cef6, rich CRRT, CrCL power + WT allometry + CRRT saturation
#   run10 = mod_cef7, ELF-extended mod_cef5, training
#   run11 = mod_cef7, ELF-extended mod_cef5, validation
#   run12 = mod_cef8, ELF-extended mod_cef6, training
#   run13 = mod_cef8, ELF-extended mod_cef6, validation
#   run14/run15 = mod_cef1_elf training/validation
#   run16/run17 = mod_cef2_elf training/validation
#   run18/run19 = mod_cef3_elf training/validation
#   run20/run21 = mod_cef4_elf training/validation
#
# Model-selection logic:
#   1. Use dat3 to screen plasma/CRRT structural candidates.
#   2. Build dat4, then split into dat5/dat6 for combined plasma + ELF modeling.
#   3. Fit ELF-extended candidates on dat5.
#   4. Validate candidates on dat6 using cycles = 0 and the corresponding training prior.
# =============================================================================

# =============================================================================
# 1. Build rich CRRT dataset for plasma/CRRT structural screening ####
# =============================================================================
# Read in data file for training base CRRT model structure ####

dat <- readxl::read_excel(file.path(PATHS$src, "FEP_comb_ELF_old_ML_combined_FDA.xlsx"),sheet="FEP_comb_ELF_old")

dat1 <- dat %>%
  # Drop any flagged rows first
  filter(EXCLUDE_DV != 1) %>%
  filter(RRT == 1, RX_include == 1) %>%
  rename(
    DOSE  = AMT,
    DUR   = TINF,
    INPUT = ADM,
    OUTEQ = DVID,
    OUT   = DV,
    LOQ   = LIMIT # change from LOQ to CENS (0 for no, 1 for yes, -1 for ALQ)
  ) %>%
  mutate(
    OUTEQ = na_if(OUTEQ, ".") |> as.numeric(),
    LOQ   = na_if(LOQ, ".")   |> as.numeric(),
    OUT   = na_if(OUT, ".")   |> as.numeric()
  ) %>%
  group_by(OUTEQ) %>%
  mutate( # this needs to translate to the observation (OUT) column as LOQ
    LOQ = case_when(
      EVID == 1 ~ NA_real_,
      EVID == 0 & OUTEQ == 1 & is.na(LOQ) ~ 0.5,
      EVID == 0 & OUTEQ == 3 & is.na(LOQ) ~ 0.5,
      EVID == 0 & OUTEQ == 4 & is.na(LOQ) ~ 0.5,
      EVID == 0 & OUTEQ == 5 & is.na(LOQ) ~ 0.5,
      TRUE ~ LOQ
    )
  ) %>%
  ungroup() %>%
  mutate(             # not sure if we will need to remap the outeqs
    OUTEQ = case_when(
      OUTEQ==1 ~ 1,
      OUTEQ==3 ~ 2,
      OUTEQ==4 ~ 3,
      OUTEQ==5 ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    OUT = case_when(
      EVID == 0 & MDV == 1  ~ -99,
      EVID == 0 & is.na(OUT) ~ -99,
      EVID == 0 & OUT == "." ~ -99,
      TRUE ~ OUT
    )
  ) %>%
  mutate(
    BAG_RESET = case_when(
      INPUT == 2 ~ 0,
      TRUE ~ BAG_RESET
    )
  ) %>%
  select(ID,EVID,TIME,DOSE,DUR,ADDL,II, INPUT,OUT,OUTEQ,everything(),
         -RX_ID,-EXCLUDE_DV,-DATE,-ALIQUOT_ID,-no_RRT,-RX_include, 
         -no_RX_exclude, -DATE, -HTIN, -CENS, -MDV, -RRT, -LOQ,
         -PLA_urea,-BAL_urea,-BAL_PK) %>%
  arrange(ID, TIME) %>%
  rename_with(tolower) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), ".", as.character(.))))

write.csv(dat1,file.path(PATHS$src, "dat_ML.csv"),row.names = F)

dat3 <- PM_data$new(file.path(PATHS$src, "dat_ML.csv"),loq=c(0,0,0,0))

#dat3$plot(outeq=1:3,tad=T,log=T)

# # --- Cefepime CVVH/CVVHDF pop PK in Pmetrics (matches your mlxtran) -----------------


# =============================================================================
# 2. Fit rich CRRT plasma/CRRT structural candidates: runs 4-9 ####
# =============================================================================

# Run 4 Base model CrCL on CL, WT on Vd, aligns with monolix FDA_run_14.mlxtran ####

mod_cef1 <- PM_model$new(
  # ---------------- Priors / bounds ----------------
  pri = list(
    V1     = ab(0.1,  30),   # L (typical central, prior wide; scaled by weight in eqn)
    K12    = ab(  0,  10),   # 1/h
    K21    = ab(  0,  10),   # 1/h
    CL1    = ab(0.1,  12),   # L/h (renal scale at CRCL=120)
    S_eff  = ab(  0,  2),    # dimensionless sieving scaler for effluent clearance
    S_post = ab(  0,  2)     # dimensionless scaler for post-filter clearance
  ),
  
  # ---------------- Covariates (names must match your dataset column headers) -----
  cov = list(
    interval_volume = interp("none"),  # L
    age             = interp(),        # years
    male            = interp("none"),  # 0/1
    ht              = interp(),        # cm or m
    wt              = interp(),        # kg
    scr             = interp(),        # mg/dL
    crcl            = interp(),        # mL/min
    bsa             = interp(),        # m^2
    crcl_bsa        = interp(),        # mL/min/1.73 m^2
    ecmo            = interp("none"),  # 0/1
    hd              = interp("none"),  # 0/1
    crrt            = interp("none"),  # 0/1
    cvvh            = interp("none"),  # 0/1
    cvvhd           = interp("none"),  # 0/1
    cvvhdf          = interp("none"),  # 0/1
    flow            = interp("none"),  # mL/h
    bfr             = interp("none"),  # mL/min (blood flow rate)
    bag_reset       = interp("none")   # 0/1
  ),
  sec = function(){
    CL2 = 7.2                          # Fixed HD clearance
    
    # Weight-scaled central volume
    V_WT = V1 * (wt/70)                # WT scaled Vd
    
    # Systemic (renal ± HD)
    CL_R   = CL1 * (crcl/120)                         # L/h
    CL_HD  = CL2 * hd                                 # L/h
    CL_sys = CL_R + CL_HD                             # L/h
    
    # CRRT clearance term from effluent flow (mL/h -> L/h)
    Qe      = (flow/1000) * crrt      # L/h; 0 when CRRT=0
    CL_CRRT = S_eff * Qe              # L/h
    
    # Total clearance and elimination from central
    CLT = CL_sys + CL_CRRT            # L/h
    Ke  = CLT / V_WT                  # 1/h
    
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    
    # mg/L pre-filter
    Cpre  = X[1]/V_WT
    # mg/L post-filter
    Cpost = Cpre * S_post * crrt
    
    Aeff = 0
    Ceff = 0
    
    # Guard against divide-by-zero
    if (X[4] > 1e-6) {
      Ceff = X[3] / X[4]
    }
    
    if (X[3] > 1e-6) {
      Aeff = X[3]
    }
  },
  # ---------------- Differential equations ----------------
  eqn = function(){
    # Two compartment primary model structure
    dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2]
    dX[2] =                  K12*X[1] - K21*X[2]
    # Bag accumulation since last dump
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    dX[3] = Qe * S_eff * (X[1]/V_WT) * Iacc  - kdump * X[3] * (1 - Iacc)
    dX[4] = Qe * Iacc                        - kdump * X[4] * (1 - Iacc)
  },
  
  # ---------------- Outputs ----------------
  out = function(){
    Y[1]  = Cpre       # DVID=1: pre/regular plasma in mg/L
    Y[2]  = Cpost      # DVID=3: post-filter plasma (0 if CRRT off)
    Y[3]  = Ceff       # DVID=4: effluent concentration (0 if not collecting)
    Y[4]  = Aeff       # DVID=5: effluent cumulative amount (mg)
  },
  
  # ---------------- Error models ----------------
  # Proportional with a floor per-assay.
  err = list(
    proportional(2, c(1, 0.15, 0, 0)),   # Y1 (Cpre)
    proportional(2, c(1, 0.15, 0, 0)),   # Y2 (Cpost)
    proportional(2, c(1, 0.15, 0, 0)),   # Y3 (Ceff_out)
    proportional(2, c(1, 0.15, 0, 0))    # Y4 (Aeff)
  )
)

run4 <- mod_cef1$fit(data = dat3,
                     cycles = 500,
                     path = PATHS$runs,
                     run = 4,
                     points = 250,
                     overwrite=TRUE)

run4 <- PM_load(4, path = PATHS$runs)

# Run 5, evaluating run4 str with bsa power scalar on clearance, aligns with monolix FDA_run_15.mlxtran ####

mod_cef2 <- PM_model$new(
  # ---------------- Priors / bounds ----------------
  pri = list(
    V1     = ab(0.1,  30),    # L (typical central, prior wide; scaled by weight in eqn)
    K12    = ab(  0,  10),    # 1/h
    K21    = ab(  0,  10),    # 1/h
    CL1    = ab(0.1,  12),    # L/h (renal scale at CRCL=120)
    S_eff  = ab(  0,  2),     # dimensionless sieving scaler for effluent clearance
    S_post = ab(  0,  2),     # dimensionless scaler for post-filter clearance
    theta1 = ab(  0,  2)      # BSA clearance scaler
  ),
  
  # ---------------- Covariates (names must match your dataset column headers) -----
  cov = list(
    interval_volume = interp("none"),  # L
    age             = interp(),        # years
    male            = interp("none"),  # 0/1
    ht              = interp(),        # cm or m
    wt              = interp(),        # kg
    scr             = interp(),        # mg/dL
    crcl            = interp(),        # mL/min
    bsa             = interp(),        # m^2
    crcl_bsa        = interp(),        # mL/min/1.73 m^2
    ecmo            = interp("none"),  # 0/1
    hd              = interp("none"),  # 0/1
    crrt            = interp("none"),  # 0/1
    cvvh            = interp("none"),  # 0/1
    cvvhd           = interp("none"),  # 0/1
    cvvhdf          = interp("none"),  # 0/1
    flow            = interp("none"),  # mL/h
    bfr             = interp("none"),  # mL/min (blood flow rate)
    bag_reset       = interp("none")   # 0/1
  ),
  sec = function(){
    CL2 = 7.2                          # fixed HD clearance
    
    # Weight-scaled central volume
    V_WT = V1 * (wt/70)                # WT scaled Vd
    
    # Systemic (renal ± HD)
    CL_R   = CL1 * (crcl/120) * (bsa/1.73)**theta1                        # L/h
    CL_HD  = CL2 * hd                                                     # L/h
    CL_sys = CL_R + CL_HD                                                 # L/h
    
    # CRRT clearance term from effluent flow (mL/h -> L/h)
    Qe      = (flow/1000) * crrt      # L/h; 0 when CRRT=0
    CL_CRRT = S_eff * Qe              # L/h
    
    # Total clearance and elimination from central
    CLT = CL_sys + CL_CRRT            # L/h
    Ke  = CLT / V_WT                  # 1/h
    
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    
    # mg/L pre-filter
    Cpre  = X[1]/V_WT
    # mg/L post-filter
    Cpost = Cpre * S_post * crrt
    
    Aeff = 0
    Ceff = 0
    
    # Guard against divide-by-zero
    if (X[4] > 1e-6) {
      Ceff = X[3] / X[4]
    }
    
    if (X[3] > 1e-6) {
      Aeff = X[3]
    }
  },
  # ---------------- Differential equations ----------------
  eqn = function(){
    # Two compartment primary model structure
    dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2]
    dX[2] =                  K12*X[1] - K21*X[2]
    # Bag accumulation since last dump
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    dX[3] = Qe * S_eff * (X[1]/V_WT) * Iacc  - kdump * X[3] * (1 - Iacc)
    dX[4] = Qe * Iacc                        - kdump * X[4] * (1 - Iacc)
  },
  
  # ---------------- Outputs ----------------
  out = function(){
    Y[1]  = Cpre       # DVID=1: pre/regular plasma in mg/L
    Y[2]  = Cpost      # DVID=3: post-filter plasma (0 if CRRT off)
    Y[3]  = Ceff       # DVID=4: effluent concentration (0 if not collecting)
    Y[4]  = Aeff       # DVID=5: effluent cumulative amount (mg)
  },
  
  # ---------------- Error models ----------------
  # Proportional with a floor per-assay.
  err = list(
    proportional(2, c(1, 0.15, 0, 0)),   # Y1 (Cpre)
    proportional(2, c(1, 0.15, 0, 0)),   # Y2 (Cpost)
    proportional(2, c(1, 0.15, 0, 0)),   # Y3 (Ceff_out)
    proportional(2, c(1, 0.15, 0, 0))    # Y4 (Aeff)
  )
)

run5 <- mod_cef2$fit(data = dat3,
                     cycles = 500,
                     path = PATHS$runs,
                     run = 5,
                     points = 250,
                     overwrite=TRUE)

run5 <- PM_load(5, path = PATHS$runs)

# Run 6, Updating model to estimate power scalar coefficient of CrCL on CL ####

mod_cef3 <- PM_model$new(
  # ---------------- Priors / bounds ----------------
  pri = list(
    V1     = ab(0.1,  30),    # L (typical central, prior wide; scaled by weight in eqn)
    K12    = ab(  0,  10),    # 1/h
    K21    = ab(  0,  10),    # 1/h
    CL1    = ab(0.1,  12),    # L/h (renal scale at CRCL=120)
    S_eff  = ab(  0,  2),     # dimensionless sieving scaler for effluent clearance
    S_post = ab(  0,  2),     # dimensionless scaler for post-filter clearance
    theta3 = ab(  0,  3)      # CRCL power scaler
  ),
  
  # ---------------- Covariates (names must match your dataset column headers) -----
  cov = list(
    interval_volume = interp("none"),  # L
    age             = interp(),        # years
    male            = interp("none"),  # 0/1
    ht              = interp(),        # cm or m
    wt              = interp(),        # kg
    scr             = interp(),        # mg/dL
    crcl            = interp(),        # mL/min
    bsa             = interp(),        # m^2
    crcl_bsa        = interp(),        # mL/min/1.73 m^2
    ecmo            = interp("none"),  # 0/1
    hd              = interp("none"),  # 0/1
    crrt            = interp("none"),  # 0/1
    cvvh            = interp("none"),  # 0/1
    cvvhd           = interp("none"),  # 0/1
    cvvhdf          = interp("none"),  # 0/1
    flow            = interp("none"),  # mL/h
    bfr             = interp("none"),  # mL/min (blood flow rate)
    bag_reset       = interp("none")   # 0/1
  ),
  sec = function(){
    CL2 = 7.2                          # fixed HD clearance
    
    # Weight-scaled central volume
    V_WT = V1 * (wt/70)                # WT scaled Vd
    
    # Systemic (renal ± HD)
    CL_R   = CL1 * (crcl/120)**theta3                               # L/h
    CL_HD  = CL2 * hd                                               # L/h
    CL_sys = CL_R + CL_HD                                           # L/h
    
    # CRRT clearance term from effluent flow (mL/h -> L/h)
    Qe      = (flow/1000) * crrt      # L/h; 0 when CRRT=0
    CL_CRRT = S_eff * Qe              # L/h
    
    # Total clearance and elimination from central
    CLT = CL_sys + CL_CRRT            # L/h
    Ke  = CLT / V_WT                  # 1/h
    
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    
    # mg/L pre-filter
    Cpre  = X[1]/V_WT
    # mg/L post-filter
    Cpost = Cpre * S_post * crrt
    
    Aeff = 0
    Ceff = 0
    
    # Guard against divide-by-zero
    if (X[4] > 1e-6) {
      Ceff = X[3] / X[4]
    }
    
    if (X[3] > 1e-6) {
      Aeff = X[3]
    }
  },
  # ---------------- Differential equations ----------------
  eqn = function(){
    # Two compartment primary model structure
    dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2]
    dX[2] =                  K12*X[1] - K21*X[2]
    # Bag accumulation since last dump
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    dX[3] = Qe * S_eff * (X[1]/V_WT) * Iacc  - kdump * X[3] * (1 - Iacc)
    dX[4] = Qe * Iacc                        - kdump * X[4] * (1 - Iacc)
  },
  
  # ---------------- Outputs ----------------
  out = function(){
    Y[1]  = Cpre       # DVID=1: pre/regular plasma in mg/L
    Y[2]  = Cpost      # DVID=3: post-filter plasma (0 if CRRT off)
    Y[3]  = Ceff       # DVID=4: effluent concentration (0 if not collecting)
    Y[4]  = Aeff       # DVID=5: effluent cumulative amount (mg)
  },
  
  # ---------------- Error models ----------------
  # Proportional with a floor per-assay.
  err = list(
    proportional(2, c(1, 0.15, 0, 0)),   # Y1 (Cpre)
    proportional(2, c(1, 0.15, 0, 0)),   # Y2 (Cpost)
    proportional(2, c(1, 0.15, 0, 0)),   # Y3 (Ceff_out)
    proportional(2, c(1, 0.15, 0, 0))    # Y4 (Aeff)
  )
)

run6 <- mod_cef3$fit(data = dat3,
                     cycles = 500,
                     path = PATHS$runs,
                     run = 6,
                     points = 250,
                     overwrite=TRUE)

run6 <- PM_load(6, path = PATHS$runs)

# Run 7, Updating with fixed allometric coefficient on CL (WT/70)**0.75 and effect of CrCL/120 on CL ####

mod_cef4 <- PM_model$new(
  # ---------------- Priors / bounds ----------------
  pri = list(
    V1     = ab(0.1,  30),    # L (typical central, prior wide; scaled by weight in eqn)
    K12    = ab(  0,  10),    # 1/h
    K21    = ab(  0,  10),    # 1/h
    CL1    = ab(0.1,  12),    # L/h (renal scale at CRCL=120)
    S_eff  = ab(  0,  2),     # dimensionless sieving scaler for effluent clearance
    S_post = ab(  0,  2)      # dimensionless scaler for post-filter clearance
  ),
  
  # ---------------- Covariates (names must match your dataset column headers) -----
  cov = list(
    interval_volume = interp("none"),  # L
    age             = interp(),        # years
    male            = interp("none"),  # 0/1
    ht              = interp(),        # cm or m
    wt              = interp(),        # kg
    scr             = interp(),        # mg/dL
    crcl            = interp(),        # mL/min
    bsa             = interp(),        # m^2
    crcl_bsa        = interp(),        # mL/min/1.73 m^2
    ecmo            = interp("none"),  # 0/1
    hd              = interp("none"),  # 0/1
    crrt            = interp("none"),  # 0/1
    cvvh            = interp("none"),  # 0/1
    cvvhd           = interp("none"),  # 0/1
    cvvhdf          = interp("none"),  # 0/1
    flow            = interp("none"),  # mL/h
    bfr             = interp("none"),  # mL/min (blood flow rate)
    bag_reset       = interp("none")   # 0/1
  ),
  sec = function(){
    CL2 = 7.2                          # fixed HD clearance
    theta2 = 0.75                      # fixed allometric scaling factor
    # Weight-scaled central volume
    V_WT = V1 * (wt/70)                # WT scaled Vd
    
    # Systemic (renal ± HD)
    CL_R   = CL1 * (crcl/120) * (wt/70)**theta2                     # L/h
    CL_HD  = CL2 * hd                                               # L/h
    CL_sys = CL_R + CL_HD                                           # L/h
    
    # CRRT clearance term from effluent flow (mL/h -> L/h)
    Qe      = (flow/1000) * crrt      # L/h; 0 when CRRT=0
    CL_CRRT = S_eff * Qe              # L/h
    
    # Total clearance and elimination from central
    CLT = CL_sys + CL_CRRT            # L/h
    Ke  = CLT / V_WT                  # 1/h
    
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    
    # mg/L pre-filter
    Cpre  = X[1]/V_WT
    # mg/L post-filter
    Cpost = Cpre * S_post * crrt
    
    Aeff = 0
    Ceff = 0
    
    # Guard against divide-by-zero
    if (X[4] > 1e-6) {
      Ceff = X[3] / X[4]
    }
    
    if (X[3] > 1e-6) {
      Aeff = X[3]
    }
  },
  # ---------------- Differential equations ----------------
  eqn = function(){
    # Two compartment primary model structure
    dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2]
    dX[2] =                  K12*X[1] - K21*X[2]
    # Bag accumulation since last dump
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    dX[3] = Qe * S_eff * (X[1]/V_WT) * Iacc  - kdump * X[3] * (1 - Iacc)
    dX[4] = Qe * Iacc                        - kdump * X[4] * (1 - Iacc)
  },
  
  # ---------------- Outputs ----------------
  out = function(){
    Y[1]  = Cpre       # DVID=1: pre/regular plasma in mg/L
    Y[2]  = Cpost      # DVID=3: post-filter plasma (0 if CRRT off)
    Y[3]  = Ceff       # DVID=4: effluent concentration (0 if not collecting)
    Y[4]  = Aeff       # DVID=5: effluent cumulative amount (mg)
  },
  
  # ---------------- Error models ----------------
  # Proportional with a floor per-assay.
  err = list(
    proportional(2, c(1, 0.15, 0, 0)),   # Y1 (Cpre)
    proportional(2, c(1, 0.15, 0, 0)),   # Y2 (Cpost)
    proportional(2, c(1, 0.15, 0, 0)),   # Y3 (Ceff_out)
    proportional(2, c(1, 0.15, 0, 0))    # Y4 (Aeff)
  )
)

run7 <- mod_cef4$fit(data = dat3,
                     cycles = 500,
                     path = PATHS$runs,
                     run = 7,
                     points = 250,
                     overwrite=TRUE)

run7 <- PM_load(7, path = PATHS$runs)

# Run 8, Updating fixed allometric coefficient: (wt/70)**0.75 and power scalar (crcl/120) effect on CL ####

mod_cef5 <- PM_model$new(
  # ---------------- Priors / bounds ----------------
  pri = list(
    V1     = ab(0.1,  30),      # L (typical central, prior wide; scaled by weight in eqn)
    K12    = ab(  0,  10),      # 1/h
    K21    = ab(  0,  10),      # 1/h
    CL1    = ab(0.1,  12),      # L/h (renal scale at CRCL=120)
    S_eff  = ab(  0,  2),       # dimensionless sieving scaler for effluent clearance
    S_post = ab(  0,  2),       # dimensionless scaler for post-filter clearance
    theta3 = ab(  0,  3)        # CRCL power scaler
  ),
  
  # ---------------- Covariates (names must match your dataset column headers) -----
  cov = list(
    interval_volume = interp("none"),  # L
    age             = interp(),        # years
    male            = interp("none"),  # 0/1
    ht              = interp(),        # cm or m
    wt              = interp(),        # kg
    scr             = interp(),        # mg/dL
    crcl            = interp(),        # mL/min
    bsa             = interp(),        # m^2
    crcl_bsa        = interp(),        # mL/min/1.73 m^2
    ecmo            = interp("none"),  # 0/1
    hd              = interp("none"),  # 0/1
    crrt            = interp("none"),  # 0/1
    cvvh            = interp("none"),  # 0/1
    cvvhd           = interp("none"),  # 0/1
    cvvhdf          = interp("none"),  # 0/1
    flow            = interp("none"),  # mL/h
    bfr             = interp("none"),  # mL/min (blood flow rate)
    bag_reset       = interp("none")   # 0/1
  ),
  sec = function(){
    CL2 = 7.2                          # fixed HD clearance
    theta2 = 0.75                      # fixed allometric scaling factor
    # Weight-scaled central volume
    V_WT = V1 * (wt/70)                # WT scaled Vd
    
    # Systemic (renal ± HD)
    CL_R   = CL1 * (crcl/120)**theta3  * (wt/70)**theta2            # L/h
    CL_HD  = CL2 * hd                                               # L/h
    CL_sys = CL_R + CL_HD                                           # L/h
    
    # CRRT clearance term from effluent flow (mL/h -> L/h)
    Qe      = (flow/1000) * crrt      # L/h; 0 when CRRT=0
    CL_CRRT = S_eff * Qe              # L/h
    
    # Total clearance and elimination from central
    CLT = CL_sys + CL_CRRT            # L/h
    Ke  = CLT / V_WT                  # 1/h
    
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    
    # mg/L pre-filter
    Cpre  = X[1]/V_WT
    # mg/L post-filter
    Cpost = Cpre * S_post * crrt
    
    Aeff = 0
    Ceff = 0
    
    # Guard against divide-by-zero
    if (X[4] > 1e-6) {
      Ceff = X[3] / X[4]
    }
    
    if (X[3] > 1e-6) {
      Aeff = X[3]
    }
  },
  # ---------------- Differential equations ----------------
  eqn = function(){
    # Two compartment primary model structure
    dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2]
    dX[2] =                  K12*X[1] - K21*X[2]
    # Bag accumulation since last dump
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    dX[3] = Qe * S_eff * (X[1]/V_WT) * Iacc  - kdump * X[3] * (1 - Iacc)
    dX[4] = Qe * Iacc                        - kdump * X[4] * (1 - Iacc)
  },
  
  # ---------------- Outputs ----------------
  out = function(){
    Y[1]  = Cpre       # DVID=1: pre/regular plasma in mg/L
    Y[2]  = Cpost      # DVID=3: post-filter plasma (0 if CRRT off)
    Y[3]  = Ceff       # DVID=4: effluent concentration (0 if not collecting)
    Y[4]  = Aeff       # DVID=5: effluent cumulative amount (mg)
  },
  
  # ---------------- Error models ----------------
  # Proportional with a floor per-assay.
  err = list(
    proportional(2, c(1, 0.15, 0, 0)),   # Y1 (Cpre)
    proportional(2, c(1, 0.15, 0, 0)),   # Y2 (Cpost)
    proportional(2, c(1, 0.15, 0, 0)),   # Y3 (Ceff_out)
    proportional(2, c(1, 0.15, 0, 0))    # Y4 (Aeff)
  )
)

run8 <- mod_cef5$fit(data = dat3,
                     cycles = 500,
                     path = PATHS$runs,
                     run = 8,
                     points = 250,
                     overwrite=TRUE)

run8 <- PM_load(8, path = PATHS$runs)

# Run once per session before exporting figures
reticulate::py_run_string('import sys')

library(plotly)

f <- list(size = 20)

a <- list(
  text = "A. Pre-filter",
  font=f,
  xref = "paper",
  yref = "paper",
  yanchor = "bottom",
  xanchor = "center",
  align = "right",
  x = 0.1,
  y = 1,
  showarrow = FALSE
)

b <- list(
  text = "B. Post-filter",
  font=f,
  xref = "paper",
  yref = "paper",
  yanchor = "bottom",
  xanchor = "center",
  align = "right",
  x = 0.1,
  y = 1,
  showarrow = FALSE
)

c <- list(
  text = "C. Effluent (mg/L)",
  font=f,
  xref = "paper",
  yref = "paper",
  yanchor = "bottom",
  xanchor = "center",
  align = "right",
  x = 0.1,
  y = 1,
  showarrow = FALSE
)

d <- list(
  text = "D. Effluent (mg)",
  font=f,
  xref = "paper",
  yref = "paper",
  yanchor = "bottom",
  xanchor = "center",
  align = "right",
  x = 0.1,
  y = 1,
  showarrow = FALSE
)


op1<-run8$op$plot(pred.type="post",outeq=1,
                  line=list(lm=list(T,color="dodgerblue")),
                  marker=list(color="goldenrod"),
                  stats = list(x = 1, y = 0.5, font = list(size = 10)))

op2<-run8$op$plot(pred.type="post",outeq=2,
                  line=list(lm=list(T,color="dodgerblue")),
                  marker=list(color="goldenrod"),
                  stats = list(x = 1, y = 0.45, font = list(size = 10)))

op3<-run8$op$plot(pred.type="post",outeq=3,
                  line=list(lm=list(T,color="dodgerblue")),
                  marker=list(color="goldenrod"),
                  stats = list(x = 1, y = 0.45, font = list(size = 10)))

op4<-run8$op$plot(pred.type="post",outeq=4,
                  line=list(lm=list(T,color="dodgerblue")),
                  marker=list(color="goldenrod"),
                  stats = list(x = 1, y = 0.15, font = list(size = 10)))

op1 <- op1 %>% layout(annotations = a)
op2 <- op2 %>% layout(annotations = b)
op3 <- op3 %>% layout(annotations = c)
op4 <- op4 %>% layout(annotations = d)

# Write to Figures folder
{
  out_file <- "FEP DEV COHORT POST FDA RENAL Plot_test.svg"
  out_path <- file.path(PATHS$fig, out_file)
  
  sub_plot(op1, op2, op3, op4,
           nrows = 2, margin = 0.05,
           titleX = TRUE, titleY = TRUE,
           shareX = TRUE, shareY = TRUE) %>%
    export_plotly(out_file, width = 2.5 * 300, height = 3.5 * 300)
  
  if (file.exists(out_file)) {
    dir.create(PATHS$fig, recursive = TRUE, showWarnings = FALSE)
    file.rename(out_file, out_path)
  } else {
    warning("Expected export not found: ", out_file)
  }
  }

# Run 9, crcl theta, fixing allometric wt on CL scalar, estimating saturation of crrt filter ####
mod_cef6 <- PM_model$new(
  # ---------------- Priors / bounds ----------------
  pri = list(
    V1     = ab(0.1,  30),   # L (typical central, prior wide; scaled by weight in eqn)
    K12    = ab(  0,  10),   # 1/h
    K21    = ab(  0,  10),   # 1/h
    CL1    = ab(0.1,  12),   # L/h (renal scale at CRCL=120)
    S_eff  = ab(  0,  2),    # dimensionless sieving scaler for effluent clearance
    S_post = ab(  0,  2),    # dimensionless scaler for post-filter clearance
    theta3 = ab(  0,  3),    # CRCL power scaler
    sat    = ab(0.5,  1)     # saturation of crrt sieving capacity
  ),
  
  # ---------------- Covariates (names must match your dataset column headers) -----
  cov = list(
    interval_volume = interp("none"),  # L
    age             = interp(),        # years
    male            = interp("none"),  # 0/1
    ht              = interp(),        # cm or m
    wt              = interp(),        # kg
    scr             = interp(),        # mg/dL
    crcl            = interp(),        # mL/min
    bsa             = interp(),        # m^2
    crcl_bsa        = interp(),        # mL/min/1.73 m^2
    ecmo            = interp("none"),  # 0/1
    hd              = interp("none"),  # 0/1
    crrt            = interp("none"),  # 0/1
    cvvh            = interp("none"),  # 0/1
    cvvhd           = interp("none"),  # 0/1
    cvvhdf          = interp("none"),  # 0/1
    flow            = interp("none"),  # mL/h
    bfr             = interp("none"),  # mL/min (blood flow rate)
    bag_reset       = interp("none")   # 0/1
  ),
  sec = function(){
    CL2 = 7.2                          # fixed HD clearance
    theta2 = 0.75                      # fixed allometric scaling factor
    # Weight-scaled central volume
    V_WT = V1 * (wt/70)                # WT scaled Vd
    # S_eff  = 0.81                    # not fixing for now
    # S_post = 0.962                   # not fixing for now
    
    
    # Systemic (renal ± HD)
    CL_R   = CL1 * (crcl/120)**theta3  * (wt/70)**theta2            # L/h
    CL_HD  = CL2 * hd                                               # L/h
    CL_sys = CL_R + CL_HD                                           # L/h
    
    # CRRT clearance term from effluent flow (mL/h -> L/h)
    Qe      = (flow/1000) * crrt      # L/h; 0 when CRRT=0
    Qe_eff  = (Qe/2.6)**sat           # saturation of filter by power rule
    CL_CRRT = S_eff * Qe_eff          # L/h
    
    # Total clearance and elimination from central
    CLT = CL_sys + CL_CRRT            # L/h
    Ke  = CLT / V_WT                  # 1/h
    
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    
    # mg/L pre-filter
    Cpre  = X[1]/V_WT
    # mg/L post-filter
    Cpost = Cpre * S_post * crrt
    
    Aeff = 0
    Ceff = 0
    
    # Guard against divide-by-zero
    if (X[4] > 1e-6) {
      Ceff = X[3] / X[4]
    }
    
    if (X[3] > 1e-6) {
      Aeff = X[3]
    }
  },
  # ---------------- Differential equations ----------------
  eqn = function(){
    # Two compartment primary model structure
    dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2]
    dX[2] =                  K12*X[1] - K21*X[2]
    # Bag accumulation since last dump
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    dX[3] = Qe_eff * S_eff * (X[1]/V_WT) * Iacc  - kdump * X[3] * (1 - Iacc)
    dX[4] = Qe_eff * Iacc                        - kdump * X[4] * (1 - Iacc)
  },
  
  # ---------------- Outputs ----------------
  out = function(){
    Y[1]  = Cpre       # DVID=1: pre/regular plasma in mg/L
    Y[2]  = Cpost      # DVID=3: post-filter plasma (0 if CRRT off)
    Y[3]  = Ceff       # DVID=4: effluent concentration (0 if not collecting)
    Y[4]  = Aeff       # DVID=5: effluent cumulative amount (mg)
  },
  
  # ---------------- Error models ----------------
  # Proportional with a floor per-assay.
  err = list(
    proportional(2, c(1, 0.15, 0, 0)),   # Y1 (Cpre)
    proportional(2, c(1, 0.15, 0, 0)),   # Y2 (Cpost)
    proportional(2, c(1, 0.15, 0, 0)),   # Y3 (Ceff_out)
    proportional(2, c(1, 0.15, 0, 0))    # Y4 (Aeff)
  )
)

run9 <- mod_cef6$fit(data = dat3,
                     cycles = 500,
                     path = PATHS$runs,
                     run = 9,
                     points = 250,
                     overwrite=TRUE)

run9 <- PM_load(9, path = PATHS$runs)

# =============================================================================
# 3. Compare rich CRRT plasma/CRRT structural candidates ####
# =============================================================================

# compare runs
PM_compare(run4, run5, run6, run7)
PM_compare(run8, run9)

cmp_47 <- PM_compare(run4, run5, run6, run7)
cmp_89 <- PM_compare(run8, run9)

compare_tbl <- tibble(
  run = c(4, 5, 6, 7, 8, 9),
  nvar = c(cmp_47$nvar, cmp_89$nvar),
  ll2  = c(cmp_47$`-2*LL`, cmp_89$`-2*LL`),
  aic  = c(cmp_47$aic, cmp_89$aic),
  notes = c(
    "base: CL1*(crcl/120) + S_eff*Qe; V1*(wt/70)",
    "add BSA power scalar to renal CL: CL1*(crcl/120)*(bsa/1.73)^theta1",
    "estimate CRCL power scalar: CL1*(crcl/120)^theta3",
    "fix WT allometry (0.75) on renal CL: CL1*(crcl/120)*(wt/70)^0.75",
    "combine: CL1*(crcl/120)^theta3*(wt/70)^0.75",
    "add effluent flow saturation: CL_CRRT via (Qe/2.6)^sat"
  )
) %>%
  mutate(
    delta_ofv_run4 = ll2 - ll2[run == 4],
    delta_aic = aic - min(aic)
  ) %>%
  arrange(run)

compare_tbl

library(knitr)
library(kableExtra)

kable(
  compare_tbl,
  digits = 1,
  caption = "Model comparison across CRRT in development",
  align = "c"
) %>%
  kable_styling(full_width = FALSE, position = "center")

compare_tbl %>%
  kable(
    digits  = 1,
    caption = "Model comparison across CRRT structural variants",
    align   = "c",
    format  = "html"
  ) %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  save_kable(file.path(PATHS$tab, "Table_S1_Development_model_comparison.html"))

write.csv(compare_tbl, file.path(PATHS$tab, "Table_S1_model_comparison.csv"), row.names = FALSE)

# =============================================================================
# 3A. Diagnostic posterior plots in the rich CRRT cohort ####
# =============================================================================

# generate comparement posterior plots in the rich sampling cohort ####
toggle_tad_plot <- function(run, tad = 0, tau = NULL, legend_rows = 2) {
  # deps
  library(dplyr); library(ggplot2)
  
  # ---------- data prep ----------
  # dose times
  dose_tbl <- run$data$data %>%
    filter(evid == 1) %>%
    select(id, dose_time = time) %>%
    arrange(id, dose_time)
  
  # observed / predicted
  obs <- run$data$data %>%
    filter(evid == 0, out != -99, !is.na(outeq)) %>%
    transmute(id = as.character(id), time, outeq, out)
  
  pred <- run$post$data %>%
    filter(icen == "median") %>%
    transmute(id = as.character(id), time, outeq, out = pred)
  
  # helper: TAD + cycle + local_tau (time to next dose for that cycle)
  add_tad_cycle <- function(df) {
    df %>%
      arrange(id, time) %>%
      group_by(id) %>%
      mutate({
        dt  <- dose_tbl$dose_time[dose_tbl$id == first(id)]
        idx <- findInterval(time, dt)                 # most recent dose index
        n   <- length(dt)
        nxt <- ifelse(idx > 0 & idx < n, idx + 1L, NA_integer_)
        tibble(
          cycle     = ifelse(idx == 0, NA_integer_, idx),
          TAD       = ifelse(idx == 0, NA_real_,  time - dt[idx]),
          local_tau = ifelse(idx == 0 | is.na(nxt), NA_real_, dt[nxt] - dt[idx])
        )
      }) %>% ungroup()
  }
  
  if (tad == 1) {
    # global display tau if not provided
    if (is.null(tau)) {
      tau <- dose_tbl %>%
        group_by(id) %>%
        summarize(tau = median(diff(dose_time), na.rm = TRUE), .groups = "drop") %>%
        summarize(tau = median(tau, na.rm = TRUE), .groups = "drop") %>%
        pull(tau)
    }
    
    obs <- add_tad_cycle(obs)  %>%
      mutate(eff_tau = ifelse(is.na(local_tau), tau, pmin(local_tau, tau))) %>%
      filter(!is.na(TAD), TAD <= eff_tau) %>%
      mutate(x = TAD, grp = interaction(id, cycle, drop = TRUE))
    
    pred <- add_tad_cycle(pred) %>%
      mutate(eff_tau = ifelse(is.na(local_tau), tau, pmin(local_tau, tau))) %>%
      filter(!is.na(TAD), TAD <= eff_tau) %>%
      mutate(x = TAD, grp = interaction(id, cycle, drop = TRUE))
    
    ggplot() +
      geom_line(data = pred, aes(x = x, y = out, group = grp, color = id),
                linewidth = 0.6, alpha = 0.9) +
      geom_point(data = obs, aes(x = x, y = out, color = id),
                 shape = 21, fill = NA, size = 2, stroke = 0.9) +
      facet_wrap(~ outeq, scales = "free_y") +
      scale_x_continuous(limits = c(0, tau),
                         breaks = scales::breaks_extended(5),
                         expand = expansion(mult = c(0.01, 0.03))) +
      theme_bw(base_size = 14) +
      labs(x = sprintf("Time After Dose (hr) — within 0–%.0f hr", tau),
           y = "Concentration (mg/L)",
           color = "Subject ID",
           title = "Observed vs Predicted by Output Equation (TAD)") +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 10),
        legend.text  = element_text(size = 9),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width  = unit(0.8, "cm"),
        legend.box.margin = margin(t = -5)
      ) +
      guides(color = guide_legend(nrow = legend_rows))
    
  } else {
    # clock time
    obs  <- obs  %>% mutate(x = time, grp = id)
    pred <- pred %>% mutate(x = time, grp = id)
    
    ggplot() +
      geom_line(data = pred, aes(x = x, y = out, group = grp, color = id),
                linewidth = 0.7, alpha = 0.8) +
      geom_point(data = obs, aes(x = x, y = out, color = id),
                 shape = 21, fill = NA, size = 2, stroke = 1) +
      facet_wrap(~ outeq, scales = "free_y") +
      theme_bw(base_size = 14) +
      labs(x = "Time (hr)", y = "Concentration (mg/L)",
           color = "Subject ID",
           title = "Observed vs Predicted by Output Equation") +
      theme(
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.title = element_text(size = 10),
        legend.text  = element_text(size = 9),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width  = unit(0.8, "cm"),
        legend.box.margin = margin(t = -5)
      ) +
      guides(color = guide_legend(nrow = legend_rows))
  }
}
# EXAMPLES:
toggle_tad_plot(run6, tad = 0)          # clock time
toggle_tad_plot(run6, tad = 1)          # TAD, auto τ from dosing
toggle_tad_plot(run6, tad = 1, tau = 6) # TAD with fixed display τ

# =============================================================================
# 3B. Exploratory CRRT clearance diagnostic ####
# =============================================================================

# Generate OP plots for observed vs predicted CRRT CL using posteriors

# -----------------------------------------------------------
# Run 9: predicted vs observed CRRT clearance
#
# Predicted (model-based):
#   Qe = flow/1000 * crrt
#   CL_CRRT_pred = Qe * s_eff
#
# Observed (Schetz hemofiltration definition):
#   S_obs = Ceff / Cpre
#   CL_CRRT_obs = Qe * S_obs
#
# Optional post-filter-corrected observed version:
#   Cpre_corr = Cpost / s_post
#   CL_CRRT_obs_postcorr = Qe * (Ceff / Cpre_corr)
# ----------------------------------------------------------

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# start from aligned OP data
op_long <- run6$op$data %>%
  as_tibble() %>%
  filter(
    pred.type == "post",
    icen == "median",
    outeq %in% c(1, 2, 3)
  ) %>%
  transmute(
    id,
    time,
    outeq,
    obs  = as.numeric(obs),
    pred = as.numeric(pred)
  )

# posterior median machine terms
cov_med <- run6$cov$data %>%
  as_tibble() %>%
  filter(icen == "median") %>%
  transmute(
    id,
    time,
    flow,
    crrt,
    s_eff,
    s_post,
    qeff_l_h = (flow / 1000) * crrt,
    cl_crrt_pred_l_h = qeff_l_h * s_eff
  ) %>%
  distinct(id, time, .keep_all = TRUE)

# split by assay
pre_dat <- op_long %>%
  filter(outeq == 1, !is.na(obs), obs > 0) %>%
  transmute(
    id,
    time_pre = time,
    obs_pre  = obs,
    pred_pre = pred
  )

post_dat <- op_long %>%
  filter(outeq == 2, !is.na(obs), obs > 0) %>%
  transmute(
    id,
    time_post = time,
    obs_post  = obs,
    pred_post = pred
  )

eff_dat <- op_long %>%
  filter(outeq == 3, !is.na(obs), obs > 0) %>%
  transmute(
    id,
    time_eff = time,
    obs_eff  = obs,
    pred_eff = pred
  )

# attach machine terms to effluent times
eff_base <- eff_dat %>%
  left_join(
    cov_med %>% rename(time_eff = time),
    by = c("id", "time_eff")
  ) %>%
  filter(crrt == 1, qeff_l_h > 0)

# helper: nearest-time join within subject
nearest_join_one <- function(eff_row, ref_df, time_col, value_cols, tol = 6) {
  sub <- ref_df %>% filter(id == eff_row$id)
  if (nrow(sub) == 0) return(as_tibble(setNames(rep(list(NA), length(value_cols) + 1),
                                                c("dt", value_cols))))
  
  dt <- abs(sub[[time_col]] - eff_row$time_eff)
  i <- which.min(dt)
  
  if (length(i) == 0 || is.infinite(dt[i]) || dt[i] > tol) {
    return(as_tibble(setNames(rep(list(NA), length(value_cols) + 1),
                              c("dt", value_cols))))
  }
  
  out <- c(dt = dt[i], as.list(sub[i, value_cols, drop = TRUE]))
  as_tibble(out)
}

# pair nearest pre and post to each effluent row
paired <- map_dfr(seq_len(nrow(eff_base)), function(i) {
  erow <- eff_base[i, ]
  
  pre_match <- nearest_join_one(
    eff_row   = erow,
    ref_df    = pre_dat,
    time_col  = "time_pre",
    value_cols = c("time_pre", "obs_pre", "pred_pre"),
    tol = 6
  )
  
  post_match <- nearest_join_one(
    eff_row   = erow,
    ref_df    = post_dat,
    time_col  = "time_post",
    value_cols = c("time_post", "obs_post", "pred_post"),
    tol = 6
  )
  
  bind_cols(erow, pre_match, post_match)
})

# calculate observed and predicted CRRT clearance
cl_pair <- paired %>%
  mutate(
    s_obs_schetz = case_when(
      !is.na(obs_pre) & obs_pre > 0 ~ obs_eff / obs_pre,
      TRUE ~ NA_real_
    ),
    cl_crrt_obs_l_h = qeff_l_h * s_obs_schetz,
    
    s_pred_schetz = case_when(
      !is.na(pred_pre) & pred_pre > 0 ~ pred_eff / pred_pre,
      TRUE ~ NA_real_
    ),
    cl_crrt_pred_from_op_l_h = qeff_l_h * s_pred_schetz,
    
    pre_from_post_obs = case_when(
      !is.na(obs_post) & !is.na(s_post) & s_post > 0 ~ obs_post / s_post,
      TRUE ~ NA_real_
    ),
    s_obs_postcorr = case_when(
      !is.na(pre_from_post_obs) & pre_from_post_obs > 0 ~ obs_eff / pre_from_post_obs,
      TRUE ~ NA_real_
    ),
    cl_crrt_obs_postcorr_l_h = qeff_l_h * s_obs_postcorr
  )

# inspect
cl_pair %>%
  select(
    id, time_eff, flow, qeff_l_h,
    s_eff, s_post,
    time_pre, obs_pre, pred_pre,
    time_post, obs_post, pred_post,
    obs_eff, pred_eff,
    cl_crrt_pred_l_h,
    cl_crrt_pred_from_op_l_h,
    cl_crrt_obs_l_h,
    cl_crrt_obs_postcorr_l_h
  ) %>%
  print(n = 30)

# plot
plot_df <- cl_pair %>%
  filter(!is.na(cl_crrt_obs_l_h), !is.na(cl_crrt_pred_l_h))

lim <- max(plot_df$cl_crrt_pred_l_h, plot_df$cl_crrt_obs_l_h, na.rm = TRUE)+1

ggplot(plot_df, aes(x = cl_crrt_pred_l_h, y = cl_crrt_obs_l_h)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = T) +
  coord_equal(xlim = c(0, lim), ylim = c(0, lim), expand = FALSE) +
  labs(
    title = "Model-predicted CRRT clearance",
    subtitle = "Observed CL using nearest paired pre-filter and effluent samples within subject",
    x = "Predicted CRRT clearance (L/h)",
    y = "Observed CRRT clearance (L/h)"
  ) +
  theme_bw(base_size = 12)

# =============================================================================
# 4. Build combined plasma + ELF dataset and create train/test split ####
# =============================================================================

# Link CRRT to sparse plasma pre-filter and ELF samples

## --- Source ---
library(readxl)
src_path  <- "src/FEP_comb_ELF_old_ML_combined_FDA.xlsx"
src_sheet <- "FEP_comb_ELF_old"

dat <- readxl::read_excel(src_path, sheet = src_sheet)

## --- Pull the single ELF row for subject 2103; force-include by clearing EXCLUDE_DV ---
row_special <- dat %>%
  mutate(
    ID   = suppressWarnings(as.numeric(ID)),
    DVID = suppressWarnings(as.numeric(DVID))
  ) %>%
  filter(ID == 2103, DVID == 2) %>%
  mutate(EXCLUDE_DV = 0)     # make sure it won't be dropped

## --- Build main dataset (drop EXCLUDE_DV==1, keep RRT==1) ---
dat4 <- dat %>%
  mutate(
    ID   = suppressWarnings(as.numeric(ID)),
    DVID = suppressWarnings(as.numeric(DVID)),
    RRT  = suppressWarnings(as.numeric(RRT)),
    EVID = suppressWarnings(as.numeric(EVID))
  ) %>%
  filter(RRT == 1, EXCLUDE_DV != 1) %>%
  # Standard names for Pmetrics
  rename(
    DOSE  = AMT,
    DUR   = TINF,
    INPUT = ADM,
    OUTEQ = DVID,
    OUT   = DV,
    LOQ   = LIMIT
  ) %>%
  # Coerce numeric where Excel might have text
  mutate(
    OUTEQ = suppressWarnings(as.numeric(OUTEQ)),
    LOQ   = suppressWarnings(as.numeric(LOQ)),
    OUT   = suppressWarnings(as.numeric(OUT)),
    TIME  = suppressWarnings(as.numeric(TIME))
  ) %>%
  # Map original DVID -> unified OUTEQ (ELF=5)
  # 1->1 (Cpre), 3->2 (Cpost), 4->3 (Ceff conc), 5->4 (Aeff), 2->5 (ELF)
  mutate(
    OUTEQ = case_when(
      OUTEQ == 1 ~ 1,
      OUTEQ == 3 ~ 2,
      OUTEQ == 4 ~ 3,
      OUTEQ == 5 ~ 4,
      OUTEQ == 2 ~ 5,
      TRUE ~ NA_real_
    ),
    # Fill LLOQ on observed rows
    LOQ = case_when(
      EVID == 1 ~ NA_real_,
      EVID == 0 & is.na(LOQ) ~ 0.5,
      TRUE ~ LOQ
    ),
    # Missing obs → -99 for observations
    OUT = case_when(
      EVID == 0 & is.na(OUT) ~ -99,
      TRUE ~ OUT
    ),
    # Set BAG_RESET = 0 on dump rows (INPUT==2), keep existing otherwise
    BAG_RESET = case_when(
      INPUT == 2 ~ 0,
      TRUE ~ BAG_RESET
    )
  )

## --- Append the special row (rename + recode to OUTEQ=5) ---
row_special_std <- row_special %>%
  rename(
    DOSE  = AMT,
    DUR   = TINF,
    INPUT = ADM,
    OUTEQ = DVID,
    OUT   = DV,
    LOQ   = LIMIT
  ) %>%
  mutate(
    OUTEQ = 5,  # ELF unified code
    LOQ   = suppressWarnings(as.numeric(LOQ)),
    OUT   = suppressWarnings(as.numeric(OUT)),
    TIME  = suppressWarnings(as.numeric(TIME))
  )

## Prevent duplicates if it somehow already exists post-mapping
dat4 <- dat4 %>%
  anti_join(row_special_std %>% select(ID, TIME, OUTEQ), by = c("ID","TIME","OUTEQ")) %>%
  bind_rows(row_special_std)


## --- Final tidy: drop clutter, order, lowercase, NA to "."
dat4 <- dat4 %>%
  select(
    ID, EVID, TIME, DOSE, DUR, ADDL, II, INPUT, OUT, OUTEQ,
    INTERVAL_VOLUME, AGE, MALE, HTIN, HT, WT, SCR, CRCL, BSA, CRCL_BSA,
    ECMO, HD, CRRT, CVVH, CVVHD, CVVHDF, FLOW, BFR, BAG_RESET,
    RX_include, no_RX_exclude, no_RRT, RX_ID, DATE, CENS, MDV, ALIQUOT_ID,
    everything()
  ) %>%
  select(
    -RX_ID, -DATE, -ALIQUOT_ID, -no_RRT, -no_RX_exclude,
    -MDV, -RX_include, -RRT, -EXCLUDE_DV,
    -PLA_urea, -BAL_urea, -BAL_PK, - LOQ
  ) %>%
  arrange(ID, TIME, desc(EVID)) %>%
  rename_with(tolower) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), ".", as.character(.))))

## --- Identify which IDs were used in dat3 (so we can stratify the 80/20 split) ---
dat3_ids <- dat %>%
  filter(EXCLUDE_DV != 1,
         RRT == 1,
         RX_include == 1) %>%
  distinct(ID) %>%
  pull()

## --- Stratified 80/20 split by group (in dat3 vs not) ---
set.seed(12345)

# Training IDs: 80% from each group
train_ids_in_dat3 <- dat4 %>%
  filter(id %in% dat3_ids) %>%
  distinct(id) %>%
  sample_frac(0.8, replace = FALSE)

train_ids_new <- dat4 %>%
  filter(!id %in% dat3_ids) %>%
  distinct(id) %>%
  sample_frac(0.8, replace = FALSE)

id_train <- bind_rows(train_ids_in_dat3, train_ids_new)

## --- Split datasets accordingly ---
dat4_train <- dat4 %>%
  semi_join(id_train, by = "id") %>%
  arrange(id, as.numeric(time))

dat4_test <- dat4 %>%
  anti_join(id_train, by = "id") %>%
  arrange(id, as.numeric(time))

## --- Write CSVs and instantiate PM_data objects ---
write.csv(dat4_train, "src/train2.csv", row.names = FALSE)
write.csv(dat4_test,  "src/test2.csv",  row.names = FALSE)

dat5 <- PM_data$new("src/train2.csv", loq = c(0,0,0,0,0))
dat6 <- PM_data$new("src/test2.csv",  loq = c(0,0,0,0,0))

# =============================================================================
# 5. Fit primary ELF-extended candidate models already present in the workflow ####
# =============================================================================

# Updating fixed allometric coefficient of  0.75 on wt/70 and theta on crcl/120 for  CL

mod_cef7 <- PM_model$new(
  # ---------------- Priors / bounds ----------------
  pri = list(
    V1     = ab(0.1,  30),     # L (typical central, prior wide; scaled by weight in eqn)
    V2     = ab(0.1, 100),     # L (ELF volume of distribution)
    K12    = ab(  0,  10),     # 1/h central to peripheral
    K21    = ab(  0,  10),     # 1/h peripheral to central
    K15    = ab(  0,  10),     # 1/h central to ELF
    K51    = ab(  0,  10),     # 1/h ELF  to central
    CL1    = ab(0.1,  12),     # L/h (renal scale at CRCL=120)
    S_eff  = ab(  0,  2),      # dimensionless sieving scaler for effluent clearance
    S_post = ab(  0,  2),      # dimensionless scaler for post-filter clearance
    theta3 = ab(  0,  3)       # CRCL power scaler
  ),
  
  # ---------------- Covariates (names must match your dataset column headers) -----
  cov = list(
    interval_volume = interp("none"),  # L
    age             = interp(),        # years
    male            = interp("none"),  # 0/1
    ht              = interp(),        # cm or m
    wt              = interp(),        # kg
    scr             = interp(),        # mg/dL
    crcl            = interp(),        # mL/min
    bsa             = interp(),        # m^2
    crcl_bsa        = interp(),        # mL/min/1.73 m^2
    ecmo            = interp("none"),  # 0/1
    hd              = interp("none"),  # 0/1
    crrt            = interp("none"),  # 0/1
    cvvh            = interp("none"),  # 0/1
    cvvhd           = interp("none"),  # 0/1
    cvvhdf          = interp("none"),  # 0/1
    flow            = interp("none"),  # mL/h
    bfr             = interp("none"),  # mL/min (blood flow rate)
    bag_reset       = interp("none")   # 0/1
  ),
  sec = function(){
    CL2 = 7.2                          # fixed HD clearance
    theta2 = 0.75                      # fixed allometric scaling factor
    # Weight-scaled central volume
    V_WT = V1 * (wt/70)                # WT scaled Vd
    
    # Systemic (renal ± HD)
    CL_R   = CL1 * (crcl/120)**theta3  * (wt/70)**theta2            # L/h
    CL_HD  = CL2 * hd                                               # L/h
    CL_sys = CL_R + CL_HD                                           # L/h
    
    # CRRT clearance term from effluent flow (mL/h -> L/h)
    Qe      = (flow/1000) * crrt      # L/h; 0 when CRRT=0
    CL_CRRT = S_eff * Qe              # L/h
    
    # Total clearance and elimination from central
    CLT = CL_sys + CL_CRRT            # L/h
    Ke  = CLT / V_WT                  # 1/h
    
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    
    # mg/L pre-filter
    Cpre  = X[1]/V_WT
    # mg/L post-filter
    Cpost = Cpre * S_post * crrt
    
    Aeff = 0
    Ceff = 0
    
    # Guard against divide-by-zero
    if (X[4] > 1e-6) {
      Ceff = X[3] / X[4]
    }
    
    if (X[3] > 1e-6) {
      Aeff = X[3]
    }
    
    Celf = X[5]/V2
  },
  # ---------------- Differential equations ----------------
  eqn = function(){
    # Three compartment primary model structure
    dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2] - K15*X[1] + K51*X[5]
    dX[2] =                  K12*X[1] - K21*X[2]
    # Bag accumulation since last dump
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    dX[3] = Qe * S_eff * (X[1]/V_WT) * Iacc  - kdump * X[3] * (1 - Iacc)
    dX[4] = Qe * Iacc                        - kdump * X[4] * (1 - Iacc)
    dX[5] =                  K15*X[1] - K51*X[5]
  },
  
  # ---------------- Outputs ----------------
  out = function(){
    Y[1]  = Cpre       # DVID=1: pre/regular plasma in mg/L
    Y[2]  = Cpost      # DVID=2: post-filter plasma (0 if CRRT off)
    Y[3]  = Ceff       # DVID=3: effluent concentration (0 if not collecting)
    Y[4]  = Aeff       # DVID=4: effluent cumulative amount (mg)
    Y[5]  = Celf       # DVID=5: ELF concentrations from urea corrected BAL method
  },
  
  # ---------------- Error models ----------------
  # Proportional with a floor per-assay.
  err = list(
    proportional(2, c(1, 0.15, 0, 0)),   # Y1 (Cpre)
    proportional(2, c(1, 0.15, 0, 0)),   # Y2 (Cpost)
    proportional(2, c(1, 0.15, 0, 0)),   # Y3 (Ceff_out)
    proportional(2, c(1, 0.15, 0, 0)),   # Y4 (Aeff)
    proportional(2, c(1, 0.15, 0, 0))    # Y5 (Celf)
  )
)

run10 <- mod_cef7$fit(data = dat5,
                      cycles = 1000,
                      path = PATHS$runs,
                      run = 10,
                      points = 300,
                      overwrite=TRUE)

run10 <- PM_load(10, path = PATHS$runs)

# model validation fit

run11 <- mod_cef7$fit(data = dat6,
                      cycles = 0,
                      path = PATHS$runs,
                      run = 11,
                      prior = 10,
                      overwrite=TRUE)

run11 <- PM_load(11, path = PATHS$runs)

mod_cef8 <- PM_model$new(
  # ---------------- Priors / bounds ----------------
  pri = list(
    V1     = ab(0.1,  30),     # L (typical central, prior wide; scaled by weight in eqn)
    V2     = ab(0.1, 100),     # L (ELF volume of distribution)
    K12    = ab(  0,  10),     # 1/h central to peripheral
    K21    = ab(  0,  10),     # 1/h peripheral to central
    K15    = ab(  0,  10),     # 1/h central to ELF
    K51    = ab(  0,  10),     # 1/h ELF  to central
    CL1    = ab(0.1,  12),     # L/h (renal scale at CRCL=120)
    S_eff  = ab(  0,  2),      # dimensionless sieving scaler for effluent clearance
    S_post = ab(  0,  2),      # dimensionless scaler for post-filter clearance
    theta3 = ab(  0,  3),      # crcl power scalar
    sat    = ab(0.5,  1)       # sat power scalar
  ),
  
  # ---------------- Covariates (names must match your dataset column headers) -----
  cov = list(
    interval_volume = interp("none"),  # L
    age             = interp(),        # years
    male            = interp("none"),  # 0/1
    ht              = interp(),        # cm or m
    wt              = interp(),        # kg
    scr             = interp(),        # mg/dL
    crcl            = interp(),        # mL/min
    bsa             = interp(),        # m^2
    crcl_bsa        = interp(),        # mL/min/1.73 m^2
    ecmo            = interp("none"),  # 0/1
    hd              = interp("none"),  # 0/1
    crrt            = interp("none"),  # 0/1
    cvvh            = interp("none"),  # 0/1
    cvvhd           = interp("none"),  # 0/1
    cvvhdf          = interp("none"),  # 0/1
    flow            = interp("none"),  # mL/h
    bfr             = interp("none"),  # mL/min (blood flow rate)
    bag_reset       = interp("none")   # 0/1
  ),
  sec = function(){
    CL2 = 7.2
    theta2 = 0.75
    # Weight-scaled central volume
    V_WT = V1 * (wt/70)
    
    # Systemic (renal ± HD)
    CL_R   = CL1 * (crcl/120)**theta3  * (wt/70)**theta2            # L/h
    CL_HD  = CL2 * hd                                               # L/h
    CL_sys = CL_R + CL_HD                                           # L/h
    
    # CRRT clearance term from effluent flow (mL/h -> L/h)
    Qe      = (flow/1000) * crrt      # L/h; 0 when CRRT=0
    Qe_eff  = (Qe/2.6)**sat
    CL_CRRT = S_eff * Qe_eff          # L/h
    
    # Total clearance and elimination from central
    CLT = CL_sys + CL_CRRT            # L/h
    Ke  = CLT / V_WT                  # 1/h
    
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    
    # mg/L pre-filter
    Cpre  = X[1]/V_WT
    # mg/L post-filter
    Cpost = Cpre * S_post * crrt
    
    Aeff = 0
    Ceff = 0
    
    # Guard against divide-by-zero
    if (X[4] > 1e-6) {
      Ceff = X[3] / X[4]
    }
    
    if (X[3] > 1e-6) {
      Aeff = X[3]
    }
    
    Celf = X[5]/V2
  },
  # ---------------- Differential equations ----------------
  eqn = function(){
    # Three compartment primary model structure
    dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2] - K15*X[1] + K51*X[5]
    dX[2] =                  K12*X[1] - K21*X[2]
    # Bag accumulation since last dump
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    dX[3] = Qe_eff * S_eff * (X[1]/V_WT) * Iacc  - kdump * X[3] * (1 - Iacc)
    dX[4] = Qe_eff * Iacc                        - kdump * X[4] * (1 - Iacc)
    dX[5] =                  K15*X[1] - K51*X[5]
  },
  
  # ---------------- Outputs ----------------
  out = function(){
    Y[1]  = Cpre       # DVID=1: pre/regular plasma in mg/L
    Y[2]  = Cpost      # DVID=2: post-filter plasma (0 if CRRT off)
    Y[3]  = Ceff       # DVID=3: effluent concentration (0 if not collecting)
    Y[4]  = Aeff       # DVID=4: effluent cumulative amount (mg)
    Y[5]  = Celf       # DVID=5: ELF concentrations from urea corrected BAL method
  },
  
  # ---------------- Error models ----------------
  # Proportional with a floor per-assay.
  err = list(
    proportional(2, c(1, 0.15, 0, 0)),   # Y1 (Cpre)
    proportional(2, c(1, 0.15, 0, 0)),   # Y2 (Cpost)
    proportional(2, c(1, 0.15, 0, 0)),   # Y3 (Ceff_out)
    proportional(2, c(1, 0.15, 0, 0)),   # Y4 (Aeff)
    proportional(2, c(1, 0.15, 0, 0))    # Y5 (Celf)
  )
)

run12 <- mod_cef8$fit(data = dat5,
                      cycles = 1000,
                      path = PATHS$runs,
                      run = 12,
                      points = 250,
                      overwrite=TRUE)

run12 <- PM_load(12, path = PATHS$runs)

PM_compare(run10,run12)

# cycles 0 hold out
run13 <- mod_cef8$fit(data = dat6,
                      cycles = 0,
                      prior = 12,
                      path = PATHS$runs,
                      run = 13,
                      overwrite=TRUE)

run13 <- PM_load(13, path = PATHS$runs)


# =============================================================================
# 6. FDA-FEP-ELF :: ELF-extended candidate models ####
# Insert after dat5/dat6 are created.
#
# Purpose:
# - Carry the original plasma/CRRT structural candidates forward into the
#   combined plasma + ELF training/test framework.
# - Fit each ELF-extended candidate on dat5.
# - Evaluate each candidate on dat6 using cycles = 0 and the corresponding prior.
#
# Notes:
# - run10/run11 already represent the ELF-extended version of mod_cef5.
# - run12/run13 already represent the ELF-extended version of mod_cef6.
# - This block adds ELF-extended versions of mod_cef1 through mod_cef4 using
#   the same explicit style as the current Analysis.R script.
# =============================================================================
#
# ELF-extension strategy:
# - mod_cef1 through mod_cef4 are explicitly rebuilt as mod_cef1_elf through mod_cef4_elf.
# - mod_cef7 is the existing ELF-extended form of mod_cef5.
# - mod_cef8 is the existing ELF-extended form of mod_cef6.
# - run10/run11 and run12/run13 are therefore retained rather than duplicated.
# =============================================================================

# =============================================================================
# ELF extension of mod_cef1 ####
# Base model: CrCL on CL, WT on Vd
# Training: run14
# Validation: run15
# =============================================================================

mod_cef1_elf <- PM_model$new(
  # ---------------- Priors / bounds ----------------
  pri = list(
    V1     = ab(0.1,  30),   # L (typical central, prior wide; scaled by weight in eqn)
    V2     = ab(0.1, 100),     # L (ELF volume of distribution)
    K12    = ab(  0,  10),   # 1/h
    K21    = ab(  0,  10),   # 1/h
    K15    = ab(  0,  10),     # 1/h central to ELF
    K51    = ab(  0,  10),     # 1/h ELF  to central
    CL1    = ab(0.1,  12),   # L/h (renal scale at CRCL=120)
    S_eff  = ab(  0,  2),    # dimensionless sieving scaler for effluent clearance
    S_post = ab(  0,  2)     # dimensionless scaler for post-filter clearance
  ),
  
  # ---------------- Covariates (names must match your dataset column headers) -----
  cov = list(
    interval_volume = interp("none"),  # L
    age             = interp(),        # years
    male            = interp("none"),  # 0/1
    ht              = interp(),        # cm or m
    wt              = interp(),        # kg
    scr             = interp(),        # mg/dL
    crcl            = interp(),        # mL/min
    bsa             = interp(),        # m^2
    crcl_bsa        = interp(),        # mL/min/1.73 m^2
    ecmo            = interp("none"),  # 0/1
    hd              = interp("none"),  # 0/1
    crrt            = interp("none"),  # 0/1
    cvvh            = interp("none"),  # 0/1
    cvvhd           = interp("none"),  # 0/1
    cvvhdf          = interp("none"),  # 0/1
    flow            = interp("none"),  # mL/h
    bfr             = interp("none"),  # mL/min (blood flow rate)
    bag_reset       = interp("none")   # 0/1
  ),
  sec = function(){
    CL2 = 7.2                          # Fixed HD clearance
    
    # Weight-scaled central volume
    V_WT = V1 * (wt/70)                # WT scaled Vd
    
    # Systemic (renal ± HD)
    CL_R   = CL1 * (crcl/120)                         # L/h
    CL_HD  = CL2 * hd                                 # L/h
    CL_sys = CL_R + CL_HD                             # L/h
    
    # CRRT clearance term from effluent flow (mL/h -> L/h)
    Qe      = (flow/1000) * crrt      # L/h; 0 when CRRT=0
    CL_CRRT = S_eff * Qe              # L/h
    
    # Total clearance and elimination from central
    CLT = CL_sys + CL_CRRT            # L/h
    Ke  = CLT / V_WT                  # 1/h
    
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    
    # mg/L pre-filter
    Cpre  = X[1]/V_WT
    # mg/L post-filter
    Cpost = Cpre * S_post * crrt
    
    Aeff = 0
    Ceff = 0
    
    # Guard against divide-by-zero
    if (X[4] > 1e-6) {
      Ceff = X[3] / X[4]
    }
    
    if (X[3] > 1e-6) {
      Aeff = X[3]
    }
    
    Celf = X[5]/V2
  },
  # ---------------- Differential equations ----------------
  eqn = function(){
    # Three compartment primary model structure
    dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2]
    dX[2] =                  K12*X[1] - K21*X[2]
    # Bag accumulation since last dump
    kdump = 1e4                       # very fast “drain” when bag_reset==0
    Iacc  = bag_reset * crrt          # 1 only when circuit ON and collecting
    dX[3] = Qe * S_eff * (X[1]/V_WT) * Iacc  - kdump * X[3] * (1 - Iacc)
    dX[4] = Qe * Iacc                        - kdump * X[4] * (1 - Iacc)
    dX[5] =                  K15*X[1] - K51*X[5]
  },
  
  # ---------------- Outputs ----------------
  out = function(){
    Y[1]  = Cpre       # DVID=1: pre/regular plasma in mg/L
    Y[2]  = Cpost      # DVID=3: post-filter plasma (0 if CRRT off)
    Y[3]  = Ceff       # DVID=4: effluent concentration (0 if not collecting)
    Y[4]  = Aeff       # DVID=5: effluent cumulative amount (mg)
    Y[5]  = Celf       # DVID=5: ELF concentrations from urea corrected BAL method
  },
  
  # ---------------- Error models ----------------
  # Proportional with a floor per-assay.
  err = list(
    proportional(2, c(1, 0.15, 0, 0)),   # Y1 (Cpre)
    proportional(2, c(1, 0.15, 0, 0)),   # Y2 (Cpost)
    proportional(2, c(1, 0.15, 0, 0)),   # Y3 (Ceff_out)
    proportional(2, c(1, 0.15, 0, 0)),    # Y4 (Aeff)
    proportional(2, c(1, 0.15, 0, 0))    # Y5 (Celf)
  )
)

run14 <- mod_cef1_elf$fit(data = dat5,
                          cycles = 1000,
                          path = PATHS$runs,
                          run = 14,
                          points = 300,
                          overwrite = TRUE)

run14 <- PM_load(14, path = PATHS$runs)

run15 <- mod_cef1_elf$fit(data = dat6,
                          cycles = 0,
                          path = PATHS$runs,
                          run = 15,
                          prior = 14,
                          overwrite = TRUE)

run15 <- PM_load(15, path = PATHS$runs)


# =============================================================================
# ELF extension of mod_cef2
# Adds BSA power scalar on renal CL
# Training: run16
# Validation: run17
# =============================================================================

mod_cef2_elf <- PM_model$new(
  # ---------------- Priors / bounds ----------------
  pri = list(
    V1     = ab(0.1,  30),
    V2     = ab(0.1, 100),
    K12    = ab(  0,  10),
    K21    = ab(  0,  10),
    K15    = ab(  0,  10),
    K51    = ab(  0,  10),
    CL1    = ab(0.1,  12),
    S_eff  = ab(  0,   2),
    S_post = ab(  0,   2),
    theta1 = ab(  0,   2)      # BSA clearance scaler
  ),
  
  # ---------------- Covariates ----------------
  cov = list(
    interval_volume = interp("none"),
    age             = interp(),
    male            = interp("none"),
    ht              = interp(),
    wt              = interp(),
    scr             = interp(),
    crcl            = interp(),
    bsa             = interp(),
    crcl_bsa        = interp(),
    ecmo            = interp("none"),
    hd              = interp("none"),
    crrt            = interp("none"),
    cvvh            = interp("none"),
    cvvhd           = interp("none"),
    cvvhdf          = interp("none"),
    flow            = interp("none"),
    bfr             = interp("none"),
    bag_reset       = interp("none")
  ),
  sec = function(){
    CL2 = 7.2
    
    # Weight-scaled central volume
    V_WT = V1 * (wt/70)
    
    # Systemic (renal +/- HD)
    CL_R   = CL1 * (crcl/120) * (bsa/1.73)**theta1
    CL_HD  = CL2 * hd
    CL_sys = CL_R + CL_HD
    
    # CRRT clearance term from effluent flow (mL/h -> L/h)
    Qe      = (flow/1000) * crrt
    CL_CRRT = S_eff * Qe
    
    # Total clearance and elimination from central
    CLT = CL_sys + CL_CRRT
    Ke  = CLT / V_WT
    
    kdump = 1e4
    Iacc  = bag_reset * crrt
    
    # Outputs / derived concentrations
    Cpre  = X[1]/V_WT
    Cpost = Cpre * S_post * crrt
    
    Aeff = 0
    Ceff = 0
    
    if (X[4] > 1e-6) {
      Ceff = X[3] / X[4]
    }
    
    if (X[3] > 1e-6) {
      Aeff = X[3]
    }
    
    Celf = X[5]/V2
  },
  eqn = function(){
    # Central + peripheral + ELF disposition model
    dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2] - K15*X[1] + K51*X[5]
    dX[2] =                  K12*X[1] - K21*X[2]
    
    # Bag accumulation since last dump
    kdump = 1e4
    Iacc  = bag_reset * crrt
    dX[3] = Qe * S_eff * (X[1]/V_WT) * Iacc  - kdump * X[3] * (1 - Iacc)
    dX[4] = Qe * Iacc                        - kdump * X[4] * (1 - Iacc)
    
    # ELF exchange
    dX[5] =                  K15*X[1] - K51*X[5]
  },
  out = function(){
    Y[1] = Cpre
    Y[2] = Cpost
    Y[3] = Ceff
    Y[4] = Aeff
    Y[5] = Celf
  },
  err = list(
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0))
  )
)

run16 <- mod_cef2_elf$fit(data = dat5,
                          cycles = 1000,
                          path = PATHS$runs,
                          run = 16,
                          points = 300,
                          overwrite = TRUE)

run16 <- PM_load(16, path = PATHS$runs)

run17 <- mod_cef2_elf$fit(data = dat6,
                          cycles = 0,
                          path = PATHS$runs,
                          run = 17,
                          prior = 16,
                          overwrite = TRUE)

run17 <- PM_load(17, path = PATHS$runs)


# =============================================================================
# ELF extension of mod_cef3
# Estimates CrCL power scalar on CL
# Training: run18
# Validation: run19
# =============================================================================

mod_cef3_elf <- PM_model$new(
  # ---------------- Priors / bounds ----------------
  pri = list(
    V1     = ab(0.1,  30),
    V2     = ab(0.1, 100),
    K12    = ab(  0,  10),
    K21    = ab(  0,  10),
    K15    = ab(  0,  10),
    K51    = ab(  0,  10),
    CL1    = ab(0.1,  12),
    S_eff  = ab(  0,   2),
    S_post = ab(  0,   2),
    theta3 = ab(  0,   3)      # CrCL power scaler
  ),
  
  # ---------------- Covariates ----------------
  cov = list(
    interval_volume = interp("none"),
    age             = interp(),
    male            = interp("none"),
    ht              = interp(),
    wt              = interp(),
    scr             = interp(),
    crcl            = interp(),
    bsa             = interp(),
    crcl_bsa        = interp(),
    ecmo            = interp("none"),
    hd              = interp("none"),
    crrt            = interp("none"),
    cvvh            = interp("none"),
    cvvhd           = interp("none"),
    cvvhdf          = interp("none"),
    flow            = interp("none"),
    bfr             = interp("none"),
    bag_reset       = interp("none")
  ),
  sec = function(){
    CL2 = 7.2
    
    # Weight-scaled central volume
    V_WT = V1 * (wt/70)
    
    # Systemic (renal +/- HD)
    CL_R   = CL1 * (crcl/120)**theta3
    CL_HD  = CL2 * hd
    CL_sys = CL_R + CL_HD
    
    # CRRT clearance term from effluent flow (mL/h -> L/h)
    Qe      = (flow/1000) * crrt
    CL_CRRT = S_eff * Qe
    
    # Total clearance and elimination from central
    CLT = CL_sys + CL_CRRT
    Ke  = CLT / V_WT
    
    kdump = 1e4
    Iacc  = bag_reset * crrt
    
    # Outputs / derived concentrations
    Cpre  = X[1]/V_WT
    Cpost = Cpre * S_post * crrt
    
    Aeff = 0
    Ceff = 0
    
    if (X[4] > 1e-6) {
      Ceff = X[3] / X[4]
    }
    
    if (X[3] > 1e-6) {
      Aeff = X[3]
    }
    
    Celf = X[5]/V2
  },
  eqn = function(){
    # Central + peripheral + ELF disposition model
    dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2] - K15*X[1] + K51*X[5]
    dX[2] =                  K12*X[1] - K21*X[2]
    
    # Bag accumulation since last dump
    kdump = 1e4
    Iacc  = bag_reset * crrt
    dX[3] = Qe * S_eff * (X[1]/V_WT) * Iacc  - kdump * X[3] * (1 - Iacc)
    dX[4] = Qe * Iacc                        - kdump * X[4] * (1 - Iacc)
    
    # ELF exchange
    dX[5] =                  K15*X[1] - K51*X[5]
  },
  out = function(){
    Y[1] = Cpre
    Y[2] = Cpost
    Y[3] = Ceff
    Y[4] = Aeff
    Y[5] = Celf
  },
  err = list(
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0))
  )
)

run18 <- mod_cef3_elf$fit(data = dat5,
                          cycles = 1000,
                          path = PATHS$runs,
                          run = 18,
                          points = 300,
                          overwrite = TRUE)

run18 <- PM_load(18, path = PATHS$runs)

run19 <- mod_cef3_elf$fit(data = dat6,
                          cycles = 0,
                          path = PATHS$runs,
                          run = 19,
                          prior = 18,
                          overwrite = TRUE)

run19 <- PM_load(19, path = PATHS$runs)


# =============================================================================
# ELF extension of mod_cef4
# Fixed allometric WT coefficient on CL: (WT/70)^0.75
# Training: run20
# Validation: run21
# =============================================================================

mod_cef4_elf <- PM_model$new(
  # ---------------- Priors / bounds ----------------
  pri = list(
    V1     = ab(0.1,  30),
    V2     = ab(0.1, 100),
    K12    = ab(  0,  10),
    K21    = ab(  0,  10),
    K15    = ab(  0,  10),
    K51    = ab(  0,  10),
    CL1    = ab(0.1,  12),
    S_eff  = ab(  0,   2),
    S_post = ab(  0,   2)
  ),
  
  # ---------------- Covariates ----------------
  cov = list(
    interval_volume = interp("none"),
    age             = interp(),
    male            = interp("none"),
    ht              = interp(),
    wt              = interp(),
    scr             = interp(),
    crcl            = interp(),
    bsa             = interp(),
    crcl_bsa        = interp(),
    ecmo            = interp("none"),
    hd              = interp("none"),
    crrt            = interp("none"),
    cvvh            = interp("none"),
    cvvhd           = interp("none"),
    cvvhdf          = interp("none"),
    flow            = interp("none"),
    bfr             = interp("none"),
    bag_reset       = interp("none")
  ),
  sec = function(){
    CL2 = 7.2
    theta2 = 0.75
    
    # Weight-scaled central volume
    V_WT = V1 * (wt/70)
    
    # Systemic (renal +/- HD)
    CL_R   = CL1 * (crcl/120) * (wt/70)**theta2
    CL_HD  = CL2 * hd
    CL_sys = CL_R + CL_HD
    
    # CRRT clearance term from effluent flow (mL/h -> L/h)
    Qe      = (flow/1000) * crrt
    CL_CRRT = S_eff * Qe
    
    # Total clearance and elimination from central
    CLT = CL_sys + CL_CRRT
    Ke  = CLT / V_WT
    
    kdump = 1e4
    Iacc  = bag_reset * crrt
    
    # Outputs / derived concentrations
    Cpre  = X[1]/V_WT
    Cpost = Cpre * S_post * crrt
    
    Aeff = 0
    Ceff = 0
    
    if (X[4] > 1e-6) {
      Ceff = X[3] / X[4]
    }
    
    if (X[3] > 1e-6) {
      Aeff = X[3]
    }
    
    Celf = X[5]/V2
  },
  eqn = function(){
    # Central + peripheral + ELF disposition model
    dX[1] = b[1] + R[1] - Ke*X[1] - K12*X[1] + K21*X[2] - K15*X[1] + K51*X[5]
    dX[2] =                  K12*X[1] - K21*X[2]
    
    # Bag accumulation since last dump
    kdump = 1e4
    Iacc  = bag_reset * crrt
    dX[3] = Qe * S_eff * (X[1]/V_WT) * Iacc  - kdump * X[3] * (1 - Iacc)
    dX[4] = Qe * Iacc                        - kdump * X[4] * (1 - Iacc)
    
    # ELF exchange
    dX[5] =                  K15*X[1] - K51*X[5]
  },
  out = function(){
    Y[1] = Cpre
    Y[2] = Cpost
    Y[3] = Ceff
    Y[4] = Aeff
    Y[5] = Celf
  },
  err = list(
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0)),
    proportional(2, c(1, 0.15, 0, 0))
  )
)

run20 <- mod_cef4_elf$fit(data = dat5,
                          cycles = 1000,
                          path = PATHS$runs,
                          run = 20,
                          points = 300,
                          overwrite = TRUE)

run20 <- PM_load(20, path = PATHS$runs)

run21 <- mod_cef4_elf$fit(data = dat6,
                          cycles = 0,
                          path = PATHS$runs,
                          run = 21,
                          prior = 20,
                          overwrite = TRUE)

run21 <- PM_load(21, path = PATHS$runs)


# =============================================================================
# Existing ELF candidates already in the current script
# =============================================================================

# mod_cef7 / run10 / run11 is the ELF-extended version of mod_cef5:
#   fixed WT allometry on CL + estimated CrCL power scalar
#   no CRRT saturation term

# mod_cef8 / run12 / run13 is the ELF-extended version of mod_cef6:
#   fixed WT allometry on CL + estimated CrCL power scalar
#   adds effluent-flow saturation term

# IMPORTANT PATCH for current script:
# After fitting run13, load run13, not run12.
# Replace this:
#   run12 <- PM_load(12, path = PATHS$runs)
# with this:
#   run13 <- PM_load(13, path = PATHS$runs)

# Optional aliases for readability in comparison tables.
mod_cef5_elf <- mod_cef7
mod_cef6_elf <- mod_cef8


# =============================================================================
# Training model comparison across ELF-extended structural candidates
# =============================================================================

library(dplyr)
library(tibble)
library(knitr)
library(kableExtra)

# Do NOT compare all 6 models in one PM_compare() call.
# PM_compare() may fail internally when assigning the "best" column across 6 rows.
# Use the same split-comparison approach used earlier for runs 4-9.

cmp_elf_1 <- PM_compare(run14, run16, run18, run20)
cmp_elf_2 <- PM_compare(run10, run12)

compare_elf_tbl <- tibble(
  model = c(
    "mod_cef1_elf",
    "mod_cef2_elf",
    "mod_cef3_elf",
    "mod_cef4_elf",
    "mod_cef5_elf / mod_cef7",
    "mod_cef6_elf / mod_cef8"
  ),
  train_run = c(14, 16, 18, 20, 10, 12),
  validation_run = c(15, 17, 19, 21, 11, 13),
  nvar = c(cmp_elf_1$nvar, cmp_elf_2$nvar),
  ll2  = c(cmp_elf_1$`-2*LL`, cmp_elf_2$`-2*LL`),
  aic  = c(cmp_elf_1$aic, cmp_elf_2$aic),
  notes = c(
    "base: CL1*(crcl/120) + S_eff*Qe; V1*(wt/70); ELF K15/K51",
    "add BSA power scalar to renal CL: CL1*(crcl/120)*(bsa/1.73)^theta1; ELF K15/K51",
    "estimate CRCL power scalar: CL1*(crcl/120)^theta3; ELF K15/K51",
    "fix WT allometry (0.75) on renal CL: CL1*(crcl/120)*(wt/70)^0.75; ELF K15/K51",
    "combine: CL1*(crcl/120)^theta3*(wt/70)^0.75; ELF K15/K51",
    "add effluent flow saturation: CL_CRRT via (Qe/2.6)^sat; ELF K15/K51"
  )
) %>%
  mutate(
    delta_ofv_best = ll2 - min(ll2, na.rm = TRUE),
    delta_aic      = aic - min(aic, na.rm = TRUE)
  ) %>%
  arrange(aic)

compare_elf_tbl

compare_elf_tbl %>%
  kable(
    digits  = 1,
    caption = "Model comparison across ELF-extended structural candidates in the training cohort",
    align   = "c",
    format  = "html"
  ) %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  save_kable(file.path(PATHS$tab, "Table_S1_ELF_extended_model_comparison.html"))

write.csv(
  compare_elf_tbl,
  file.path(PATHS$tab, "Table_S1_ELF_extended_model_comparison.csv"),
  row.names = FALSE
)


# =============================================================================
# Validation performance summary for dat6 candidates ####
# =============================================================================

summarize_op <- function(run, run_label = NA_character_) {
  
  op_df <- run$op$data %>%
    as_tibble() %>%
    filter(
      pred.type == "post",
      icen == "median",
      outeq %in% c(1, 2, 5),
      !is.na(obs),
      !is.na(pred),
      obs > 0
    ) %>%
    mutate(
      pe      = pred - obs,
      ape     = abs(pred - obs),
      rpe     = 100 * (pred - obs) / obs,
      ape_pct = 100 * abs(pred - obs) / obs
    )
  
  op_df %>%
    group_by(outeq) %>%
    summarise(
      n = n(),
      mpe = median(pe, na.rm = TRUE),
      mape = median(ape, na.rm = TRUE),
      rMPE = median(rpe, na.rm = TRUE),
      rMAPE = median(ape_pct, na.rm = TRUE),
      r2 = {
        d <- cur_data()
        if (
          nrow(d) >= 2 &&
          length(unique(d$obs)) > 1 &&
          length(unique(d$pred)) > 1
        ) {
          suppressWarnings(summary(lm(obs ~ pred, data = d))$r.squared)
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    ) %>%
    mutate(
      validation_run = run_label,
      matrix = case_when(
        outeq == 1 ~ "Pre-filter plasma",
        outeq == 2 ~ "Post-filter plasma",
        outeq == 5 ~ "ELF",
        TRUE ~ paste0("OUTEQ ", outeq)
      )
    ) %>%
    select(validation_run, matrix, outeq, n, mpe, mape, rMPE, rMAPE, r2)
}

validation_elf_tbl <- bind_rows(
  summarize_op(run15, "run15: mod_cef1_elf"),
  summarize_op(run17, "run17: mod_cef2_elf"),
  summarize_op(run19, "run19: mod_cef3_elf"),
  summarize_op(run21, "run21: mod_cef4_elf"),
  summarize_op(run11, "run11: mod_cef5_elf / mod_cef7"),
  summarize_op(run13, "run13: mod_cef6_elf / mod_cef8")
)
# Model-selection note:
# - Training AIC favored mod_cef5_elf / mod_cef7 in the current results.
# - Held-out ELF validation favored mod_cef1_elf in the current results.
# - Final model selection should consider training fit, plasma validation,
#   ELF validation, physiologic plausibility, and the small ELF validation sample size.

validation_elf_tbl

# View(validation_elf_tbl)  # interactive-only; keep commented for reproducible script

validation_elf_tbl %>%
  kable(
    digits  = 2,
    caption = "Validation performance across ELF-extended structural candidates",
    align   = "c",
    format  = "html"
  ) %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  save_kable(file.path(PATHS$tab, "Table_S3_ELF_extended_validation_performance.html"))

write.csv(
  validation_elf_tbl,
  file.path(PATHS$tab, "Table_S3_ELF_extended_validation_performance.csv"),
  row.names = FALSE
)


library(DT)

DT::datatable(compare_elf_tbl)
DT::datatable(validation_elf_tbl)
# =============================================================================
# 8. Sensitivity diagnostic: paired ELF validation comparison
# =============================================================================

elf_val_compare <- bind_rows(
  run15$op$data %>%
    as_tibble() %>%
    filter(pred.type == "post", icen == "median", outeq == 5, !is.na(obs), !is.na(pred), obs > 0) %>%
    mutate(model = "run15: mod_cef1_elf"),
  
  run11$op$data %>%
    as_tibble() %>%
    filter(pred.type == "post", icen == "median", outeq == 5, !is.na(obs), !is.na(pred), obs > 0) %>%
    mutate(model = "run11: mod_cef5_elf / mod_cef7")
) %>%
  mutate(
    pe = pred - obs,
    rpe = 100 * (pred - obs) / obs,
    ape_pct = 100 * abs(pred - obs) / obs
  ) %>%
  select(model, id, time, obs, pred, pe, rpe, ape_pct) %>%
  arrange(id, time, model)

elf_val_compare


# 0) Fresh start (recommended if you've been experimenting)
# Restart R session, then run everything below.

# 1) Load reticulate
library(reticulate)

# 2) Create a clean conda env for kaleido
#(pin Python; 3.9 or 3.10 both fine with kaleido 0.1.0)
conda_create("r-reticulate", packages = "python=3.10")

#3) Install kaleido *0.1.0* from conda-forge into that env
conda_install(
  envname  = "r-reticulate",
  packages = c("python-kaleido=0.1.0", "plotly=5.19.*"),
  channel  = "conda-forge"
)

# 4) Point reticulate at this env for the current R session
use_condaenv("r-reticulate", required = TRUE)

# 5) (Important) Nudge reticulate to fully init Python
reticulate::py_run_string('import sys')

# 6) Sanity checks — confirm Python and kaleido version
reticulate::py_config()
reticulate::py_run_string('import kaleido, sys; print("kaleido:", kaleido.__version__);import plotly, sys; print("plotly:", plotly.__version__); print("python:", sys.executable)')

# plot the OP results for the training model run10
library(plotly)

f <- list(size = 20)

a <- list(
  text = "A. Population",
  font=f,
  xref = "paper",
  yref = "paper",
  yanchor = "bottom",
  xanchor = "center",
  align = "right",
  x = 0.075,
  y = 1,
  showarrow = FALSE
)

b <- list(
  text = "B. Posterior",
  font=f,
  xref = "paper",
  yref = "paper",
  yanchor = "bottom",
  xanchor = "center",
  align = "right",
  x = 0.05,
  y = 1,
  showarrow = FALSE
)

op1<-run10$op$plot(pred.type="pop",outeq=1,
                   line=list(lm=list(T,color="dodgerblue")),
                   marker=list(color="goldenrod"),
                   stats = list(x = 0.75, y = 0.25, font = list(size = 10)))

op2<-run10$op$plot(pred.type="post",outeq=1,
                   line=list(lm=list(T,color="dodgerblue")),
                   marker=list(color="goldenrod"),
                   stats = list(x = 0.95, y = 0.25, font = list(size = 10)))

op1 <- op1 %>% layout(annotations = a)
op2 <- op2 %>% layout(annotations = b)

sub_plot(op1, op2, nrows = 1, 
         margin=0.05, titleX=T, titleY=T,shareX=T,shareY=T) %>%
  export_plotly("FEP PRE-FILTER FDA RENAL Plot_test.svg", 
                width = 5.5 * 300, 
                height = 3.5 * 300, 
                scale=2)


op1<-run10$op$plot(pred.type="pop",outeq=2,
                   line=list(lm=list(T,color="dodgerblue")),
                   marker=list(color="goldenrod"),
                   stats = list(x = 0.75, y = 0.25, font = list(size = 10)))

op2<-run10$op$plot(pred.type="post",outeq=2,
                   line=list(lm=list(T,color="dodgerblue")),
                   marker=list(color="goldenrod"),
                   stats = list(x = 0.95, y = 0.25, font = list(size = 10)))

op1 <- op1 %>% layout(annotations = a)
op2 <- op2 %>% layout(annotations = b)

sub_plot(op1, op2, nrows = 1, margin=0.05, titleX=T, titleY=T,shareX=T,shareY=T) %>%
  export_plotly("FEP POST FILTER FDA RENAL Plot test.svg", width = 3.5 * 300, height = 1.5 * 300)


op1<-run10$op$plot(pred.type="pop",outeq=3,
                   line=list(lm=list(T,color="dodgerblue")),
                   marker=list(color="goldenrod"),
                   stats = list(x = 0.75, y = 0.25, font = list(size = 10)))

op2<-run10$op$plot(pred.type="post",outeq=3,
                   line=list(lm=list(T,color="dodgerblue")),
                   marker=list(color="goldenrod"),
                   stats = list(x = 0.95, y = 0.25, font = list(size = 10)))

op1 <- op1 %>% layout(annotations = a)
op2 <- op2 %>% layout(annotations = b)

sub_plot(op1, op2, nrows = 1, margin=0.05, titleX=T, titleY=T,shareX=T,shareY=T) %>%
  export_plotly("FEP EFFLUENT CONC FDA RENAL Plot test.svg", width = 3.5 * 300, height = 1.5 * 300)

op1<-run10$op$plot(pred.type="pop",outeq=4,
                   line=list(lm=list(T,color="dodgerblue")),
                   marker=list(color="goldenrod"),
                   stats = list(x = 0.75, y = 0.25, font = list(size = 10)))

op2<-run10$op$plot(pred.type="post",outeq=4,
                   line=list(lm=list(T,color="dodgerblue")),
                   marker=list(color="goldenrod"),
                   stats = list(x = 0.95, y = 0.25, font = list(size = 10)))

op1 <- op1 %>% layout(annotations = a)
op2 <- op2 %>% layout(annotations = b)

sub_plot(op1, op2, nrows = 1, margin=0.05, titleX=T, titleY=T,shareX=T,shareY=T) %>%
  export_plotly("FEP EFFLUENT AMT FDA RENAL Plot test.svg", width = 3.5 * 300, height = 1.5 * 300)

op1<-run10$op$plot(pred.type="pop",outeq=5,
                   line=list(lm=list(T,color="dodgerblue")),
                   marker=list(color="goldenrod"),
                   stats = list(x = 0.75, y = 0.25, font = list(size = 10)))

op2<-run10$op$plot(pred.type="post",outeq=5,
                   line=list(lm=list(T,color="dodgerblue")),
                   marker=list(color="goldenrod"),
                   stats = list(x = 0.95, y = 0.25, font = list(size = 10)))

op1 <- op1 %>% layout(annotations = a)
op2 <- op2 %>% layout(annotations = b)

sub_plot(op1, op2, nrows = 1, margin=0.05, titleX=T, titleY=T,shareX=T,shareY=T) %>%
  export_plotly("FEP ELF CONC FDA RENAL Plot test.svg", width = 3.5 * 300, height = 1.5 * 300)

# # Make the post-hoc run 11 using run 10 prior plots plotly sub plot is sucking
# f <- list(size = 18)
# 
# a <- list(
#   text = "A. Pre-filter Plasma",
#   font=f,
#   xref = "paper",
#   yref = "paper",
#   yanchor = "bottom",
#   xanchor = "center",
#   align = "right",
#   x = 0.1,
#   y = 1,
#   showarrow = FALSE
# )
# 
# b <- list(
#   text = "B. Post-filter Plasma",
#   font=f,
#   xref = "paper",
#   yref = "paper",
#   yanchor = "bottom",
#   xanchor = "center",
#   align = "right",
#   x = 0.1,
#   y = 1,
#   showarrow = FALSE
# )
# 
# c <- list(
#   text = "C. ELF",
#   font=f,
#   xref = "paper",
#   yref = "paper",
#   yanchor = "bottom",
#   xanchor = "center",
#   align = "right",
#   x = 0.1,
#   y = 1,
#   showarrow = FALSE
# )
# 
# op1<-run11$op$plot(pred.type="post",outeq=1,
#                    line=list(lm=list(T,color="dodgerblue")),
#                    marker=list(color="goldenrod"),
#                    stats = list(x = 0.55, y = 0.18, font = list(size = 10)))
# 
# op2<-run11$op$plot(pred.type="post",outeq=2,
#                    line=list(lm=list(T,color="dodgerblue")),
#                    marker=list(color="goldenrod"),
#                    stats = list(x = 0.85, y = 0.18, font = list(size = 10)))
# 
# op3<-run11$op$plot(pred.type="post",outeq=5,
#                    line=list(lm=list(T,color="dodgerblue")),
#                    marker=list(color="goldenrod"),
#                    stats = list(x = 1, y = 0.18, font = list(size = 10)))
# 
# op1 <- op1 %>% layout(annotations = a)
# op2 <- op2 %>% layout(annotations = b)
# op3 <- op3 %>% layout(annotations = c)
# 
# 
# sub_plot(op1, op2, op3, 
#          nrows = 1,
#          margin=0.05, 
#          titleX=T, 
#          titleY=T,
#          shareX=T,
#          shareY=T) %>%
#   export_plotly("FEP DEVELOPMENT FITS PRE POST ELF.svg", 
#                 width = 5 * 300, 
#                 height = 1.5 * 300,
#                 scale = 3)
# 
# # Make the post-hoc development plot using run10 fit
# # this method is trash for these plots because it suggest bias where it does not exist
# op1 <- run10$op$plot(
#   pred.type = "post", outeq = 1,
#   line   = list(lm = list(TRUE, color = "dodgerblue")),
#   marker = list(color = "goldenrod"),
#   stats  = list(x = 0.65, y = 0.25, font = list(size = 10))
# )
# 
# op2 <- run10$op$plot(
#   pred.type = "post", outeq = 2,
#   line   = list(lm = list(TRUE, color = "dodgerblue")),
#   marker = list(color = "goldenrod"),
#   stats  = list(x = 0.75, y = 0.25, font = list(size = 10))
# )
# 
# op3 <- run10$op$plot(
#   pred.type = "post", outeq = 5,
#   line   = list(lm = list(TRUE, color = "dodgerblue")),
#   marker = list(color = "goldenrod"),
#   stats  = list(x = 0.95, y = 0.25, font = list(size = 10))
# )
# 
# op1 <- op1 %>% layout(annotations = a)
# op2 <- op2 %>% layout(annotations = b)
# op3 <- op3 %>% layout(annotations = c)
# 
# comb <- subplot(
#   op1, op2, op3,
#   nrows  = 1,
#   margin = 0.05,
#   titleX = TRUE,
#   titleY = TRUE,
#   shareX = FALSE,
#   shareY = FALSE
# )
# 
# ## --- PATCH THE SCALEANCHORs SO EACH PANEL IS 1:1 ---
# 
# # Find all xaxes in the combined layout
# x_axes <- grep("^xaxis", names(comb$x$layout), value = TRUE)
# 
# for (x_name in x_axes) {
#   # corresponding yaxis name: xaxis -> yaxis, xaxis2 -> yaxis2, etc.
#   y_name <- sub("xaxis", "yaxis", x_name)
#   
#   # anchor that yaxis to *its own* x
#   # e.g. xaxis2 -> scaleanchor = "x2"
#   anchor_id <- sub("xaxis", "x", x_name)
#   
#   comb$x$layout[[y_name]]$scaleanchor <- anchor_id
#   comb$x$layout[[y_name]]$scaleratio  <- 1
# }
# 
# comb %>%
#   export_plotly(
#     "FEP DEVELOPMENT FITS PRE POST ELF.svg",
#     width  = 5.0 * 300,
#     height = 1.5 * 300,
#     scale = 3
#   )

# Make OPs for compartments of interest for development and validation cohorts ####
library(dplyr)
library(ggplot2)
library(broom)
library(patchwork)

# Grab raw op data (adjust this if your structure is slightly different)
op_raw <- run10$op$data

# Tidy it up: one row per observation
op_df <- op_raw %>%
  filter(outeq== 1 | outeq==2 | outeq==5,icen=="median", pred.type=="post") %>% 
  mutate(
    Pred = pred,   # change to your actual column names
    Obs  = obs
  )

max_val <- max(op_df$Pred, op_df$Obs, na.rm = TRUE)
# round up to nearest 10 for nice axes
lim_val <- ceiling(max_val / 10) * 10


make_op_panel <- function(dat, eq, panel_title, lim) {
  d <- dat %>%
    filter(outeq == eq) %>%
    mutate(
      Censor = if_else(cens == "bloq", "BLOQ", "Observed")
    )
  
  # Fit regression for stats
  fit        <- lm(Obs ~ Pred, data = d)
  fit_tidy   <- broom::tidy(fit)
  fit_glance <- broom::glance(fit)
  
  slope <- fit_tidy$estimate[fit_tidy$term == "Pred"]
  int   <- fit_tidy$estimate[fit_tidy$term == "(Intercept)"]
  r2    <- fit_glance$r.squared
  
  stat_text <- sprintf(
    "R² = %.3f\nIntercept = %.2f\nSlope = %.3f",
    r2, int, slope
  )
  
  ggplot(d, aes(x = Pred, y = Obs)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    geom_point(
      aes(fill = Censor, shape = Censor),
      color  = "black",
      size   = 2.8,
      stroke = 0.4,
      alpha  = 0.9
    ) +
    geom_smooth(
      method = "lm", se = TRUE,
      colour = "dodgerblue",
      fill   = "lightblue"
    ) +
    annotate(
      "text",
      x = Inf, y = -Inf,
      label = stat_text,
      hjust = 1.1, vjust = -0.1,
      size = 3.2
    ) +
    scale_fill_manual(
      values = c(
        "Observed" = "goldenrod",
        "BLOQ"     = "seagreen3"
      ),
      name = NULL
    ) +
    scale_shape_manual(
      values = c(
        "Observed" = 21,  # filled circle
        "BLOQ"     = 24   # filled triangle
      ),
      name = NULL
    ) +
    labs(
      title = panel_title,
      x = "Predicted",
      y = "Observed"
    ) +
    coord_equal(
      xlim   = c(0, lim),
      ylim   = c(0, lim),
      expand = FALSE
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0),
      panel.border    = element_rect(fill = NA, colour = "black"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background= element_rect(fill = "white", colour = NA),
      legend.position = "none"
    )
}


p_pre  <- make_op_panel(op_df, eq = 1, panel_title = "A. Pre-filter Plasma",  lim = lim_val)
p_post <- make_op_panel(op_df, eq = 2, panel_title = "B. Post-filter Plasma", lim = lim_val)
p_elf  <- make_op_panel(op_df, eq = 5, panel_title = "C. ELF",                lim = lim_val)

three_panel <- p_pre + p_post + p_elf +
  plot_layout(nrow = 1) +
  plot_annotation(
    title    = "Cefepime Model Fits In Development Cohort (n=40)",
    subtitle = "Observed vs predicted for pre-filter plasma, post-filter plasma, and ELF"
  )

three_panel

ggsave(
  "Manuscript/Figures/FEP_DEVELOPMENT_FITS_PRE_POST_ELF.svg",
  plot   = three_panel,
  width  = 11,   # tweak to taste
  height = 3.5,
  units  = "in"
)

ggsave(
  "Manuscript/Figures/FEP_DEVELOPMENT_FITS_PRE_POST_ELF.tiff",
  plot   = three_panel,
  width  = 11,
  height = 3.5,
  units  = "in",
  dpi    = 600,
  compression = "lzw"
)

# Now do validation cohort #
# Grab raw op data (adjust this if your structure is slightly different)
op_raw <- run11$op$data

# Tidy it up: one row per observation
op_df <- op_raw %>%
  filter(outeq== 1 | outeq==2 | outeq==5,icen=="median", pred.type=="post") %>% 
  mutate(
    Pred = pred,   # change to your actual column names
    Obs  = obs
  )

max_val <- max(op_df$Pred, op_df$Obs, na.rm = TRUE)
# round up to nearest 10 for nice axes
lim_val <- ceiling(max_val / 10) * 10


make_op_panel <- function(dat, eq, panel_title, lim) {
  d <- dat %>%
    filter(outeq == eq) %>%
    mutate(
      Censor = if_else(cens == "bloq", "BLOQ", "Observed")
    )
  
  # Fit regression for stats
  fit        <- lm(Obs ~ Pred, data = d)
  fit_tidy   <- broom::tidy(fit)
  fit_glance <- broom::glance(fit)
  
  slope <- fit_tidy$estimate[fit_tidy$term == "Pred"]
  int   <- fit_tidy$estimate[fit_tidy$term == "(Intercept)"]
  r2    <- fit_glance$r.squared
  
  stat_text <- sprintf(
    "R² = %.3f\nIntercept = %.2f\nSlope = %.3f",
    r2, int, slope
  )
  
  ggplot(d, aes(x = Pred, y = Obs)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    geom_point(
      aes(fill = Censor, shape = Censor),
      color  = "black",
      size   = 2.8,
      stroke = 0.4,
      alpha  = 0.9
    ) +
    geom_smooth(
      method = "lm", se = TRUE,
      colour = "dodgerblue",
      fill   = "lightblue"
    ) +
    annotate(
      "text",
      x = Inf, y = -Inf,
      label = stat_text,
      hjust = 1.1, vjust = -0.1,
      size = 3.2
    ) +
    scale_fill_manual(
      values = c(
        "Observed" = "goldenrod",
        "BLOQ"     = "seagreen3"
      ),
      name = NULL
    ) +
    scale_shape_manual(
      values = c(
        "Observed" = 21,  # filled circle
        "BLOQ"     = 24   # filled triangle
      ),
      name = NULL
    ) +
    labs(
      title = panel_title,
      x = "Predicted",
      y = "Observed"
    ) +
    coord_equal(
      xlim   = c(0, lim),
      ylim   = c(0, lim),
      expand = FALSE
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0),
      panel.border    = element_rect(fill = NA, colour = "black"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background= element_rect(fill = "white", colour = NA),
      legend.position = "none"
    )
}


p_pre  <- make_op_panel(op_df, eq = 1, panel_title = "A. Pre-filter Plasma",  lim = lim_val)
p_post <- make_op_panel(op_df, eq = 2, panel_title = "B. Post-filter Plasma", lim = lim_val)
p_elf  <- make_op_panel(op_df, eq = 5, panel_title = "C. ELF",                lim = lim_val)

three_panel <- p_pre + p_post + p_elf +
  plot_layout(nrow = 1) +
  plot_annotation(
    title    = "Cefepime Model Fits In Validation Cohort (n=11)",
    subtitle = "Observed vs predicted for pre-filter plasma, post-filter plasma, and ELF"
  )

three_panel

ggsave(
  "Manuscript/Figures/FEP_VALIDATION_FITS_PRE_POST_ELF.svg",
  plot   = three_panel,
  width  = 11,   # tweak to taste
  height = 3.5,
  units  = "in"
)

ggsave(
  "Manuscript/Figures/FEP_VALIDATION_FITS_PRE_POST_ELF.tiff",
  plot   = three_panel,
  width  = 11,
  height = 3.5,
  units  = "in",
  dpi    = 600,
  compression = "lzw"
)



# Table 1 generation code ####
library(dplyr)
library(purrr)
library(stringr)
library(gt)
library(scales)

# ---- combined data ----
dat <- run10$data$standard_data
val <- run11$data$standard_data
comb <- rbind(dat, val)

# ---- config ----
numeric_order <- c("age","ht","wt","scr","crcl","bsa","crcl_bsa")
binary_order  <- c("male","ecmo","hd","crrt","cvvh","cvvhd","cvvhdf")

numeric_vars <- intersect(numeric_order, names(comb))
binary_vars  <- intersect(binary_order,  names(comb))

labels <- c(
  age="Age", 
  ht="Height", 
  wt="Weight",
  scr="Serum creatinine", 
  crcl="Creatinine clearance",
  bsa="Body surface area", 
  crcl_bsa="CrCl (per 1.73 m²)",
  male="Male",
  ecmo="ECMO", 
  hd="Hemodialysis", 
  crrt="CRRT",
  cvvh="CVVH", 
  cvvhd="CVVHD", 
  cvvhdf="CVVHDF"
)
units <- c(
  age="years", 
  ht="cm", 
  wt="kg",
  scr="mg/dL", 
  crcl="mL/min",
  bsa="m²", 
  crcl_bsa="mL/min/1.73 m²"
)
lab <- function(v) ifelse(v %in% names(labels), labels[[v]], str_to_title(gsub("_"," ",v)))
uni <- function(v) ifelse(v %in% names(units),  units[[v]], "")

# ---- baseline subjects: time == 0 or earliest ----
baseline <- comb %>% filter(time == 0) %>% group_by(id) %>% slice_head(n = 1) %>% ungroup()
if (nrow(baseline) == 0) {
  baseline <- comb %>% group_by(id) %>% arrange(time, .by_group = TRUE) %>% slice_head(n = 1) %>% ungroup()
}
ids0   <- unique(baseline$id)
n_subj <- length(ids0)

# ---- helper: median (IQR) formatter ----
fmt_median_iqr <- function(x, digits = 2) {
  x <- x[is.finite(x)]
  if (!length(x)) return("-")
  q <- stats::quantile(x, c(.25,.5,.75), na.rm = TRUE)
  sprintf(paste0("%.",digits,"f (%.",digits,"f–%.",digits,"f)"), q[2], q[1], q[3])
}

# ---- numerics @ baseline ----
num_tbl <- map_dfr(numeric_vars, function(v){
  tibble(
    Variable = lab(v),
    Unit     = uni(v),
    Summary  = fmt_median_iqr(baseline[[v]])
  )
})

# ---- binaries EVER (over full data) among baseline subjects ----
bin_max_all <- comb %>%
  group_by(id) %>%
  summarise(across(all_of(binary_vars), ~ as.integer(any(.x %in% c(1, TRUE), na.rm = TRUE))), .groups = "drop")
bin_max_ids0 <- bin_max_all %>% filter(id %in% ids0)

bin_tbl <- map_dfr(binary_vars, function(v){
  s <- sum(bin_max_ids0[[v]], na.rm = TRUE)
  tibble(
    Variable = lab(v),
    Unit     = "",
    Summary  = sprintf("%d (%s)", s, scales::percent(s / n_subj, accuracy = 0.1))
  )
})

# ---- NEW: HD / CRRT / flow-first-CRRT summary ----
flow_summary <- comb %>%
  group_by(id) %>%
  summarise(
    ever_hd   = as.integer(any(hd == 1, na.rm = TRUE)),
    ever_crrt = as.integer(any(crrt == 1, na.rm = TRUE)),
    flow_first_crrt = {
      tmp <- flow[which(crrt == 1)[1]]
      if (length(tmp)) tmp / 1000 else NA_real_
    },
    .groups = "drop"
  )

# summarised overall tallies
n_hd   <- sum(flow_summary$ever_hd, na.rm = TRUE)
n_crrt <- sum(flow_summary$ever_crrt, na.rm = TRUE)

# ---- flow table builder ----
flow_tbl <- tibble(
  Variable = c(
    "Ever on Hemodialysis",
    "Ever on CRRT",
    "Effluent Flow (L/hr) at initial CRRT"
  ),
  Unit     = c("", "", ""),
  Summary  = c(
    sprintf("%d (%s)", n_hd,   scales::percent(n_hd   / n_subj, accuracy = 0.1)),
    sprintf("%d (%s)", n_crrt, scales::percent(n_crrt / n_subj, accuracy = 0.1)),
    fmt_median_iqr(flow_summary$flow_first_crrt, digits = 2)
  )
)

# ---- assemble final table ----
table1_df <- bind_rows(
  num_tbl %>% mutate(.ord = match(Variable, sapply(numeric_vars, lab))),
  bin_tbl %>% mutate(.ord = length(numeric_vars) + match(Variable, sapply(binary_vars, lab))),
  flow_tbl %>% mutate(.ord = 999 + row_number())
) %>%
  arrange(.ord) %>%
  select(Variable, Unit, Summary)

library(gt)
# ---- GT output ----
gt_tab1 <- table1_df |>
  gt() |>
  tab_header(
    title    = md(sprintf("**Table 1.** Baseline characteristics (n = %d)", n_subj)),
    subtitle = md("Continuous variables shown as median (IQR); Categorical as n (%)")
  ) |>
  cols_label(Variable = "Variable", Unit = "Unit", Summary = "Summary") |>
  cols_width(
    Variable ~ px(360),
    Unit     ~ px(150),
    Summary  ~ px(520)
  ) |>
  tab_options(table.font.names = "Arial", table.align = "center")

save_gt <- function(gt_obj, filename, row_count = NULL, subdir = "Manuscript/Tables") {
  # optional footer note
  if (!is.null(row_count)) {
    gt_obj <- gt_obj |>
      tab_source_note(md(sprintf("*Rows: %d*", row_count)))
  }
  
  # make sure subdir exists (relative to getwd())
  if (!dir.exists(subdir)) {
    dir.create(subdir, recursive = TRUE)
  }
  
  # build file paths in that subdir
  html_file <- file.path(subdir, paste0(filename, ".html"))
  png_file  <- file.path(subdir, paste0(filename, ".png"))
  
  # save both formats
  gt::gtsave(gt_obj, html_file)
  gt::gtsave(gt_obj, png_file)
  
  message("✅ Saved table to:\n  ",
          normalizePath(html_file, mustWork = FALSE),
          "\n  ",
          normalizePath(png_file, mustWork = FALSE))
}


save_gt(gt_tab1, "Table1_baseline", row_count = nrow(table1_df))

# Pop parameter summary table ####
run10$final$summary()


library(dplyr)
library(gt)

# Define prospective IDs
dat3_ids <- c(2058, 2060, 2076, 2103, 2126, 9, 15, 20)

# Label meanings for outeq
outeq_labels <- c(
  "1" = "Pre-filter",
  "2" = "Post-filter",
  "3" = "Effluent concentration",
  "4" = "Effluent mg (volume × conc)",
  "5" = "ELF"
)
samples_per_id_outeq <- comb %>%
  filter(!is.na(outeq)) %>%
  group_by(id, outeq) %>%
  summarise(`n(samples)` = sum(!is.na(out)), .groups = "drop") %>%
  mutate(
    Cohort = case_when(
      id %in% dat3_ids            ~ "Prospective (Richly Sampled)",
      id > 100 & id < 1700        ~ "Historical",
      id > 1700                   ~ "Salvaged Blood Samples",
      TRUE                        ~ "Other / Unknown"
    ),
    Method = recode(as.character(outeq), !!!outeq_labels)
  ) %>%
  arrange(match(Cohort,
                c("Prospective (Richly Sampled)",
                  "Historical",
                  "Salvaged Blood Samples")),
          id, outeq)

gt_S1 <- gt(samples_per_id_outeq, groupname_col = "Cohort") %>%
  tab_header(
    title = md("**Table S1.** Number of samples per subject per method (all cohorts)"),
    subtitle = md("Prospective subjects shown in bold; cohorts distinguished by sampling origin.")
  ) %>%
  cols_label(
    id = "ID",
    outeq = "outeq",
    Method = "Sampling method",
    `n(samples)` = "n(samples)"
  ) %>%
  tab_style(
    style = list(cell_text(weight = "bold")),
    locations = cells_body(rows = id %in% dat3_ids)
  ) %>%
  cols_width(
    id ~ px(120),
    outeq ~ px(100),
    Method ~ px(260),
    `n(samples)` ~ px(140)
  ) %>%
  tab_options(table.font.names = "Arial")

save_gt(gt_S1, "TableS1_samples_per_id", row_count = nrow(samples_per_id_outeq))

table2_rollup <- samples_per_id_outeq %>%
  group_by(Cohort, Method) %>%
  summarise(
    `Total samples` = sum(`n(samples)`, na.rm = TRUE),
    `Subjects contributing` = n_distinct(id),
    .groups = "drop"
  ) %>%
  arrange(match(Cohort,
                c("Prospective (Richly Sampled)",
                  "Historical",
                  "Salvaged Blood Samples")),
          Method)

gt_tab2 <- gt(table2_rollup, groupname_col = "Cohort") %>%
  tab_header(
    title = md("**Table 2.** Sample counts by cohort and sampling method"),
    subtitle = md("Summarised totals (not per subject) across prospective, historical, and salvaged cohorts.")
  ) %>%
  cols_label(
    Method = "Sampling method",
    `Subjects contributing` = "Subjects (n)",
    `Total samples` = "Samples (n)"
  ) %>%
  cols_width(
    Method ~ px(260),
    `Subjects contributing` ~ px(160),
    `Total samples` ~ px(160)
  ) %>%
  tab_options(table.font.names = "Arial")

save_gt(gt_tab2, "Table2_samples_by_cohort_method", row_count = nrow(table2_rollup))

# Monte Carlo simulation methods ####

# simulate plasma with various dose regimens and effluent rates
simdat <- PM_data$new(data="Sim/sim3.csv", loq=c(0,0,0,0,0))

names(simdat$standard_data)

sim1 <- PM_sim$new(
  poppar = run10$final,
  data="Sim/sim3.csv",
  model=mod_cef7,
  seed = 12345,
  limits=c(0,1),
  nsim=1000,
  predInt = c(23.9, 48, 0.1)
)

simlabels <- c(
  "2g q12 EI + LD; 60kg",
  "1g q12 II + LD; 60kg",
  "1g q12 EI + LD; 60kg",
  "2g q12 II - LD; 60kg",
  
  "2g q12 EI + LD; 80kg",
  "1g q12 II + LD; 80kg",
  "1g q12 EI + LD; 80kg",
  "2g q12 II - LD; 80kg",
  
  "2g q12 EI + LD; 120kg",
  "1g q12 II + LD; 120kg",
  "1g q12 EI + LD; 120kg",
  "2g q12 II - LD; 120kg"
)

library(plotly)
library(stringr)
library(purrr)

# Shared fonts
panel_title_font <- list(size = 20)
axis_font        <- list(family = "Arial", size = 20)

xlab_common <- list(
  text = "Time",
  bold = TRUE,
  font = axis_font
)

ylab_common <- list(
  text = "Cefepime (mg/L)",
  bold = TRUE,
  font = axis_font
)

# Parse one simlabel like "2g q12 EI + LD; 120kg"
parse_simlabel <- function(label) {
  parts     <- str_split_fixed(label, "\\s*;\\s*", 2)
  reg_part  <- str_squish(parts[, 1])  # "2g q12 EI + LD"
  wt_part   <- str_squish(parts[, 2])  # "120kg"
  
  tokens <- str_split(reg_part, "\\s+")[[1]]
  # tokens = c("2g","q12","EI","+","LD")
  
  list(
    dose     = tokens[1],                            # "2g" / "1g"
    interval = tokens[2],                            # "q12"
    strategy = tokens[3],                            # "EI" / "II"
    ld_flag  = paste(tokens[4:5], collapse = " "),   # "+ LD" / "- LD"
    weight   = wt_part                               # "60kg"/"80kg"/"120kg"
  )
}

# Turn a simlabel into e.g. "2g IV every 12hr over 0.5 hr"
# You can tweak the EI/II → infusion-time mapping here
format_dose_text <- function(label) {
  info <- parse_simlabel(label)
  
  interval_hr  <- gsub("^q", "", info$interval)      # "q12" -> "12"
  interval_txt <- paste0("every ", interval_hr, "hr")
  
  # Hard-coded mapping: II = 0.5 hr, EI = 4 hr
  inf_dur <- ifelse(info$strategy == "EI", "4 hr", "0.5 hr")
  
  paste(info$dose, "IV", interval_txt, "over", inf_dur)
}


make_panel_annotation <- function(letter, compartment, dose_text) {
  list(
    text   = sprintf("<b>%s. %s: %s</b>", letter, compartment, dose_text),
    font   = panel_title_font,
    xref   = "paper",
    yref   = "paper",
    yanchor = "bottom",
    xanchor = "left",
    align  = "left",
    x      = 0.01,
    y      = 1,
    showarrow = FALSE
  )
}

# Build Plasma + ELF panels for a single id.
# simlabels is used to auto-generate the panel title text.
make_panel_pair_from_id <- function(sim_obj,
                                    simlabels,
                                    id,
                                    letters = c("A", "D"),
                                    xlim = c(24, 48)) {
  label     <- simlabels[id]
  dose_text <- format_dose_text(label)
  
  ## Plasma (outeq = 1)
  ann_plasma <- make_panel_annotation(
    letter      = letters[1],
    compartment = "Plasma",
    dose_text   = dose_text
  )
  
  p_plasma <- sim_obj$plot(
    include = id,
    outeq   = 1,
    log     = FALSE,
    xlim    = xlim,
    ci      = 0,
    binSize = 0.1,
    line    = list(color = "dodgerblue", width = 2),
    xlab    = xlab_common,
    ylab    = ylab_common,
    title   = NULL
  )$p %>%
    layout(annotations = ann_plasma,
           xaxis = list(titlefont = axis_font, tickfont = axis_font),
           yaxis = list(titlefont = axis_font, tickfont = axis_font))
  
  ## ELF (outeq = 5)
  ann_elf <- make_panel_annotation(
    letter      = letters[2],
    compartment = "ELF",
    dose_text   = dose_text
  )
  
  p_elf <- sim_obj$plot(
    include = id,
    outeq   = 5,
    log     = FALSE,
    xlim    = xlim,
    ci      = 0,
    binSize = 0.1,
    line    = list(color = "indianred", width = 2),
    xlab    = xlab_common,
    ylab    = ylab_common,
    title   = NULL
  )$p %>%
    layout(annotations = ann_elf,
           xaxis = list(titlefont = axis_font, tickfont = axis_font),
           yaxis = list(titlefont = axis_font, tickfont = axis_font))
  
  list(plasma = p_plasma, elf = p_elf)
}

# ids must be length 3: one per regimen you want shown as columns
# By default:
#   first id  -> A (Plasma), D (ELF)
#   second id -> B (Plasma), E (ELF)
#   third id  -> C (Plasma), F (ELF)
make_6panel_from_ids <- function(sim_obj,
                                 simlabels,
                                 ids,
                                 file = NULL,
                                 xlim = c(24, 48),
                                 width = 4.5 * 400,
                                 height = 3   * 400,
                                 scale = 3) {
  stopifnot(length(ids) == 3)
  
  letter_pairs <- list(c("A", "D"), c("B", "E"), c("C", "F"))
  
  panel_list <- map2(
    ids,
    letter_pairs,
    ~ make_panel_pair_from_id(
      sim_obj   = sim_obj,
      simlabels = simlabels,
      id        = .x,
      letters   = .y,
      xlim      = xlim
    )
  )
  
  combined <- subplot(
    # top row: Plasma
    panel_list[[1]]$plasma,
    panel_list[[2]]$plasma,
    panel_list[[3]]$plasma,
    # bottom row: ELF
    panel_list[[1]]$elf,
    panel_list[[2]]$elf,
    panel_list[[3]]$elf,
    nrows  = 2,
    margin = 0.01,
    titleX = TRUE,
    titleY = TRUE,
    shareX = TRUE,
    shareY = TRUE
  )
  
  if (!is.null(file)) {
    combined %>%
      export_plotly(
        file,
        width  = width,
        height = height,
        scale  = scale
      )
  }
  
  combined
}

combined_plot <- make_6panel_from_ids(
  sim_obj   = sim1,
  simlabels = simlabels,
  ids       = c(12, 9, 11),
  file      = "MEDIAN_SIM_FEP_PLA_ELF_1-2q12_II_EI_120kg.svg"
)

combined_plot_60 <- make_6panel_from_ids(
  sim_obj   = sim1,
  simlabels = simlabels,
  ids       = c(4, 1, 3),   # Ref II, High EI, Low EI at 60kg
  file      = "MEDIAN_SIM_FEP_PLA_ELF_1-2q12_II_EI_60kg.svg"
)


combined_plot_80 <- make_6panel_from_ids(
  sim_obj   = sim1,
  simlabels = simlabels,
  ids       = c(8, 5, 7),   # Ref II, High EI, Low EI at 80kg
  file      = "MEDIAN_SIM_FEP_PLA_ELF_1-2q12_II_EI_80kg.svg"
)

# Number of regimens
n_reg <- 12   # or derive programmatically if you prefer

# Viridis palette for 12 regimens
regimen_colors <- viridisLite::viridis(n_reg, option = "D")  # "D" is default option

marker_style <- list(
  color  = regimen_colors,
  symbol = rep(
    c("circle", "square", "diamond", "cross",
      "triangle-up", "triangle-down", "star", "x"),
    length.out = n_reg
  )
)

# Option 1: colored lines matching markers
line_style <- list(
  color = regimen_colors,
  dash  = rep(
    c("solid", "dot", "dash", "longdash", "dashdot", "longdashdot"),
    length.out = n_reg
  )
)

## PTA output for plasma and ELF compartments ####

# Plasma PTA

target = list(c(0.25, 0.5, 1, 2, 4, 8, 16, 32))

target_type = c("time")

success = c(1)

pta1 <- PM_pta$new(
  simdata=sim1,
  simlabels=simlabels,
  target=target,
  target_type=target_type,
  success = success,
  outeq = 1,
  free_fraction = 0.8
)

pta1$plot(legend=F)

std_reg <- pta1$plot(include=c(4,8,12),legend=list(orientation="h"))

high_reg <- pta1$plot(include=c(1,5,9),legend=list(orientation="h"))

low_reg <- pta1$plot(include=c(3,7,11),legend=list(orientation="h"))

# ELF PTA
target = list(c(0.25, 0.5, 1, 2, 4, 8, 16, 32))

target_type = c("time")

success = c(1)

pta2 <- PM_pta$new(
  simdata=sim1,
  simlabels=simlabels,
  target=target,
  target_type=target_type,
  success = success,
  outeq = 5,
  free_fraction = 0.8
)

## Plotting PTA as a function of MIC ####

## Assume pta1 and pta2 are already created as in your code above

# FONT FOR PANEL TITLES
f <- list(size = 16)

## ---------- PANEL A: Plasma PTA ----------

pta.1 <- pta1$plot(
  ylab   = "Proportion Achieving Goal (%)",
  xlab   = "FEP MIC (mg/L)",
  grid   = TRUE,
  marker = marker_style,
  line   = line_style,
  legend = list(orientation = "h")
)

pta.1 <- plotly::layout(
  pta.1,
  yaxis = list(
    tickvals = seq(0, 1, 0.2),
    ticktext = paste0(seq(0, 100, 20), "%")
  ),
  shapes = list(
    list(
      type = "line",
      x0   = 0,      # in paper coordinates (0–1)
      x1   = 1,
      xref = "paper",
      y0   = 0.9,
      y1   = 0.9,
      yref = "y",
      line = list(color = "black", dash = "dash", width = 2)
    )
  )
)

la <- list(
  text   = "<b>A. Plasma PTA 100% fT>MIC</b>",
  font   = f,
  xref   = "paper",
  yref   = "paper",
  yanchor = "bottom",
  xanchor = "left",
  align  = "left",
  x      = 0.01,
  y      = 1,
  showarrow = FALSE
)

pta.1 <- pta.1 %>% layout(annotations = la)


## ---------- PANEL B: ELF PTA ----------


pta.2 <- pta2$plot(
  ylab   = "Proportion Achieving Goal (%)",
  xlab   = "FEP MIC (mg/L)",
  grid   = TRUE,
  marker = marker_style,
  line   = line_style,
  legend = list(orientation = "h")
)

pta.2 <- plotly::layout(
  pta.2,
  yaxis = list(
    tickvals = seq(0, 1, 0.2),
    ticktext = paste0(seq(0, 100, 20), "%")
  ),
  shapes = list(
    list(
      type = "line",
      x0   = 0,
      x1   = 1,
      xref = "paper",
      y0   = 0.9,
      y1   = 0.9,
      yref = "y",
      line = list(color = "black", dash = "dash", width = 2)
    )
  )
)

lb <- list(
  text   = "<b>B. ELF PTA 100% fT>MIC</b>",
  font   = f,
  xref   = "paper",
  yref   = "paper",
  yanchor = "bottom",
  xanchor = "left",
  align  = "left",
  x      = 0.01,
  y      = 1,
  showarrow = FALSE
)

pta.2 <- pta.2 %>% layout(annotations = lb)

pta.2 <- pta.2 %>% layout(
  legend = list(orientation = "h"),  # Horizontal legend
  showlegend = TRUE                # Ensure legend is shown
)

pta.2 <- pta.2 %>% style(showlegend = FALSE)


## ---------- COMBINED SUBPLOT + LEGEND ----------

# combined_plot <- subplot(
#   pta.1, pta.2,
#   nrows  = 1,
#   margin = 0.05,
#   titleX = TRUE,
#   titleY = TRUE,
#   shareX = TRUE,
#   shareY = TRUE
# ) %>% layout(
#   showlegend = TRUE,
#   legend = list(
#     orientation = "h",
#     x = 0.5,
#     xanchor = "center",
#     y = -0.2
#   )
# )

combined_plot <- subplot(
  pta.1, pta.2,
  nrows  = 1,
  margin = 0.05,
  titleX = TRUE,
  titleY = TRUE,
  shareX = TRUE,
  shareY = TRUE
) %>% layout(
  showlegend = TRUE,
  legend = list(
    orientation = "h",
    x = 0.5,
    xanchor = "center",
    y = -0.2
  ),
  
  xaxis = list(
    type = "log",
    tickvals = c(0.25, 0.5, 1, 2, 4, 8, 16, 32),
    ticktext = c("0.25", "0.5", "1", "2", "4", "8", "16", "32"),
    title = "FEP MIC (mg/L)"
  ),
  
  xaxis2 = list(
    type = "log",
    tickvals = c(0.25, 0.5, 1, 2, 4, 8, 16, 32),
    ticktext = c("0.25", "0.5", "1", "2", "4", "8", "16", "32"),
    title = "FEP MIC (mg/L)"
  ),
  
  yaxis = list(
    tickvals = seq(0, 1, 0.2),
    ticktext = paste0(seq(0, 100, 20), "%"),
    showline = TRUE,
    linecolor = "black",
    linewidth = 1
  ),
  
  yaxis2 = list(
    tickvals = seq(0, 1, 0.2),
    ticktext = paste0(seq(0, 100, 20), "%"),
    showline = TRUE,
    linecolor = "black",
    linewidth = 1
  )
)

combined_plot


# Export (if you have export_plotly defined as in your first script)
combined_plot %>%
  export_plotly("Figure_PTA_FEP_PLA_ELF.svg",
                width = 4 * 300,
                height =2 * 300,
                scale = 2)

## Individual PTA plots ####
library(dplyr)
library(stringr)

# simlabels you defined (order = regimen 1–12)
simlabels <- c(
  "2g q12 EI + LD; 60kg",
  "1g q12 II + LD; 60kg",
  "1g q12 EI + LD; 60kg",
  "2g q12 II - LD; 60kg",
  
  "2g q12 EI + LD; 80kg",
  "1g q12 II + LD; 80kg",
  "1g q12 EI + LD; 80kg",
  "2g q12 II - LD; 80kg",
  
  "2g q12 EI + LD; 120kg",
  "1g q12 II + LD; 120kg",
  "1g q12 EI + LD; 120kg",
  "2g q12 II - LD; 120kg"
)

# Wilson 95% CI for binomial proportion
binom_wilson <- function(p_hat, n, z = 1.96) {
  denom      <- 1 + z^2 / n
  center     <- (p_hat + z^2 / (2 * n)) / denom
  half_width <- z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2)) / denom
  
  tibble(
    PTA_lower = pmax(0, center - half_width),
    PTA_upper = pmin(1, center + half_width)
  )
}

library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

build_pta_df <- function(pta_obj, compartment_name, simlabels) {
  df_raw <- pta_obj$data$data %>%
    mutate(
      RegimenID  = match(label, simlabels),
      Label_new  = label,
      # how many sims & how many successes per row
      n_sims     = map_int(success, length),
      k_success  = map_int(success, sum),
      PTA        = k_success / n_sims
    ) %>%
    # attach Wilson CI
    bind_cols(
      binom_wilson(p_hat = .$PTA, n = .$n_sims)
    )
  
  # Split "<regimen>; <weight>"
  parts    <- str_split_fixed(df_raw$Label_new, "\\s*;\\s*", 2)
  reg_part <- str_squish(parts[, 1])
  wt_part  <- str_squish(parts[, 2])
  
  # "2g q12 EI + LD" -> c("2g","q12","EI","+","LD")
  reg_tokens <- str_split_fixed(reg_part, "\\s+", 5)
  
  df <- df_raw %>%
    transmute(
      target,
      PTA,
      PTA_lower,
      PTA_upper,
      RegimenID,
      Label       = Label_new,
      Dose_amount = reg_tokens[, 1],
      Interval    = reg_tokens[, 2],
      Strategy    = reg_tokens[, 3],
      LD_flag     = paste(reg_tokens[, 4], reg_tokens[, 5]),
      Weight      = wt_part
    ) %>%
    mutate(
      Compartment = compartment_name,
      Strategy    = factor(Strategy, levels = c("II", "EI")),
      Interval    = factor(Interval, levels = c("q12")),
      Weight      = factor(Weight, levels = c("60kg", "80kg", "120kg"))
    ) %>%
    mutate(
      RegimenClass = case_when(
        RegimenID %in% c(4, 8, 12) ~ "Ref II",
        RegimenID %in% c(3, 7, 11) ~ "Low EI",
        RegimenID %in% c(2, 6, 10) ~ "Low II",
        RegimenID %in% c(1, 5, 9)  ~ "High EI",
        TRUE                       ~ NA_character_
      ),
      RegimenClass = factor(RegimenClass,
                            levels = c("Low EI", "Low II", "High EI", "Ref II"))
    ) %>%
    filter(RegimenClass != "Low II")
  
  df
}

# apply to regimens
plasma_df <- build_pta_df(pta1, "Plasma", simlabels)
elf_df    <- build_pta_df(pta2, "ELF",    simlabels)

all_reg <- bind_rows(plasma_df, elf_df) %>%
  filter(!is.na(RegimenClass)) %>%
  mutate(Compartment = factor(Compartment, levels = c("ELF", "Plasma")))

# quick sanity check:
all_reg %>% count(Compartment, Weight, RegimenClass, target)

# sanity check on PTA vs extracted results
# Example: ELF, 80kg, High EI, MIC 8 mg/L
all_reg %>%
  filter(
    Compartment  == "ELF",
    Weight       == "60kg",
    RegimenClass == "Low EI",
    target       == 8
  ) %>%
  select(Label, PTA)

pta2$data$data %>%
  filter(label == "1g q12 EI + LD; 60kg", target == 8) %>%
  select(prop_success)
# both same

all_reg %>%
  filter(
    Compartment  == "Plasma",
    Weight       == "120kg",
    RegimenClass == "Low EI",
    target       == 8
  ) %>%
  select(Label, PTA)


pta1$data$data %>%
  filter(label == "1g q12 EI + LD; 120kg", target == 8) %>%
  select(prop_success)
# both same

regimen_cols <- c(
  "Low EI"  = "seagreen3",
  "High EI" = "indianred",
  "Ref II"  = "dodgerblue"
)

plot_mic_regimen <- function(dat, mic_val) {
  dat %>%
    filter(target == mic_val) %>%
    ggplot(
      aes(
        x    = RegimenClass,
        y    = PTA,
        fill = RegimenClass
      )
    ) +
    geom_col(width = 0.7) +
    geom_errorbar(
      aes(ymin = PTA_lower, ymax = PTA_upper),
      width = 0.2
    ) +
    facet_grid(Compartment ~ Weight) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    scale_fill_manual(values = regimen_cols, drop = FALSE) +
    geom_hline(yintercept = 0.9, linetype = 2) +
    labs(
      title    = paste0("PTA at MIC = ", mic_val, " mg/L"),
      subtitle = "Comparison of Plasma vs ELF by regimen (Low EI / Ref II / High EI)",
      x = "Regimen class",
      y = "Probability of Target Attainment",
      fill = "Regimen"
    ) +
    theme_minimal() +
    theme(
      axis.text.x   = element_text(face = "bold"),
      strip.text    = element_text(face = "bold"),
      panel.border  = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position  = "none"
    )
}

# regenerate
p1 <- plot_mic_regimen(all_reg, 1)
p2 <- plot_mic_regimen(all_reg, 2)
p4 <- plot_mic_regimen(all_reg, 4)
p8 <- plot_mic_regimen(all_reg, 8)

library(patchwork)
main_combined <- (p4 + p8) +
  plot_layout(guides = "collect")

main_combined <- main_combined + plot_annotation(
  title    = "Cefepime PTA by MIC, Compartment, Weight, and Regimen Class",
  subtitle = "Low EI / High EI / Ref II every 12 hr; outcome = 100% fT>MIC"
)

supp_combined <- (p1 + p2) +
  plot_layout(guides = "collect")

supp_combined <- supp_combined + plot_annotation(
  title    = "Cefepime PTA at MIC 1–2 mg/L by Compartment, Weight, and Regimen Class",
  subtitle = "Low EI / High EI / Ref II every 12 hr; outcome = 100% fT>MIC"
)

## ---------- ALL FOUR MIC PLOTS IN ONE LARGE FIGURE ----------

all_mic_combined <- (p1 + p2) /
  (p4 + p8) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title    = "Cefepime PTA by MIC, Compartment, Weight, and Regimen Class",
    subtitle = "Low EI / High EI / Ref II every 12 hr; outcome = 100% fT>MIC"
  )

all_mic_combined

save_combined <- function(plot_obj, file_base,
                          dpi = 1200,
                          w_px = 13 * 1000,
                          h_px = 13 * 500) {
  w_in <- w_px / dpi
  h_in <- h_px / dpi
  
  # SVG
  ggsave(
    filename = paste0(file_base, ".svg"),
    plot     = plot_obj,
    width    = w_in, height = h_in, units = "in", dpi = dpi,
    device   = "svg", bg = "white"
  )
  
  # TIFF
  ggsave(
    filename = paste0(file_base, ".tiff"),
    plot     = plot_obj,
    width    = w_in, height = h_in, units = "in", dpi = dpi,
    device   = "tiff", compression = "lzw", bg = "white"
  )
}

# # Save main and supplemental
# save_combined(main_combined, "Pmetrics_PTA_PLA_ELF_MIC_4_8_regimen")
# save_combined(supp_combined, "Pmetrics_PTA_PLA_ELF_MIC_1_2_regimen")

save_combined(
  all_mic_combined,
  "Pmetrics_PTA_PLA_ELF_MIC_1_2_4_8_regimen",
  w_px = 13 * 1000,
  h_px = 13 * 1000
)

# Monte Carlo simulation methods: sim4 4g/day administration strategy ####

# Analysis goal:
#   Compare the same cefepime maintenance dose across administration strategies:
#     1) II = intermittent infusion: 2g q12h over 0.5 hr
#     2) EI = extended infusion:    2g q12h over 4 hr
#     3) CI = continuous infusion:  2g LD, then 4g/24 hr continuous infusion
#
# Figure:
#   Compare II vs EI vs CI by bodyweight strata.
# 50, 68, 100% targets 5/6/26 - distill down (not for everyone) look at 1 regimen II, EI, CI
# focus on 1 for main manuscript EI cefepime (rest supplement) 
#Flow table will tie together (body weight falls once stratify with flow)
#PLA, ELF BP for (4, 8) 
# Supplemental violin (target - 100%) 80 kg and flow rate (MMopt FEP example) CI vs EI vs II

## Setup ####

library(plotly)
library(stringr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

# Shared fonts
panel_title_font4 <- list(size = 20)
axis_font4        <- list(family = "Arial", size = 20)

xlab_common4 <- list(
  text = "Time",
  bold = TRUE,
  font = axis_font4
)

ylab_common4 <- list(
  text = "Cefepime (mg/L)",
  bold = TRUE,
  font = axis_font4
)


## Simulate plasma and ELF with 4g/day strategy regimens ####

simdat4 <- PM_data$new(data="Sim/sim4.csv", loq=c(0,0,0,0,0))

names(simdat4$standard_data)

sim4 <- PM_sim$new(
  poppar = run10$final,
  data="Sim/sim4.csv",
  model=mod_cef7,
  seed = 12345,
  limits=c(0,1),
  nsim=1000,
  predInt = c(23.9, 48, 0.1)
)


## simlabels4 ####

# simlabels4 order should match sim4.csv.
# All regimens represent 4g/24 hr maintenance dosing.
#
# CI label:
#   "4g q24 CI + LD" = 2g LD followed by 4g/24 hr continuous infusion.
#
# II/EI labels:
#   "2g q12 II + LD" and "2g q12 EI + LD" = 2g q12h maintenance dosing,
#   administered as intermittent or extended infusion.

simlabels4 <- c(
  "2g q12 II + LD; 60kg",
  "2g q12 EI + LD; 60kg",
  "4g q24 CI + LD; 60kg",
  
  "2g q12 II + LD; 80kg",
  "2g q12 EI + LD; 80kg",
  "4g q24 CI + LD; 80kg",
  
  "2g q12 II + LD; 120kg",
  "2g q12 EI + LD; 120kg",
  "4g q24 CI + LD; 120kg"
)


## Helper functions for labels ####

# Parse one sim4 label like "2g q12 EI + LD; 120kg"
parse_simlabel4 <- function(label) {
  parts     <- str_split_fixed(label, "\\s*;\\s*", 2)
  reg_part  <- str_squish(parts[, 1])  # "2g q12 EI + LD"
  wt_part   <- str_squish(parts[, 2])  # "120kg"
  
  tokens <- str_split(reg_part, "\\s+")[[1]]
  # tokens = c("2g","q12","EI","+","LD")
  
  list(
    dose     = tokens[1],                            # "2g" / "4g"
    interval = tokens[2],                            # "q12" / "q24"
    strategy = tokens[3],                            # "II" / "EI" / "CI"
    ld_flag  = paste(tokens[4:5], collapse = " "),   # "+ LD"
    weight   = wt_part                               # "60kg"/"80kg"/"120kg"
  )
}

# Turn a sim4 label into readable panel text
format_dose_text4 <- function(label) {
  info <- parse_simlabel4(label)
  
  interval_hr  <- gsub("^q", "", info$interval)
  interval_txt <- paste0("every ", interval_hr, "hr")
  
  if (info$strategy == "CI") {
    paste("2g LD then", info$dose, "IV continuous infusion over 24 hr")
  } else {
    inf_dur <- ifelse(info$strategy == "EI", "4 hr", "0.5 hr")
    paste("2g LD then", info$dose, "IV", interval_txt, "over", inf_dur)
  }
}

make_panel_annotation4 <- function(letter, compartment, dose_text) {
  list(
    text    = sprintf("<b>%s. %s: %s</b>", letter, compartment, dose_text),
    font    = panel_title_font4,
    xref    = "paper",
    yref    = "paper",
    yanchor = "bottom",
    xanchor = "left",
    align   = "left",
    x       = 0.01,
    y       = 1,
    showarrow = FALSE
  )
}


## Median concentration-time plots ####

# Build Plasma + ELF panels for a single sim4 id.
# simlabels4 is used to auto-generate the panel title text.
make_panel_pair_from_id4 <- function(sim_obj,
                                     simlabels,
                                     id,
                                     letters = c("A", "D"),
                                     xlim = c(24, 48)) {
  label     <- simlabels[id]
  dose_text <- format_dose_text4(label)
  
  ## Plasma (outeq = 1)
  ann_plasma <- make_panel_annotation4(
    letter      = letters[1],
    compartment = "Plasma",
    dose_text   = dose_text
  )
  
  p_plasma <- sim_obj$plot(
    include = id,
    outeq   = 1,
    log     = FALSE,
    xlim    = xlim,
    ci      = 0,
    binSize = 0.1,
    line    = list(color = "dodgerblue", width = 2),
    xlab    = xlab_common4,
    ylab    = ylab_common4,
    title   = NULL
  )$p %>%
    layout(
      annotations = ann_plasma,
      xaxis = list(titlefont = axis_font4, tickfont = axis_font4),
      yaxis = list(titlefont = axis_font4, tickfont = axis_font4)
    )
  
  ## ELF (outeq = 5)
  ann_elf <- make_panel_annotation4(
    letter      = letters[2],
    compartment = "ELF",
    dose_text   = dose_text
  )
  
  p_elf <- sim_obj$plot(
    include = id,
    outeq   = 5,
    log     = FALSE,
    xlim    = xlim,
    ci      = 0,
    binSize = 0.1,
    line    = list(color = "indianred", width = 2),
    xlab    = xlab_common4,
    ylab    = ylab_common4,
    title   = NULL
  )$p %>%
    layout(
      annotations = ann_elf,
      xaxis = list(titlefont = axis_font4, tickfont = axis_font4),
      yaxis = list(titlefont = axis_font4, tickfont = axis_font4)
    )
  
  list(plasma = p_plasma, elf = p_elf)
}

# ids must be length 3: one per regimen you want shown as columns
# By default:
#   first id  -> A (Plasma), D (ELF)
#   second id -> B (Plasma), E (ELF)
#   third id  -> C (Plasma), F (ELF)
make_6panel_from_ids4 <- function(sim_obj,
                                  simlabels,
                                  ids,
                                  file = NULL,
                                  xlim = c(24, 48),
                                  width = 4.5 * 400,
                                  height = 3   * 400,
                                  scale = 3) {
  stopifnot(length(ids) == 3)
  
  letter_pairs <- list(c("A", "D"), c("B", "E"), c("C", "F"))
  
  panel_list <- map2(
    ids,
    letter_pairs,
    ~ make_panel_pair_from_id4(
      sim_obj   = sim_obj,
      simlabels = simlabels,
      id        = .x,
      letters   = .y,
      xlim      = xlim
    )
  )
  
  combined <- subplot(
    # top row: Plasma
    panel_list[[1]]$plasma,
    panel_list[[2]]$plasma,
    panel_list[[3]]$plasma,
    
    # bottom row: ELF
    panel_list[[1]]$elf,
    panel_list[[2]]$elf,
    panel_list[[3]]$elf,
    
    nrows  = 2,
    margin = 0.01,
    titleX = TRUE,
    titleY = TRUE,
    shareX = TRUE,
    shareY = TRUE
  )
  
  if (!is.null(file)) {
    combined %>%
      export_plotly(
        file,
        width  = width,
        height = height,
        scale  = scale
      )
  }
  
  combined
}

combined_plot4_60 <- make_6panel_from_ids4(
  sim_obj   = sim4,
  simlabels = simlabels4,
  ids       = c(1, 2, 3),   # II, EI, CI at 60kg
  file      = "MEDIAN_SIM4_FEP_PLA_ELF_4GDAY_II_EI_CI_60kg.svg"
)

combined_plot4_80 <- make_6panel_from_ids4(
  sim_obj   = sim4,
  simlabels = simlabels4,
  ids       = c(4, 5, 6),   # II, EI, CI at 80kg
  file      = "MEDIAN_SIM4_FEP_PLA_ELF_4GDAY_II_EI_CI_80kg.svg"
)

combined_plot4_120 <- make_6panel_from_ids4(
  sim_obj   = sim4,
  simlabels = simlabels4,
  ids       = c(7, 8, 9),   # II, EI, CI at 120kg
  file      = "MEDIAN_SIM4_FEP_PLA_ELF_4GDAY_II_EI_CI_120kg.svg"
)


## PTA output for plasma and ELF compartments: sim4 ####

# Plasma PTA
target4 = list(c(0.25, 0.5, 1, 2, 4, 8, 16, 32))

target_type4 = c("time")

success4 = c(0.68)

pta4.1 <- PM_pta$new(
  simdata=sim4,
  simlabels=simlabels4,
  target=target4,
  target_type=target_type4,
  success = success4,
  outeq = 1,
  free_fraction = 0.8
)

pta4.1$plot(legend=F)

ii_reg4 <- pta4.1$plot(include=c(1,4,7), legend=list(orientation="h"))

ei_reg4 <- pta4.1$plot(include=c(2,5,8), legend=list(orientation="h"))

ci_reg4 <- pta4.1$plot(include=c(3,6,9), legend=list(orientation="h"))

# ELF PTA
target4 = list(c(0.25, 0.5, 1, 2, 4, 8, 16, 32))

target_type4 = c("time")

success4 = c(0.68)

pta4.2 <- PM_pta$new(
  simdata=sim4,
  simlabels=simlabels4,
  target=target4,
  target_type=target_type4,
  success = success4,
  outeq = 5,
  free_fraction = 0.8
)


## Plotting sim4 PTA as a function of MIC ####

# Number of sim4 regimens
n_reg4 <- 9

# Viridis palette for 9 regimens
regimen_colors4 <- viridisLite::viridis(n_reg4, option = "D")

marker_style4 <- list(
  color  = regimen_colors4,
  symbol = rep(
    c("circle", "square", "diamond", "cross",
      "triangle-up", "triangle-down", "star", "x"),
    length.out = n_reg4
  )
)

# Option 1: colored lines matching markers
line_style4 <- list(
  color = regimen_colors4,
  dash  = rep(
    c("solid", "dot", "dash", "longdash", "dashdot", "longdashdot"),
    length.out = n_reg4
  )
)

# FONT FOR PANEL TITLES
f4 <- list(size = 16)

## ---------- PANEL A: Plasma PTA ----------

pta4.plasma <- pta4.1$plot(
  ylab   = "Proportion Achieving Goal (%)",
  xlab   = "FEP MIC (mg/L)",
  grid   = TRUE,
  marker = marker_style4,
  line   = line_style4,
  legend = list(orientation = "h")
)

pta4.plasma <- plotly::layout(
  pta4.plasma,
  yaxis = list(
    tickvals = seq(0, 1, 0.2),
    ticktext = paste0(seq(0, 100, 20), "%")
  ),
  shapes = list(
    list(
      type = "line",
      x0   = 0,
      x1   = 1,
      xref = "paper",
      y0   = 0.9,
      y1   = 0.9,
      yref = "y",
      line = list(color = "black", dash = "dash", width = 2)
    )
  )
)

la4 <- list(
  text    = "<b>A. Plasma PTA 100% fT>MIC</b>",
  font    = f4,
  xref    = "paper",
  yref    = "paper",
  yanchor = "bottom",
  xanchor = "left",
  align   = "left",
  x       = 0.01,
  y       = 1,
  showarrow = FALSE
)

pta4.plasma <- pta4.plasma %>% layout(annotations = la4)

## ---------- PANEL B: ELF PTA ----------

pta4.elf <- pta4.2$plot(
  ylab   = "Proportion Achieving Goal (%)",
  xlab   = "FEP MIC (mg/L)",
  grid   = TRUE,
  marker = marker_style4,
  line   = line_style4,
  legend = list(orientation = "h")
)

pta4.elf <- plotly::layout(
  pta4.elf,
  yaxis = list(
    tickvals = seq(0, 1, 0.2),
    ticktext = paste0(seq(0, 100, 20), "%")
  ),
  shapes = list(
    list(
      type = "line",
      x0   = 0,
      x1   = 1,
      xref = "paper",
      y0   = 0.9,
      y1   = 0.9,
      yref = "y",
      line = list(color = "black", dash = "dash", width = 2)
    )
  )
)

lb4 <- list(
  text    = "<b>B. ELF PTA 100% fT>MIC</b>",
  font    = f4,
  xref    = "paper",
  yref    = "paper",
  yanchor = "bottom",
  xanchor = "left",
  align   = "left",
  x       = 0.01,
  y       = 1,
  showarrow = FALSE
)

pta4.elf <- pta4.elf %>% layout(annotations = lb4)

pta4.elf <- pta4.elf %>%
  layout(
    legend = list(orientation = "h"),
    showlegend = TRUE
  )

pta4.elf <- pta4.elf %>% style(showlegend = FALSE)

## ---------- COMBINED SUBPLOT + LEGEND ----------

mic_ticks4 <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32)

combined_plot4 <- subplot(
  pta4.plasma, pta4.elf,
  nrows  = 1,
  margin = 0.05,
  titleX = TRUE,
  titleY = TRUE,
  shareX = TRUE,
  shareY = TRUE
) %>%
  layout(
    showlegend = TRUE,
    legend = list(
      orientation = "h",
      x = 0.5,
      xanchor = "center",
      y = -0.2
    ),
    
    ## Panel A: Plasma x-axis
    xaxis = list(
      type = "log",
      tickvals = mic_ticks4,
      ticktext = as.character(mic_ticks4),
      title = "FEP MIC (mg/L)",
      showline = TRUE,
      linecolor = "black",
      linewidth = 1
    ),
    
    ## Panel B: ELF x-axis
    xaxis2 = list(
      type = "log",
      tickvals = mic_ticks4,
      ticktext = as.character(mic_ticks4),
      title = "FEP MIC (mg/L)",
      showline = TRUE,
      linecolor = "black",
      linewidth = 1
    ),
    
    ## Panel A y-axis
    yaxis = list(
      tickvals = seq(0, 1, 0.2),
      ticktext = paste0(seq(0, 100, 20), "%"),
      showline = TRUE,
      linecolor = "black",
      linewidth = 1
    ),
    
    ## Panel B y-axis
    yaxis2 = list(
      tickvals = seq(0, 1, 0.2),
      ticktext = paste0(seq(0, 100, 20), "%"),
      showline = TRUE,
      linecolor = "black",
      linewidth = 1
    )
  )

combined_plot4

combined_plot4 %>%
  export_plotly(
    "Figure_PTA_SIM4_FEP_PLA_ELF_4GDAY_II_EI_CI.svg",
    width = 4 * 300,
    height = 2 * 300,
    scale = 2
  )


## Individual PTA plots: sim4 ####

# Wilson 95% CI for binomial proportion
binom_wilson4 <- function(p_hat, n, z = 1.96) {
  denom      <- 1 + z^2 / n
  center     <- (p_hat + z^2 / (2 * n)) / denom
  half_width <- z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2)) / denom
  
  tibble(
    PTA_lower = pmax(0, center - half_width),
    PTA_upper = pmin(1, center + half_width)
  )
}

build_pta_df4 <- function(pta_obj, compartment_name, simlabels) {
  df_raw <- pta_obj$data$data %>%
    mutate(
      RegimenID  = match(label, simlabels),
      Label_new  = label,
      # how many sims & how many successes per row
      n_sims     = map_int(success, length),
      k_success  = map_int(success, sum),
      PTA        = k_success / n_sims
    ) %>%
    # attach Wilson CI
    bind_cols(
      binom_wilson4(p_hat = .$PTA, n = .$n_sims)
    )
  
  # Split "<regimen>; <weight>"
  parts    <- str_split_fixed(df_raw$Label_new, "\\s*;\\s*", 2)
  reg_part <- str_squish(parts[, 1])
  wt_part  <- str_squish(parts[, 2])
  
  # "2g q12 EI + LD" -> c("2g","q12","EI","+","LD")
  reg_tokens <- str_split_fixed(reg_part, "\\s+", 5)
  
  df <- df_raw %>%
    transmute(
      target,
      PTA,
      PTA_lower,
      PTA_upper,
      RegimenID,
      Label       = Label_new,
      Dose_amount = reg_tokens[, 1],
      Interval    = reg_tokens[, 2],
      Strategy    = reg_tokens[, 3],
      LD_flag     = paste(reg_tokens[, 4], reg_tokens[, 5]),
      Weight      = wt_part
    ) %>%
    mutate(
      Compartment = compartment_name,
      Strategy    = factor(Strategy, levels = c("II", "EI", "CI")),
      Interval    = factor(Interval, levels = c("q12", "q24")),
      Weight      = factor(Weight, levels = c("60kg", "80kg", "120kg"))
    ) %>%
    mutate(
      RegimenClass = case_when(
        RegimenID %in% c(1, 4, 7) ~ "II",
        RegimenID %in% c(2, 5, 8) ~ "EI",
        RegimenID %in% c(3, 6, 9) ~ "CI",
        TRUE                      ~ NA_character_
      ),
      RegimenClass = factor(
        RegimenClass,
        levels = c("II", "EI", "CI")
      )
    )
  
  df
}

# apply to regimens
plasma_df4 <- build_pta_df4(pta4.1, "Plasma", simlabels4)
elf_df4    <- build_pta_df4(pta4.2, "ELF",    simlabels4)

all_reg4 <- bind_rows(plasma_df4, elf_df4) %>%
  filter(!is.na(RegimenClass)) %>%
  mutate(Compartment = factor(Compartment, levels = c("ELF", "Plasma")))

# quick sanity check:
all_reg4 %>% count(Compartment, Weight, RegimenClass, target)

# sanity check on PTA vs extracted results
# Example: ELF, 80kg, EI, MIC 8 mg/L
all_reg4 %>%
  filter(
    Compartment  == "ELF",
    Weight       == "80kg",
    RegimenClass == "EI",
    target       == 8
  ) %>%
  select(Label, PTA)

pta4.2$data$data %>%
  filter(label == "2g q12 EI + LD; 80kg", target == 8) %>%
  select(prop_success)
# both same

all_reg4 %>%
  filter(
    Compartment  == "Plasma",
    Weight       == "120kg",
    RegimenClass == "CI",
    target       == 8
  ) %>%
  select(Label, PTA)

pta4.1$data$data %>%
  filter(label == "4g q24 CI + LD; 120kg", target == 8) %>%
  select(prop_success)
# both same

regimen_cols4 <- c(
  "II" = "dodgerblue",
  "EI" = "seagreen3",
  "CI" = "indianred"
)

plot_mic_regimen4 <- function(dat, mic_val) {
  dat %>%
    filter(target == mic_val) %>%
    ggplot(
      aes(
        x    = RegimenClass,
        y    = PTA,
        fill = RegimenClass
      )
    ) +
    geom_col(width = 0.7) +
    geom_errorbar(
      aes(ymin = PTA_lower, ymax = PTA_upper),
      width = 0.2
    ) +
    facet_grid(Compartment ~ Weight) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    scale_fill_manual(values = regimen_cols4, drop = FALSE) +
    geom_hline(yintercept = 0.9, linetype = 2) +
    labs(
      title    = paste0("PTA at MIC = ", mic_val, " mg/L"),
      subtitle = "4g/24 hr maintenance dosing",
      x = "Dosing strategy",
      y = "Probability of Target Attainment",
      fill = "Strategy"
    ) +
    theme_minimal() +
    theme(
      axis.text.x      = element_text(face = "bold"),
      strip.text       = element_text(face = "bold"),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position  = "none"
    )
}

# regenerate
p1_4 <- plot_mic_regimen4(all_reg4, 1)
p2_4 <- plot_mic_regimen4(all_reg4, 2)
p4_4 <- plot_mic_regimen4(all_reg4, 4)
p8_4 <- plot_mic_regimen4(all_reg4, 8)

main_combined4 <- patchwork::wrap_plots(
  p4_4,
  p8_4,
  nrow = 1,
  guides = "collect"
)

main_combined4 <- main_combined4 + patchwork::plot_annotation(
  title    = "Cefepime PTA by MIC, Compartment, Weight, and Dosing Strategy",
  subtitle = "II vs EI vs CI; 4g/24 hr maintenance dosing; outcome = 100% fT>MIC"
)

main_combined4

supp_combined4 <- patchwork::wrap_plots(
  p1_4,
  p2_4,
  nrow = 1,
  guides = "collect"
)

supp_combined4 <- supp_combined4 + patchwork::plot_annotation(
  title    = "Cefepime PTA at MIC 1-2 mg/L by Compartment, Weight, and Dosing Strategy",
  subtitle = "II vs EI vs CI; 4g/24 hr maintenance dosing; outcome = 100% fT>MIC"
)

supp_combined4

## Combine all four MIC plots into one large figure ####

all_mic_combined4 <- patchwork::wrap_plots(
  p1_4,
  p2_4,
  p4_4,
  p8_4,
  ncol = 2,
  nrow = 2,
  guides = "collect"
)

all_mic_combined4 <- all_mic_combined4 + patchwork::plot_annotation(
  title    = "Cefepime PTA by MIC, Compartment, Weight, and Dosing Strategy",
  subtitle = "II vs EI vs CI; 4g/24 hr maintenance dosing; outcome = 100% fT>MIC",
  tag_levels = "A"
)

all_mic_combined4


## Save figures ####

save_combined <- function(plot_obj, file_base,
                          dpi = 1200,
                          w_px = 13 * 1000,
                          h_px = 13 * 500) {
  w_in <- w_px / dpi
  h_in <- h_px / dpi
  
  # SVG
  ggsave(
    filename = paste0(file_base, ".svg"),
    plot     = plot_obj,
    width    = w_in,
    height   = h_in,
    units    = "in",
    dpi      = dpi,
    device   = "svg",
    bg       = "white"
  )
  
  # TIFF
  ggsave(
    filename    = paste0(file_base, ".tiff"),
    plot        = plot_obj,
    width       = w_in,
    height      = h_in,
    units       = "in",
    dpi         = dpi,
    device      = "tiff",
    compression = "lzw",
    bg          = "white"
  )
}

# save_combined(main_combined4, "Pmetrics_SIM4_PTA_PLA_ELF_MIC_4_8_4GDAY_II_EI_CI")
# save_combined(supp_combined4, "Pmetrics_SIM4_PTA_PLA_ELF_MIC_1_2_4GDAY_II_EI_CI")

save_combined(
  all_mic_combined4,
  "Pmetrics_SIM4_PTA_PLA_ELF_MIC_1_2_4_8_4GDAY_II_EI_CI",
  h_px = 13 * 1000
)

# Monte Carlo simulation methods: sim4 flow-rate PTA table ####

# Separate analysis goal:
#   Evaluate PTA across CRRT effluent flow rates while keeping the same
#   dosing strategy comparison: II vs EI vs CI at 4g/24 hr maintenance dosing.
#
# Recommended separate simulation file:
#   Sim/sim4_flow.csv
#
# The flow file should contain the same II/EI/CI regimen structure across:
#   bodyweight strata: 60kg, 80kg, 120kg
#   flow rates: 1, 2, 3, 4 L/hr
#
# This creates 36 simulated profiles:
#   3 dosing strategies x 3 bodyweight strata x 4 flow rates


## Simulate plasma and ELF with flow-rate scenarios ####


simdat4_flow <- PM_data$new(data="Sim/sim4_flow.csv", loq=c(0,0,0,0,0))

names(simdat4_flow$standard_data)

sim4_flow <- PM_sim$new(
  poppar = run10$final,
  data="Sim/sim4_flow.csv",
  model=mod_cef7,
  seed = 12345,
  limits=c(0,1),
  nsim=1000,
  predInt = c(23.9, 48, 0.1)
)


## simlabels4_flow ####

simlabels4_flow <- c(
  # FLOW 1 L/hr
  "2g q12 II + LD; 60kg; FLOW 1L/hr",
  "2g q12 EI + LD; 60kg; FLOW 1L/hr",
  "4g q24 CI + LD; 60kg; FLOW 1L/hr",
  "2g q12 II + LD; 80kg; FLOW 1L/hr",
  "2g q12 EI + LD; 80kg; FLOW 1L/hr",
  "4g q24 CI + LD; 80kg; FLOW 1L/hr",
  "2g q12 II + LD; 120kg; FLOW 1L/hr",
  "2g q12 EI + LD; 120kg; FLOW 1L/hr",
  "4g q24 CI + LD; 120kg; FLOW 1L/hr",
  
  # FLOW 2 L/hr
  "2g q12 II + LD; 60kg; FLOW 2L/hr",
  "2g q12 EI + LD; 60kg; FLOW 2L/hr",
  "4g q24 CI + LD; 60kg; FLOW 2L/hr",
  "2g q12 II + LD; 80kg; FLOW 2L/hr",
  "2g q12 EI + LD; 80kg; FLOW 2L/hr",
  "4g q24 CI + LD; 80kg; FLOW 2L/hr",
  "2g q12 II + LD; 120kg; FLOW 2L/hr",
  "2g q12 EI + LD; 120kg; FLOW 2L/hr",
  "4g q24 CI + LD; 120kg; FLOW 2L/hr",
  
  # FLOW 3 L/hr
  "2g q12 II + LD; 60kg; FLOW 3L/hr",
  "2g q12 EI + LD; 60kg; FLOW 3L/hr",
  "4g q24 CI + LD; 60kg; FLOW 3L/hr",
  "2g q12 II + LD; 80kg; FLOW 3L/hr",
  "2g q12 EI + LD; 80kg; FLOW 3L/hr",
  "4g q24 CI + LD; 80kg; FLOW 3L/hr",
  "2g q12 II + LD; 120kg; FLOW 3L/hr",
  "2g q12 EI + LD; 120kg; FLOW 3L/hr",
  "4g q24 CI + LD; 120kg; FLOW 3L/hr",
  
  # FLOW 4 L/hr
  "2g q12 II + LD; 60kg; FLOW 4L/hr",
  "2g q12 EI + LD; 60kg; FLOW 4L/hr",
  "4g q24 CI + LD; 60kg; FLOW 4L/hr",
  "2g q12 II + LD; 80kg; FLOW 4L/hr",
  "2g q12 EI + LD; 80kg; FLOW 4L/hr",
  "4g q24 CI + LD; 80kg; FLOW 4L/hr",
  "2g q12 II + LD; 120kg; FLOW 4L/hr",
  "2g q12 EI + LD; 120kg; FLOW 4L/hr",
  "4g q24 CI + LD; 120kg; FLOW 4L/hr"
)


## PTA output for flow-rate analysis ####

target4_flow = list(c(0.25, 0.5, 1, 2, 4, 8, 16, 32))

target_type4_flow = c("time")

success4_flow = c(1)

# Plasma PTA
pta4_flow.1 <- PM_pta$new(
  simdata=sim4_flow,
  simlabels=simlabels4_flow,
  target=target4_flow,
  target_type=target_type4_flow,
  success = success4_flow,
  outeq = 1,
  free_fraction = 0.8
)

# ELF PTA
pta4_flow.2 <- PM_pta$new(
  simdata=sim4_flow,
  simlabels=simlabels4_flow,
  target=target4_flow,
  target_type=target_type4_flow,
  success = success4_flow,
  outeq = 5,
  free_fraction = 0.8
)


## Build flow-rate PTA dataframe ####

parse_simlabel4_flow <- function(label) {
  parts     <- str_split_fixed(label, "\\s*;\\s*", 3)
  reg_part  <- str_squish(parts[, 1])
  wt_part   <- str_squish(parts[, 2])
  flow_part <- str_squish(parts[, 3])
  
  tokens <- str_split(reg_part, "\\s+")[[1]]
  
  list(
    dose     = tokens[1],
    interval = tokens[2],
    strategy = tokens[3],
    ld_flag  = paste(tokens[4:5], collapse = " "),
    weight   = wt_part,
    flow     = str_remove(flow_part, "^FLOW\\s+")
  )
}

build_pta_df4_flow <- function(pta_obj, compartment_name, simlabels) {
  df_raw <- pta_obj$data$data %>%
    mutate(
      RegimenID  = match(label, simlabels),
      Label_new  = label,
      # how many sims & how many successes per row
      n_sims     = map_int(success, length),
      k_success  = map_int(success, sum),
      PTA        = k_success / n_sims
    ) %>%
    # attach Wilson CI
    bind_cols(
      binom_wilson4(p_hat = .$PTA, n = .$n_sims)
    )
  
  #We want weight column to go away
  
  # Split "<regimen>; <weight>; <flow>"
  parts     <- str_split_fixed(df_raw$Label_new, "\\s*;\\s*", 3)
  reg_part  <- str_squish(parts[, 1])
  wt_part   <- str_squish(parts[, 2])
  flow_part <- str_squish(parts[, 3])
  
  # "2g q12 EI + LD" -> c("2g","q12","EI","+","LD")
  reg_tokens <- str_split_fixed(reg_part, "\\s+", 5)
  
  df <- df_raw %>%
    transmute(
      target,
      PTA,
      PTA_lower,
      PTA_upper,
      RegimenID,
      Label       = Label_new,
      Dose_amount = reg_tokens[, 1],
      Interval    = reg_tokens[, 2],
      Strategy    = reg_tokens[, 3],
      LD_flag     = paste(reg_tokens[, 4], reg_tokens[, 5]),
      Weight      = wt_part,
      Flow        = str_remove(flow_part, "^FLOW\\s+")
    ) %>%
    mutate(
      Compartment = compartment_name,
      Strategy    = factor(Strategy, levels = c("II", "EI", "CI")),
      Interval    = factor(Interval, levels = c("q12", "q24")),
      Weight      = factor(Weight, levels = c("60kg", "80kg", "120kg")),
      Flow        = factor(Flow, levels = c("1L/hr", "2L/hr", "3L/hr", "4L/hr"))
    ) %>%
    mutate(
      RegimenClass = case_when(
        Strategy == "II" ~ "II",
        Strategy == "EI" ~ "EI",
        Strategy == "CI" ~ "CI",
        TRUE             ~ NA_character_
      ),
      RegimenClass = factor(
        RegimenClass,
        levels = c("II", "EI", "CI")
      )
    )
  
  df
}

# apply to flow-rate regimens
plasma_df4_flow <- build_pta_df4_flow(pta4_flow.1, "Plasma", simlabels4_flow)
elf_df4_flow    <- build_pta_df4_flow(pta4_flow.2, "ELF",    simlabels4_flow)

all_reg4_flow <- bind_rows(plasma_df4_flow, elf_df4_flow) %>%
  filter(!is.na(RegimenClass)) %>%
  mutate(
    Compartment  = factor(Compartment, levels = c("ELF", "Plasma")),
    RegimenClass = factor(RegimenClass, levels = c("II", "EI", "CI")),
    Weight       = factor(Weight, levels = c("60kg", "80kg", "120kg")),
    Flow         = factor(Flow, levels = c("1L/hr", "2L/hr", "3L/hr", "4L/hr"))
  )

# quick sanity check:
all_reg4_flow %>% count(Compartment, Weight, Flow, RegimenClass, target)


## Flow-rate PTA table ####

# This table summarizes PTA across CRRT effluent flow rates of 1-4 L/hr.
# Use target %in% c(4, 8) for a concise main table, or change to any MIC(s)
# of interest.

flow_pta_table4 <- all_reg4_flow %>%
  filter(target %in% c(4, 8)) %>%
  mutate(
    MIC_mg_L          = target,
    PTA_percent       = round(PTA * 100, 1),
    PTA_lower_percent = round(PTA_lower * 100, 1),
    PTA_upper_percent = round(PTA_upper * 100, 1),
    PTA_95CI = paste0(
      PTA_percent,
      " (",
      PTA_lower_percent,
      "-",
      PTA_upper_percent,
      ")"
    )
  ) %>%
  select(
    Compartment,
    MIC_mg_L,
    Weight,
    Flow,
    RegimenClass,
    PTA_95CI
  ) %>%
  pivot_wider(
    names_from  = RegimenClass,
    values_from = PTA_95CI
  ) %>%
  arrange(
    Compartment,
    MIC_mg_L,
    Weight,
    Flow
  )

flow_pta_table4

write.csv(
  flow_pta_table4,
  "Pmetrics_SIM4_PTA_flow_table_MIC_4_8_4GDAY_II_EI_CI.csv",
  row.names = FALSE
)


library(dplyr)
library(readr)
library(stringr)
library(gt)

flow_tab <- read_csv(
  "Pmetrics_SIM4_PTA_flow_table_MIC_4_8_4GDAY_II_EI_CI.csv",
  show_col_types = FALSE
) %>%
  mutate(
    Compartment = factor(Compartment, levels = c("Plasma", "ELF")),
    MIC_mg_L    = factor(MIC_mg_L, levels = c(4, 8)),
    Weight      = factor(Weight, levels = c("60kg", "80kg", "120kg")),
    Flow        = factor(Flow, levels = c("1L/hr", "2L/hr", "3L/hr", "4L/hr")),
    
    # numeric versions for conditional formatting
    II_num = as.numeric(str_extract(II, "^[0-9.]+")),
    EI_num = as.numeric(str_extract(EI, "^[0-9.]+")),
    CI_num = as.numeric(str_extract(CI, "^[0-9.]+"))
  ) %>%
  arrange(Compartment, MIC_mg_L, Weight, Flow)

make_flow_gt <- function(dat, compartment_name) {
  
  dat_sub <- dat %>%
    filter(Compartment == compartment_name) %>%
    select(MIC_mg_L, Weight, Flow, II, EI, CI, II_num, EI_num, CI_num)
  
  gt_tab <- dat_sub %>%
    gt(groupname_col = "MIC_mg_L") %>%
    
    cols_label(
      MIC_mg_L = "MIC (mg/L)",
      Weight   = "Weight",
      Flow     = md("CRRT effluent flow"),
      II       = md("**2 g q12 II + LD**"),
      EI       = md("**2 g q12 EI + LD**"),
      CI       = md("**4 g q24 CI + LD**")
    ) %>%
    
    tab_spanner(
      label = md("**PTA, % (95% CI)**"),
      columns = c(II, EI, CI)
    ) %>%
    
    tab_header(
      title = md(paste0("**", compartment_name, " PTA across CRRT effluent flow rates**")),
      subtitle = md("Target = 100% fT>MIC")
    ) %>%
    
    cols_align(
      align = "center",
      columns = c(Weight, Flow, II, EI, CI)
    ) %>%
    
    fmt_missing(
      columns = everything(),
      missing_text = "—"
    ) %>%
    
    # Bold values meeting the typical PTA threshold
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(columns = II, rows = II_num >= 90)
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(columns = EI, rows = EI_num >= 90)
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(columns = CI, rows = CI_num >= 90)
    ) %>%
    
    cols_hide(columns = c(II_num, EI_num, CI_num)) %>%
    
    tab_source_note(
      source_note = md("Values are presented as **PTA % (95% Wilson CI)**. **Bold** indicates PTA ≥ 90%.")
    ) %>%
    
    opt_row_striping() %>%
    
    tab_options(
      table.font.size = 11,
      heading.align = "left",
      data_row.padding = px(4),
      row_group.font.weight = "bold",
      column_labels.font.weight = "bold"
    )
  
  gt_tab
}

tab_plasma <- make_flow_gt(flow_tab, "Plasma")
tab_elf    <- make_flow_gt(flow_tab, "ELF")

tab_plasma
tab_elf

gtsave(tab_plasma, "Table_SIM4_flow_Plasma.html")
gtsave(tab_elf,    "Table_SIM4_flow_ELF.html")

gtsave(tab_plasma, "Table_SIM4_flow_Plasma.rtf")
gtsave(tab_elf,    "Table_SIM4_flow_ELF.rtf")

## Optional: flow-rate table function ####

make_flow_pta_table4 <- function(dat,
                                 mic_vals = c(4, 8),
                                 compartment_vals = c("ELF", "Plasma")) {
  dat %>%
    filter(
      target %in% mic_vals,
      Compartment %in% compartment_vals
    ) %>%
    mutate(
      MIC_mg_L          = target,
      PTA_percent       = round(PTA * 100, 1),
      PTA_lower_percent = round(PTA_lower * 100, 1),
      PTA_upper_percent = round(PTA_upper * 100, 1),
      PTA_95CI = paste0(
        PTA_percent,
        " (",
        PTA_lower_percent,
        "-",
        PTA_upper_percent,
        ")"
      )
    ) %>%
    select(
      Compartment,
      MIC_mg_L,
      Weight,
      Flow,
      RegimenClass,
      PTA_95CI
    ) %>%
    pivot_wider(
      names_from  = RegimenClass,
      values_from = PTA_95CI
    ) %>%
    arrange(
      Compartment,
      MIC_mg_L,
      Weight,
      Flow
    )
}




## Let's now graph the PDI distributions as a violin ####

library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)

build_pdi_df <- function(pta_obj, compartment_name, simlabels) {
  df_raw <- pta_obj$data$data %>%
    mutate(
      RegimenID = match(label, simlabels),
      Label_new = label
    )
  
  # Split "<regimen>; <weight>"
  parts    <- stringr::str_split_fixed(df_raw$Label_new, "\\s*;\\s*", 2)
  reg_part <- str_squish(parts[, 1])
  wt_part  <- str_squish(parts[, 2])
  
  # "2g q12 EI + LD" -> c("2g","q12","EI","+","LD")
  reg_tokens <- stringr::str_split_fixed(reg_part, "\\s+", 5)
  
  df <- df_raw %>%
    transmute(
      target,
      RegimenID,
      Label       = Label_new,
      Dose_amount = reg_tokens[, 1],
      Interval    = reg_tokens[, 2],
      Strategy    = reg_tokens[, 3],
      LD_flag     = paste(reg_tokens[, 4], reg_tokens[, 5]),
      Weight      = wt_part,
      pdi,
      success
    ) %>%
    # unnest sim-level PDI and success
    tidyr::unnest(cols = c(pdi, success)) %>%
    group_by(target, RegimenID, Label, Dose_amount, Interval,
             Strategy, LD_flag, Weight) %>%
    mutate(sim_id = row_number()) %>%  # 1..nsim per regimen/MIC
    ungroup() %>%
    mutate(
      Compartment = compartment_name,
      Strategy    = factor(Strategy, levels = c("II", "EI")),
      Interval    = factor(Interval, levels = c("q12")),
      Weight      = factor(Weight, levels = c("60kg", "80kg", "120kg"))
    ) %>%
    mutate(
      RegimenClass = case_when(
        RegimenID %in% c(4, 8, 12) ~ "Ref II",
        RegimenID %in% c(3, 7, 11) ~ "Low EI",
        RegimenID %in% c(2, 6, 10) ~ "Low II",
        RegimenID %in% c(1, 5, 9)  ~ "High EI",
        TRUE                       ~ NA_character_
      ),
      RegimenClass = factor(RegimenClass,
                            levels = c("Low EI", "Low II", "High EI", "Ref II"))
    ) %>%
    filter(RegimenClass != "Low II")
  
  df
}

# Build PDI datasets for Plasma and ELF
plasma_pdi_df <- build_pdi_df(pta1, "Plasma", simlabels)
elf_pdi_df    <- build_pdi_df(pta2, "ELF",    simlabels)

all_pdi <- bind_rows(plasma_pdi_df, elf_pdi_df) %>%
  filter(!is.na(RegimenClass)) %>%
  mutate(
    Compartment = factor(Compartment, levels = c("ELF", "Plasma"))
  )

#4/1/26 AV 
pdi_summary <- all_pdi %>%
  group_by(target, Compartment, Weight, RegimenClass) %>%
  summarise(
    pdi_median = median(pdi, na.rm = TRUE),
    pdi_q1     = quantile(pdi, 0.25, na.rm = TRUE),
    pdi_q3     = quantile(pdi, 0.75, na.rm = TRUE),
    pdi_mean   = mean(pdi, na.rm = TRUE),
    .groups = "drop"
  )

View(pdi_summary)


regimen_cols <- c(
  "Low EI"  = "seagreen3",
  "High EI" = "indianred",
  "Ref II"  = "dodgerblue"
)

plot_pdi_violin <- function(dat, mic_val) {
  dat %>%
    filter(target == mic_val) %>%
    ggplot(
      aes(
        x    = RegimenClass,
        y    = pdi,
        fill = RegimenClass
      )
    ) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    # median PDI per group
    stat_summary(
      fun   = median,
      geom  = "point",
      size  = 1.8,
      colour = "black"
    ) +
    facet_grid(Compartment ~ Weight) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.2)
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = regimen_cols, drop = FALSE) +
    geom_hline(yintercept = 1, linetype = 2) +
    labs(
      title    = paste0("PDI distribution at MIC = ", mic_val, " mg/L"),
      subtitle = "Violin = PDI; point = median; dashed line = target (PDI = 1)",
      x = "Regimen class",
      y = "PDI (fraction of interval with C > MIC)",
      fill = "Regimen"
    ) +
    theme_minimal() +
    theme(
      axis.text.x   = element_text(face = "bold"),
      strip.text    = element_text(face = "bold"),
      panel.border  = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position  = "none"
    )
}

# MIC 8, 4, 1, 2 similarly
pdi_v8 <- plot_pdi_violin(all_pdi, 8)
pdi_v4 <- plot_pdi_violin(all_pdi, 4)
pdi_v2 <- plot_pdi_violin(all_pdi, 2)
pdi_v1 <- plot_pdi_violin(all_pdi, 1)

# Ultimate figure, show probs and pdi distributions that underpin them

library(patchwork)

save_combined <- function(plot_obj, file_base,
                          dpi = 1200,
                          w_px = 13 * 1000,
                          h_px = 13 * 1000) {
  w_in <- w_px / dpi
  h_in <- h_px / dpi
  
  # SVG
  ggsave(
    filename = paste0(file_base, ".svg"),
    plot     = plot_obj,
    width    = w_in, height = h_in, units = "in", dpi = dpi,
    device   = "svg", bg = "white"
  )
  
  # TIFF
  ggsave(
    filename = paste0(file_base, ".tiff"),
    plot     = plot_obj,
    width    = w_in, height = h_in, units = "in", dpi = dpi,
    device   = "tiff", compression = "lzw", bg = "white"
  )
}

# MIC = 8
pta_p8    <- plot_mic_regimen(all_reg, 8)           # your bar plot
pdi_v8    <- plot_pdi_violin(all_pdi, 8)            # or all_pdi_clamped if you made it

combined_p8 <- pta_p8 / pdi_v8 +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    title    = "Cefepime PTA and PDI at MIC = 8 mg/L",
    subtitle = "Top: PTA by regimen / weight / compartment (error bars = 95% binomial CI)\n" %+%
      "Bottom: PDI distributions (violin = PDI, dot = median, dashed line = target PDI = 1)",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11)
    )
  )

combined_p8

# MIC = 4
pta_p4    <- plot_mic_regimen(all_reg, 4)           # your bar plot
pdi_v4    <- plot_pdi_violin(all_pdi, 4)            # or all_pdi_clamped if you made it

combined_p4 <- pta_p4 / pdi_v4 +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    title    = "Cefepime PTA and PDI at MIC = 4 mg/L",
    subtitle = "Top: PTA by regimen / weight / compartment (error bars = 95% binomial CI)\n" %+%
      "Bottom: PDI distributions (violin = PDI, dot = median, dashed line = target PDI = 1)",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11)
    )
  )

combined_p4

# MIC = 2
pta_p2    <- plot_mic_regimen(all_reg, 2)           # your bar plot
pdi_v2    <- plot_pdi_violin(all_pdi, 2)            # or all_pdi_clamped if you made it

combined_p2 <- pta_p2 / pdi_v2 +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(
    title    = "Cefepime PTA and PDI at MIC = 2 mg/L",
    subtitle = "Top: PTA by regimen / weight / compartment (error bars = 95% binomial CI)\n" %+%
      "Bottom: PDI distributions (violin = PDI, dot = median, dashed line = target PDI = 1)",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 11)
    )
  )

combined_p2

# Save the plots
save_combined(combined_p8, "Pmetrics_PTA_PDI_PLA_ELF_MIC_8_regimen")
save_combined(combined_p4, "Pmetrics_PTA_PDI_PLA_ELF_MIC_4_regimen")
save_combined(combined_p2, "Pmetrics_PTA_PDI_PLA_ELF_MIC_2_regimen")

# Model schematic plotting ####
library(DiagrammeR)
library(DiagrammeR)
library(DiagrammeRsvg)
library(magick)

save_grviz <- function(grviz_obj,
                       filename,
                       out_dir  = "Manuscript/Figures",
                       dpi      = 800,
                       width_in = 5) {   # physical width for the TIFF only
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  svg_file  <- file.path(out_dir, paste0(filename, ".svg"))
  #tiff_file <- file.path(out_dir, paste0(filename, ".tiff"))
  
  ## --- SVG (vector, no dpi) ---
  svg_txt <- DiagrammeRsvg::export_svg(grviz_obj)
  writeLines(svg_txt, svg_file)
  
  # ## --- TIFF (raster, dpi matters) ---
  # img <- magick::image_read_svg(charToRaw(svg_txt))
  # 
  # # scale to requested physical width at given dpi
  # if (!is.null(width_in) && !is.null(dpi)) {
  #   width_px <- as.integer(width_in * dpi)
  #   img <- magick::image_scale(img, paste0(width_px))
  # }
  # 
  # magick::image_write(
  #   img,
  #   path       = tiff_file,
  #   format     = "tiff",
  #   compression = "lzw",
  #   density    = paste0(dpi, "x", dpi)  # dpi only relevant here
  # )
  # 
  message(
    "Saved figure to:\n  ",
    normalizePath(svg_file,  mustWork = FALSE)#, "\n  ",
    #normalizePath(tiff_file, mustWork = FALSE)
  )
  
  invisible(list(svg = svg_file#, 
                 #     tiff = tiff_file
  )
  )
}

library(DiagrammeR)

viz_cef_struct_minimal <- grViz("
  digraph cef_struct_minimal {

    graph [
      layout  = neato
      overlap = false
      splines = true
    ]

    node [
      shape    = box
      fontsize = 12
      fontname = 'Arial'
      style    = 'filled'
      penwidth = 1.1
    ]

    edge [
      fontsize   = 10
      fontname   = 'Arial'
      arrowsize  = 0.7
      penwidth   = 1
    ]

    # ---- compartments (fixed positions) ----
    dose [label='Dose\\n(IV infusion)',
          style='rounded,filled',
          fillcolor='#D3ECE4',
          pos='0,2.6!']

    central [label='Pre-filter\\nplasma\\nObs: C_pre',
             fillcolor='#FADADD',
             pos='0,0!']

    peripheral [label='Peripheral\\n(no Obs)',
                fillcolor='#D4E4E8',
                pos='-2.6,1.6!']

    elf [label='ELF\\nObs: C_elf',
         fillcolor='#D4E4E8',
         pos='-2.6,-1.6!']

    post [label='Post-filter\\nplasma\\nObs: C_post',
          fillcolor='#FADADD',
          pos='2.6,1.6!']

    eff [label='Effluent\\nObs: C_eff',
         fillcolor='#FFF2C6',
         pos='0,-2.6!']

    out [label='Non-CRRT\\nElimination',
         fillcolor='#D4E4E8',
         pos='2.6,-1.6!']

    # ---- label nodes (hacky but intentional) ----
    node [shape=plaintext, style='', fontsize=10]

    # NOTE: lab_rate and lab_clcrrt are deliberately
    # left without positions/labels; we use edge labels
    # for those, which helps keep the vertical arrows
    # nice and straight in neato.
    lab_rate   [label='']
    lab_clcrrt [label='']

    lab_k12    [label='k12, k21',    pos='-1.4, 1.1!']
    lab_k15    [label='k15, k51',    pos='-1.4,-1.1!']
    lab_spost  [label='S_post',      pos=' 1.4, 1.1!']
    lab_clsys  [label='CL_systemic', pos=' 1.4,-1.1!']

    # ---- edges (labels handled via label nodes + edge labels) ----

    dose:s      -> central:n [label='Rate_in']

    central:nw  -> peripheral:se [dir=both]
    central:sw  -> elf:ne        [dir=both]
    central:ne  -> post:sw       [dir=both]

    central:s   -> eff:n         [label='CL_CRRT']
    central:se  -> out:nw
  }
  ")

viz_cef_struct_minimal

save_grviz(viz_cef_struct_minimal, "Figure2_mod_struct_flow_diagram")

# Cohort development validation plotting ####
# attempt at high-level orientation → 80/20 split → simulations flow diagram
# install if needed:
# install.packages(c("DiagrammeRsvg", "magick"))

library(DiagrammeR)
library(DiagrammeRsvg)
library(magick)

save_grviz <- function(grviz_obj,
                       filename,
                       out_dir  = "Manuscript/Figures",
                       dpi      = 800,
                       width_in = 5) {   # physical width for the TIFF only
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  svg_file  <- file.path(out_dir, paste0(filename, ".svg"))
  #tiff_file <- file.path(out_dir, paste0(filename, ".tiff"))
  
  ## --- SVG (vector, no dpi) ---
  svg_txt <- DiagrammeRsvg::export_svg(grviz_obj)
  writeLines(svg_txt, svg_file)
  
  # ## --- TIFF (raster, dpi matters) ---
  # img <- magick::image_read_svg(charToRaw(svg_txt))
  # 
  # # scale to requested physical width at given dpi
  # if (!is.null(width_in) && !is.null(dpi)) {
  #   width_px <- as.integer(width_in * dpi)
  #   img <- magick::image_scale(img, paste0(width_px))
  # }
  # 
  # magick::image_write(
  #   img,
  #   path       = tiff_file,
  #   format     = "tiff",
  #   compression = "lzw",
  #   density    = paste0(dpi, "x", dpi)  # dpi only relevant here
  # )
  # 
  message(
    "Saved figure to:\n  ",
    normalizePath(svg_file,  mustWork = FALSE)#, "\n  ",
    #normalizePath(tiff_file, mustWork = FALSE)
  )
  
  invisible(list(svg = svg_file#, 
                 #     tiff = tiff_file
  )
  )
}


flow_diagram <- grViz("
digraph cefepime_flow {

  graph [rankdir = TB, fontsize = 10]

  # global defaults
  node [
    shape    = box
    fontsize = 10
    fontname = 'Arial'
    penwidth = 1.1
  ]

  edge [
    arrowsize = 0.7
    penwidth  = 1
  ]

  # ---- color groups (soft viridis-inspired pastels) ----
  # data / cohorts:   #D4E4E8  (pale blue)
  # combined / split: #FFF2C6  (soft yellow)
  # model / sims:     #D3ECE4  (pale mint)

  # data sources + individual cohorts (blue, square corners)
  data_sources [
    label     = 'Cefepime PK data sources\nProspective, historical, salvaged',
    style     = 'filled',
    fillcolor = '#D4E4E8'
  ]

  prospective [
    label     = 'Prospective (richly sampled)\nELF, effluent, pre/post-filter',
    style     = 'filled',
    fillcolor = '#D4E4E8'
  ]

  historical [
    label     = 'Historical\nELF, pre-filter plasma',
    style     = 'filled',
    fillcolor = '#D4E4E8'
  ]

  salvaged [
    label     = 'Salvaged blood samples\nPre-filter plasma',
    style     = 'filled',
    fillcolor = '#D4E4E8'
  ]

  # combined dataset + 80/20 split (yellow, square corners)
  combined [
    label     = 'Combined PK dataset\nAll cohorts merged',
    style     = 'filled',
    fillcolor = '#FFF2C6'
  ]

  dev [
    label     = 'Model development cohort (80%)\nRandom 80% of subjects',
    style     = 'filled',
    fillcolor = '#FFF2C6'
  ]

  val [
    label     = 'Internal validation cohort (20%)\nRemaining 20% of subjects',
    style     = 'filled',
    fillcolor = '#FFF2C6'
  ]

  # model + simulations (mint, rounded corners)
  model [
    label     = 'Validated population PK model',
    style     = 'rounded,filled',
    fillcolor = '#D3ECE4'
  ]

  sims [
    label     = 'Monte Carlo simulations\nPre-filter plasma and ELF\nprobability of target attainment',
    style     = 'rounded,filled',
    fillcolor = '#D3ECE4'
  ]

  # layout hints
  { rank = same; prospective; historical; salvaged }
  { rank = same; dev; val }

  # edges
  data_sources -> prospective
  data_sources -> historical
  data_sources -> salvaged

  prospective -> combined
  historical  -> combined
  salvaged    -> combined

  combined -> dev
  combined -> val

  dev   -> model
  val   -> model

  model -> sims
}
")

save_grviz(flow_diagram, "Figure1_flow_diagram")

## Testing zone ####

library(dplyr)
library(tidyr)
library(ggplot2)

## ---- 1. Drill in: 1g EI, 60kg, MIC 8 ----

row_1g60_mic8 <- pta1$data$data %>%
  filter(label == "1g q12 EI + LD; 60kg",
         target == 8)

raw_1g60_mic8 <- row_1g60_mic8 %>%
  select(reg_num, label, target, success_ratio, prop_success, success, pdi) %>%
  unnest(cols = c(success, pdi)) %>%   # 1000 rows
  mutate(sim_id = row_number())

# PTA sanity check for this row
raw_1g60_mic8 %>%
  summarise(
    PTA_stored = first(prop_success),
    PTA_calc   = mean(success),
    n          = n()
  )

# PDI summary
raw_1g60_mic8 %>%
  summarise(
    min_pdi = min(pdi),
    q25     = quantile(pdi, 0.25),
    median  = median(pdi),
    q75     = quantile(pdi, 0.75),
    max_pdi = max(pdi)
  )

# Success vs pdi >= success_ratio
raw_1g60_mic8 %>%
  mutate(criterion = pdi >= success_ratio) %>%
  count(success, criterion)


## ---- 2. Compare 1g vs 2g EI, 60kg, MIC 8 ----

raw_60_mic8 <- pta1$data$data %>%
  filter(label %in% c("1g q12 EI + LD; 60kg",   # 1g EI
                      "2g q12 EI + LD; 60kg"),  # 2g EI
         target == 8) %>%
  select(label, success_ratio, prop_success, success, pdi) %>%
  unnest(cols = c(success, pdi)) %>%
  mutate(sim_id = row_number())

# Numeric comparison
raw_60_mic8 %>%
  group_by(label) %>%
  summarise(
    PTA_stored = first(prop_success),
    PTA_calc   = mean(success),
    median_pdi = median(pdi),
    q05_pdi    = quantile(pdi, 0.05),
    q95_pdi    = quantile(pdi, 0.95),
    .groups    = "drop"
  )

# Density plot
raw_60_mic8 %>%
  ggplot(aes(x = pdi, fill = label)) +
  geom_density(alpha = 0.4) +
  geom_vline(
    data = distinct(raw_60_mic8, label, success_ratio),
    aes(xintercept = success_ratio, color = label),
    linetype = 2
  ) +
  labs(
    title = "PDI distributions at MIC 8, 60kg (1g vs 2g EI)",
    x     = "PDI",
    y     = "Density"
  )


library(dplyr)
library(tidyr)
library(ggplot2)

## ---- 1. Drill in: 1g EI, 120kg, MIC 8 ----

row_1g120_mic8 <- pta1$data$data %>%
  filter(label == "1g q12 EI + LD; 120kg",
         target == 8)

raw_1g120_mic8 <- row_1g120_mic8 %>%
  select(reg_num, label, target, success_ratio, prop_success, success, pdi) %>%
  unnest(cols = c(success, pdi)) %>%   # 1000 rows
  mutate(sim_id = row_number())

# PTA sanity check for this row
raw_1g120_mic8 %>%
  summarise(
    PTA_stored = first(prop_success),
    PTA_calc   = mean(success),
    n          = n()
  )

# PDI summary
raw_1g120_mic8 %>%
  summarise(
    min_pdi = min(pdi),
    q25     = quantile(pdi, 0.25),
    median  = median(pdi),
    q75     = quantile(pdi, 0.75),
    max_pdi = max(pdi)
  )

# Success vs pdi >= success_ratio
raw_1g120_mic8 %>%
  mutate(criterion = pdi >= success_ratio) %>%
  count(success, criterion)


## ---- 2. Compare 1g vs 2g EI, 120kg, MIC 8 ----

raw_120_mic8 <- pta1$data$data %>%
  filter(label %in% c("1g q12 EI + LD; 120kg",   # 1g EI
                      "2g q12 EI + LD; 120kg"),  # 2g EI
         target == 8) %>%
  select(label, success_ratio, prop_success, success, pdi) %>%
  unnest(cols = c(success, pdi)) %>%
  mutate(sim_id = row_number())

# Numeric comparison
raw_120_mic8 %>%
  group_by(label) %>%
  summarise(
    PTA_stored = first(prop_success),
    PTA_calc   = mean(success),
    median_pdi = median(pdi),
    q05_pdi    = quantile(pdi, 0.05),
    q95_pdi    = quantile(pdi, 0.95),
    .groups    = "drop"
  )

# Density plot like the 60kg one
raw_120_mic8 %>%
  ggplot(aes(x = pdi, fill = label)) +
  geom_density(alpha = 0.4) +
  geom_vline(
    data = distinct(raw_120_mic8, label, success_ratio),
    aes(xintercept = success_ratio, color = label),
    linetype = 2
  ) +
  labs(
    title = "PDI distributions at MIC 8, 120kg (1g vs 2g EI)",
    x     = "PDI",
    y     = "Density"
  )

## VPC explorations for the holdout group ####

sim2 <- PM_sim$new(
  poppar = run10$final,
  data="Pmetrics/src/test2.csv",
  model=mod_cef6,
  seed = 12345,
  limits=c(0,1),
  nsim=1000,
  predInt = c(0.1)
)

# # External model validation fit
# 
# dat7 <- read.csv("Pmetrics/src/data.csv")
# 
# glimpse(dat7)
# 
# dat7 <- dat7 %>% select(-DIALYSATE,-UFR,-FILTER,-STUDY)
# 
# write.csv(dat7,"Pmetrics/src/ext.csv",row.names=F)
# 
# dat8 <- PM_data$new("Pmetrics/src/ext.csv",loq=c(0,0,0,0,0))
# 
# run13 <- mod_cef6$fit(data = dat8,
#                       cycles = 0,
#                       path = PATHS$runs,
#                       run = 13,
#                       prior = 10,
#                       overwrite=TRUE)
# 
# run13 <- PM_load(13, path = PATHS$runs)
# 

#### Average sampling times #### AV 
# ============================================================
# Confirm cefepime sampling times from manuscript data files
# ============================================================

library(dplyr)
library(readr)
library(readxl)
library(tidyr)
library(stringr)
library(purrr)

# ------------------------------------------------------------
# User settings: update these paths if needed
# ------------------------------------------------------------
src_dir   <- "../src"
xlsx_path <- file.path(src_dir, "FEP_comb_ELF_old_ML_combined_FDA.xlsx")
xlsx_sheet <- "FEP_comb_ELF_old"

train_path <- file.path(src_dir, "train2.csv")
test_path  <- file.path(src_dir, "test2.csv")
datml_path <- file.path(src_dir, "dat_ML.csv")

# Prospective IDs used in your workflow
prospective_ids <- c(2058, 2060, 2076, 2103, 2126, 9, 15, 20)

# OUTEQ labels after your remapping
outeq_labels <- c(
  "1" = "Pre-filter plasma",
  "2" = "Post-filter plasma",
  "3" = "Effluent concentration",
  "4" = "Effluent amount",
  "5" = "ELF"
)

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

clean_pm_csv <- function(path) {
  read_csv(path, show_col_types = FALSE, na = c(".", "NA", "")) %>%
    mutate(
      id         = suppressWarnings(as.numeric(id)),
      evid       = suppressWarnings(as.numeric(evid)),
      time       = suppressWarnings(as.numeric(time)),
      dose       = suppressWarnings(as.numeric(dose)),
      dur        = suppressWarnings(as.numeric(dur)),
      addl       = suppressWarnings(as.numeric(addl)),
      ii         = suppressWarnings(as.numeric(ii)),
      input      = suppressWarnings(as.numeric(input)),
      out        = suppressWarnings(as.numeric(out)),
      outeq      = suppressWarnings(as.numeric(outeq)),
      bag_reset  = suppressWarnings(as.numeric(bag_reset))
    ) %>%
    arrange(id, time, desc(evid))
}

add_cohort_labels <- function(df, prospective_ids = prospective_ids) {
  df %>%
    mutate(
      Cohort = case_when(
        id %in% prospective_ids ~ "Prospective (Richly Sampled)",
        id > 100 & id < 1700    ~ "Historical",
        id > 1700               ~ "Salvaged Blood Samples",
        TRUE                    ~ "Other / Unknown"
      )
    )
}

# Anchors each row to the most recent POSITIVE dose event
# Dose = 0 rows are not treated as actual doses
add_sampling_times <- function(df) {
  df %>%
    arrange(id, time, desc(evid)) %>%
    group_by(id) %>%
    group_modify(~{
      dat_id <- .x %>% arrange(time, desc(evid))
      
      dose_rows <- dat_id %>%
        filter(evid == 1, !is.na(dose), dose > 0) %>%
        transmute(
          dose_time = time,
          dose_end  = time + if_else(is.na(dur), 0, dur)
        ) %>%
        distinct(dose_time, .keep_all = TRUE) %>%
        arrange(dose_time)
      
      if (nrow(dose_rows) == 0) {
        dat_id$last_dose_time <- NA_real_
        dat_id$last_dose_end  <- NA_real_
        dat_id$tad_start      <- NA_real_
        dat_id$tad_end        <- NA_real_
        return(dat_id)
      }
      
      idx <- findInterval(dat_id$time, dose_rows$dose_time)
      
      last_dose_time <- rep(NA_real_, nrow(dat_id))
      last_dose_end  <- rep(NA_real_, nrow(dat_id))
      
      valid <- idx > 0
      last_dose_time[valid] <- dose_rows$dose_time[idx[valid]]
      last_dose_end[valid]  <- dose_rows$dose_end[idx[valid]]
      
      dat_id %>%
        mutate(
          last_dose_time = last_dose_time,
          last_dose_end  = last_dose_end,
          tad_start      = time - last_dose_time,
          tad_end        = time - last_dose_end
        )
    }) %>%
    ungroup()
}

# Keep only usable observation rows
keep_valid_obs <- function(df) {
  df %>%
    filter(
      evid == 0,
      !is.na(outeq),
      !is.na(out),
      out != -99,
      !is.na(last_dose_time)
    ) %>%
    mutate(
      Matrix = recode(as.character(outeq), !!!outeq_labels)
    )
}

round_timing_table <- function(df, digits = 2) {
  df %>%
    mutate(across(where(is.numeric), ~ round(.x, digits)))
}

summarise_sampling_times <- function(df, dataset_name = "dataset") {
  df %>%
    add_sampling_times() %>%
    keep_valid_obs() %>%
    group_by(Matrix) %>%
    summarise(
      n_obs             = n(),
      n_subjects        = n_distinct(id),
      median_tad_start  = median(tad_start, na.rm = TRUE),
      q1_tad_start      = quantile(tad_start, 0.25, na.rm = TRUE),
      q3_tad_start      = quantile(tad_start, 0.75, na.rm = TRUE),
      min_tad_start     = min(tad_start, na.rm = TRUE),
      max_tad_start     = max(tad_start, na.rm = TRUE),
      median_tad_end    = median(tad_end, na.rm = TRUE),
      q1_tad_end        = quantile(tad_end, 0.25, na.rm = TRUE),
      q3_tad_end        = quantile(tad_end, 0.75, na.rm = TRUE),
      min_tad_end       = min(tad_end, na.rm = TRUE),
      max_tad_end       = max(tad_end, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(dataset = dataset_name) %>%
    select(dataset, everything())
}

summarise_sampling_times_by_cohort <- function(df, dataset_name = "dataset") {
  df %>%
    add_cohort_labels() %>%
    add_sampling_times() %>%
    keep_valid_obs() %>%
    group_by(Cohort, Matrix) %>%
    summarise(
      n_obs             = n(),
      n_subjects        = n_distinct(id),
      median_tad_start  = median(tad_start, na.rm = TRUE),
      q1_tad_start      = quantile(tad_start, 0.25, na.rm = TRUE),
      q3_tad_start      = quantile(tad_start, 0.75, na.rm = TRUE),
      median_tad_end    = median(tad_end, na.rm = TRUE),
      q1_tad_end        = quantile(tad_end, 0.25, na.rm = TRUE),
      q3_tad_end        = quantile(tad_end, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(dataset = dataset_name) %>%
    select(dataset, Cohort, Matrix, everything())
}

make_manuscript_timing_table <- function(df, dataset_name = "dataset") {
  summarise_sampling_times(df, dataset_name) %>%
    transmute(
      dataset,
      Matrix,
      n = n_obs,
      `Subjects (n)` = n_subjects,
      `Median h after dose start (IQR)` = paste0(
        round(median_tad_start, 2), " (",
        round(q1_tad_start, 2), "–",
        round(q3_tad_start, 2), ")"
      ),
      `Median h after end of infusion (IQR)` = paste0(
        round(median_tad_end, 2), " (",
        round(q1_tad_end, 2), "–",
        round(q3_tad_end, 2), ")"
      )
    )
}

inspect_matrix_rows <- function(df, outeq_value, n = 25, tail = FALSE) {
  dat <- df %>%
    add_sampling_times() %>%
    keep_valid_obs() %>%
    filter(outeq == outeq_value) %>%
    select(id, time, Matrix, out, last_dose_time, last_dose_end, tad_start, tad_end)
  
  if (tail) {
    dat %>% arrange(desc(tad_start)) %>% head(n)
  } else {
    dat %>% arrange(tad_start) %>% head(n)
  }
}

dose_duration_summary <- function(df) {
  df %>%
    filter(evid == 1, !is.na(dose), dose > 0) %>%
    count(dur, sort = TRUE)
}

count_valid_obs <- function(df) {
  df %>%
    add_sampling_times() %>%
    keep_valid_obs() %>%
    count(Matrix)
}

# ------------------------------------------------------------
# Optional: rebuild pooled dat4 directly from Excel source
# using the same logic as your workflow
# ------------------------------------------------------------
rebuild_dat4_from_excel <- function(xlsx_path, xlsx_sheet) {
  dat <- read_excel(xlsx_path, sheet = xlsx_sheet)
  
  row_special <- dat %>%
    mutate(
      ID   = suppressWarnings(as.numeric(ID)),
      DVID = suppressWarnings(as.numeric(DVID))
    ) %>%
    filter(ID == 2103, DVID == 2) %>%
    mutate(EXCLUDE_DV = 0)
  
  dat4 <- dat %>%
    mutate(
      ID   = suppressWarnings(as.numeric(ID)),
      DVID = suppressWarnings(as.numeric(DVID)),
      RRT  = suppressWarnings(as.numeric(RRT)),
      EVID = suppressWarnings(as.numeric(EVID))
    ) %>%
    filter(RRT == 1, EXCLUDE_DV != 1) %>%
    rename(
      DOSE  = AMT,
      DUR   = TINF,
      INPUT = ADM,
      OUTEQ = DVID,
      OUT   = DV,
      LOQ   = LIMIT
    ) %>%
    mutate(
      OUTEQ = suppressWarnings(as.numeric(OUTEQ)),
      LOQ   = suppressWarnings(as.numeric(LOQ)),
      OUT   = suppressWarnings(as.numeric(OUT)),
      TIME  = suppressWarnings(as.numeric(TIME)),
      DOSE  = suppressWarnings(as.numeric(DOSE)),
      DUR   = suppressWarnings(as.numeric(DUR)),
      BAG_RESET = suppressWarnings(as.numeric(BAG_RESET))
    ) %>%
    mutate(
      OUTEQ = case_when(
        OUTEQ == 1 ~ 1,  # pre-filter
        OUTEQ == 3 ~ 2,  # post-filter
        OUTEQ == 4 ~ 3,  # effluent concentration
        OUTEQ == 5 ~ 4,  # effluent amount
        OUTEQ == 2 ~ 5,  # ELF
        TRUE ~ NA_real_
      ),
      OUT = case_when(
        EVID == 0 & is.na(OUT) ~ -99,
        TRUE ~ OUT
      ),
      BAG_RESET = case_when(
        INPUT == 2 ~ 0,
        TRUE ~ BAG_RESET
      )
    ) %>%
    select(ID, EVID, TIME, DOSE, DUR, ADDL, II, INPUT, OUT, OUTEQ, everything()) %>%
    rename_with(tolower)
  
  row_special_std <- row_special %>%
    rename(
      DOSE  = AMT,
      DUR   = TINF,
      INPUT = ADM,
      OUTEQ = DVID,
      OUT   = DV
    ) %>%
    mutate(
      ID        = suppressWarnings(as.numeric(ID)),
      EVID      = suppressWarnings(as.numeric(EVID)),
      TIME      = suppressWarnings(as.numeric(TIME)),
      DOSE      = suppressWarnings(as.numeric(DOSE)),
      DUR       = suppressWarnings(as.numeric(DUR)),
      INPUT     = suppressWarnings(as.numeric(INPUT)),
      OUT       = suppressWarnings(as.numeric(OUT)),
      OUTEQ     = 5,
      BAG_RESET = suppressWarnings(as.numeric(BAG_RESET))
    ) %>%
    rename_with(tolower)
  
  dat4_full <- dat4 %>%
    anti_join(
      row_special_std %>% select(id, time, outeq),
      by = c("id", "time", "outeq")
    ) %>%
    bind_rows(row_special_std) %>%
    arrange(id, time, desc(evid))
  
  dat4_full
}

# ------------------------------------------------------------
# Load the actual files used in the manuscript workflow
# ------------------------------------------------------------
train2 <- clean_pm_csv(train_path)
test2  <- clean_pm_csv(test_path)
dat_ml <- clean_pm_csv(datml_path)

final_comb <- bind_rows(train2, test2) %>%
  arrange(id, time, desc(evid))

# Optional rebuilt pooled source directly from Excel
dat4_full <- rebuild_dat4_from_excel(xlsx_path, xlsx_sheet)

# ------------------------------------------------------------
# 1) Main timing summaries
# ------------------------------------------------------------
final_timing <- summarise_sampling_times(final_comb, "train2 + test2") %>%
  round_timing_table()

datml_timing <- summarise_sampling_times(dat_ml, "dat_ML") %>%
  round_timing_table()

dat4full_timing <- summarise_sampling_times(dat4_full, "rebuilt dat4_full") %>%
  round_timing_table()

cat("\n================ FINAL POOLED MANUSCRIPT DATASET ================\n")
print(final_timing)

cat("\n================ BASE CRRT STRUCTURAL DATASET (dat_ML) ================\n")
print(datml_timing)

cat("\n================ REBUILT DAT4 FROM EXCEL SOURCE ================\n")
print(dat4full_timing)

# ------------------------------------------------------------
# 2) By-cohort timing summaries
# ------------------------------------------------------------
final_timing_by_cohort <- summarise_sampling_times_by_cohort(final_comb, "train2 + test2") %>%
  round_timing_table()

datml_timing_by_cohort <- summarise_sampling_times_by_cohort(dat_ml, "dat_ML") %>%
  round_timing_table()

cat("\n================ FINAL DATASET BY COHORT ================\n")
print(final_timing_by_cohort)

cat("\n================ dat_ML BY COHORT ================\n")
print(datml_timing_by_cohort)

# ------------------------------------------------------------
# 3) Manuscript-friendly timing table
# ------------------------------------------------------------
manuscript_timing <- make_manuscript_timing_table(final_comb, "train2 + test2")

cat("\n================ MANUSCRIPT-FRIENDLY TIMING TABLE ================\n")
print(manuscript_timing)

# ------------------------------------------------------------
# 4) Confirm valid observation counts contributing to summaries
# ------------------------------------------------------------
cat("\n================ VALID OBSERVATION COUNTS: FINAL DATASET ================\n")
print(count_valid_obs(final_comb))

cat("\n================ VALID OBSERVATION COUNTS: dat_ML ================\n")
print(count_valid_obs(dat_ml))

# ------------------------------------------------------------
# 5) Confirm infusion durations
# ------------------------------------------------------------
cat("\n================ DOSE DURATION SUMMARY: FINAL DATASET ================\n")
print(dose_duration_summary(final_comb))

cat("\n================ DOSE DURATION SUMMARY: dat_ML ================\n")
print(dose_duration_summary(dat_ml))

# ------------------------------------------------------------
# 6) Inspect raw rows if a matrix looks odd
#    Change outeq_value as needed:
#      1 = pre-filter
#      2 = post-filter
#      3 = effluent concentration
#      4 = effluent amount
#      5 = ELF
# ------------------------------------------------------------
cat("\n================ FIRST 25 ELF ROWS SORTED BY TAD_START ================\n")
print(inspect_matrix_rows(final_comb, outeq_value = 5, n = 25, tail = FALSE))

cat("\n================ TOP 25 PRE-FILTER ROWS WITH LARGEST TAD_START ================\n")
print(inspect_matrix_rows(final_comb, outeq_value = 1, n = 25, tail = TRUE))

# ------------------------------------------------------------
# 7) Optional: export results to CSV for review
# ------------------------------------------------------------
write_csv(final_timing, file.path(src_dir, "timing_summary_final_combined.csv"))
write_csv(datml_timing, file.path(src_dir, "timing_summary_dat_ML.csv"))
write_csv(dat4full_timing, file.path(src_dir, "timing_summary_rebuilt_dat4_full.csv"))
write_csv(final_timing_by_cohort, file.path(src_dir, "timing_summary_final_by_cohort.csv"))
write_csv(manuscript_timing, file.path(src_dir, "timing_table_manuscript_ready.csv"))

cat("\nCSV outputs written to Pmetrics/src/\n")

# ------------------------------------------------------------
# 8) Optional: quick cross-check between final pooled files and rebuilt dat4
# ------------------------------------------------------------
compare_final_vs_rebuilt <- final_timing %>%
  select(Matrix, final_n = n_obs,
         final_med_start = median_tad_start,
         final_med_end = median_tad_end) %>%
  full_join(
    dat4full_timing %>%
      select(Matrix, rebuilt_n = n_obs,
             rebuilt_med_start = median_tad_start,
             rebuilt_med_end = median_tad_end),
    by = "Matrix"
  )

cat("\n================ FINAL FILES VS REBUILT DAT4 CROSS-CHECK ================\n")
print(compare_final_vs_rebuilt)

# ------------------------------------------------------------
# End
# ------------------------------------------------------------
