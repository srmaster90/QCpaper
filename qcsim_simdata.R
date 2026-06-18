# S Master, Children's Hospital of Philadelphia
# masters@chop.edu
# last edit: 2/6/26
#
# generate QC simulation data

library(qcsim)
cores <- parallel::detectCores() - 1
# cores <- 6

my.seed <- 153
######  Figures for paper

iters <- 101
maxs <- 5
maxs_x <- 7 # (for figures with 5s curves)

# Simulate 3S rule (bias) at 1/10/50/100-plex, 1 level of QC
pres.bias.3s.1 <- qcsim(maxs=maxs, xsize=iters, plex=1, cores=cores, seed=my.seed)
pres.bias.3s.10 <- qcsim(maxs=maxs, xsize=iters, plex=10, cores=cores, seed=my.seed)
pres.bias.3s.50 <- qcsim(maxs=maxs, xsize=iters, plex=50, cores=cores, seed=my.seed)
pres.bias.3s.100 <- qcsim(maxs=maxs, xsize=iters, plex=100, cores=cores, seed=my.seed)

# Simulate 3S rule (imprecision) at 1/10/50/100-plex, 1 level of QC
pres.prec.3s.1 <- qcsim(maxs=maxs, xsize=iters, plex=1, precision=T, cores=cores, seed=my.seed)
pres.prec.3s.10 <- qcsim(maxs=maxs, xsize=iters, plex=10, precision=T, cores=cores, seed=my.seed)
pres.prec.3s.50 <- qcsim(maxs=maxs, xsize=iters, plex=50, precision=T, cores=cores, seed=my.seed)
pres.prec.3s.100 <- qcsim(maxs=maxs, xsize=iters, plex=100, precision=T, cores=cores, seed=my.seed)

# Simulate 2S rule (bias) at 1/10/50/100-plex, 1 level of QC
pres.bias.2s.1 <- qcsim(maxs=maxs, xsize=iters, plex=1, cores=cores, seed=my.seed, sig=2)
pres.bias.2s.10 <- qcsim(maxs=maxs, xsize=iters, plex=10, cores=cores, seed=my.seed, sig=2)
pres.bias.2s.50 <- qcsim(maxs=maxs, xsize=iters, plex=50, cores=cores, seed=my.seed, sig=2)
pres.bias.2s.100 <- qcsim(maxs=maxs, xsize=iters, plex=100, cores=cores, seed=my.seed, sig=2)

# Simulate 2S rule (imprecision) at 1/10/50/100-plex, 1 level of QC
pres.prec.2s.1 <- qcsim(maxs=maxs, xsize=iters, plex=1, precision=T, cores=cores, seed=my.seed, sig=2)
pres.prec.2s.10 <- qcsim(maxs=maxs, xsize=iters, plex=10, precision=T, cores=cores, seed=my.seed, sig=2)
pres.prec.2s.50 <- qcsim(maxs=maxs, xsize=iters, plex=50, precision=T, cores=cores, seed=my.seed, sig=2)
pres.prec.2s.100 <- qcsim(maxs=maxs, xsize=iters, plex=100, precision=T, cores=cores, seed=my.seed, sig=2)

save.image(file="Results/qcdata_interim.RData")

# Simulate 3S rule (bias) at 1/10/50/100-plex, 2 levels of QC
pres.bias.3s2.1 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=1, cores=cores, seed=my.seed)
pres.bias.3s2.10 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=10, cores=cores, seed=my.seed)
pres.bias.3s2.50 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=50, cores=cores, seed=my.seed)
pres.bias.3s2.100 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=100, cores=cores, seed=my.seed)

# Simulate 3S rule (imprecision) at 1/10/50/100-plex, 2 levels of QC
pres.prec.3s2.1 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=1, precision=T, cores=cores, seed=my.seed)
pres.prec.3s2.10 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=10, precision=T, cores=cores, seed=my.seed)
pres.prec.3s2.50 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=50, precision=T, cores=cores, seed=my.seed)
pres.prec.3s2.100 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=100, precision=T, cores=cores, seed=my.seed)

# see above.  NB: "_x" versions have extended range for comparioson with 5s
pres.bias.3s2_x.1 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=1, cores=cores, seed=my.seed)
pres.bias.3s2_x.10 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=10, cores=cores, seed=my.seed)
pres.bias.3s2_x.50 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=50, cores=cores, seed=my.seed)
pres.bias.3s2_x.100 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=100, cores=cores, seed=my.seed)

# see above.  NB: "_x" versions have extended range for comparioson with 5s
pres.prec.3s2_x.1 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=1, precision=T, cores=cores, seed=my.seed)
pres.prec.3s2_x.10 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=10, precision=T, cores=cores, seed=my.seed)
pres.prec.3s2_x.50 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=50, precision=T, cores=cores, seed=my.seed)
pres.prec.3s2_x.100 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=100, precision=T, cores=cores, seed=my.seed)

# ...now back to regular range.  # Simulate 2S rule (bias) at 1/10/50/100-plex, 2 levels of QC
pres.bias.2s2.1 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=1, cores=cores, seed=my.seed, sig=2)
pres.bias.2s2.10 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=10, cores=cores, seed=my.seed, sig=2)
pres.bias.2s2.50 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=50, cores=cores, seed=my.seed, sig=2)
pres.bias.2s2.100 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=100, cores=cores, seed=my.seed, sig=2)

# Simulate 2S rule (imprecision) at 1/10/50/100-plex, 2 levels of QC
pres.prec.2s2.1 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=1, precision=T, cores=cores, seed=my.seed, sig=2)
pres.prec.2s2.10 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=10, precision=T, cores=cores, seed=my.seed, sig=2)
pres.prec.2s2.50 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=50, precision=T, cores=cores, seed=my.seed, sig=2)
pres.prec.2s2.100 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=100, precision=T, cores=cores, seed=my.seed, sig=2)

save.image(file="Results/qcdata_interim.RData")


# Simulate Westgard rules (1-3S, 2-2S, R-4S, 4-1S, 10X) (bias) at 1/10/50/100-plex, 2 levels of QC
pres.bias.wg2.1 <- qcsim(maxs=maxs, xsize=iters, plex=1,
                         levels=2,
                         rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                         cores=cores, seed=my.seed)
pres.bias.wg2.10 <- qcsim(maxs=maxs, xsize=iters, plex=10,
                          levels=2,
                          rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                          cores=cores, seed=my.seed)
pres.bias.wg2.50 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                          levels=2,
                          rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                          cores=cores, seed=my.seed)
pres.bias.wg2.100 <- qcsim(maxs=maxs, xsize=iters, plex=100,
                           levels=2,
                           rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                           cores=cores, seed=my.seed)

# Simulate Westgard rules (1-3S, 2-2S, R-4S, 4-1S, 10X) (imprecision) at 1/10/50/100-plex, 2 levels of QC
pres.prec.wg2.1 <- qcsim(maxs=maxs, xsize=iters, plex=1,
                         levels=2,
                         rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                         precision=T, cores=cores, seed=my.seed)
pres.prec.wg2.10 <- qcsim(maxs=maxs, xsize=iters, plex=10,
                          levels=2,
                          rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                          precision=T, cores=cores, seed=my.seed)
pres.prec.wg2.50 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                          levels=2,
                          rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                          precision=T, cores=cores, seed=my.seed)
pres.prec.wg2.100 <- qcsim(maxs=maxs, xsize=iters, plex=100,
                           levels=2,
                           rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                           precision=T, cores=cores, seed=my.seed)

# Simulate Westgard rules (1-3S, 2-2S, R-4S, 4-1S, 10X) (bias) at 1/10/50/100-plex, 1 level of QC
pres.bias.wg1.1 <- qcsim(maxs=maxs, xsize=iters, plex=1,
                         levels=1,
                         rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                         cores=cores, seed=my.seed)
pres.bias.wg1.10 <- qcsim(maxs=maxs, xsize=iters, plex=10,
                          levels=1,
                          rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                          cores=cores, seed=my.seed)
pres.bias.wg1.50 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                          levels=1,
                          rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                          cores=cores, seed=my.seed)
pres.bias.wg1.100 <- qcsim(maxs=maxs, xsize=iters, plex=100,
                           levels=1,
                           rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                           cores=cores, seed=my.seed)

# Simulate Westgard rules (1-3S, 2-2S, R-4S, 4-1S, 10X) (imprecision) at 1/10/50/100-plex, 1 level of QC
pres.prec.wg1.1 <- qcsim(maxs=maxs, xsize=iters, plex=1,
                         levels=1,
                         rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                         precision=T, cores=cores, seed=my.seed)
pres.prec.wg1.10 <- qcsim(maxs=maxs, xsize=iters, plex=10,
                          levels=1,
                          rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                          precision=T, cores=cores, seed=my.seed)
pres.prec.wg1.50 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                          levels=1,
                          rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                          precision=T, cores=cores, seed=my.seed)
pres.prec.wg1.100 <- qcsim(maxs=maxs, xsize=iters, plex=100,
                           levels=1,
                           rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                           precision=T, cores=cores, seed=my.seed)

save.image(file="Results/qcdata_interim.RData")

# Simulate 4S rule (bias) at 1/10/50/100-plex, 1 level of QC
pres.bias.4s.1 <- qcsim(maxs=maxs_x, xsize=iters, plex=1, cores=cores, seed=my.seed, sig=4)
pres.bias.4s.10 <- qcsim(maxs=maxs_x, xsize=iters, plex=10, cores=cores, seed=my.seed, sig=4)
pres.bias.4s.50 <- qcsim(maxs=maxs_x, xsize=iters, plex=50, cores=cores, seed=my.seed, sig=4)
pres.bias.4s.100 <- qcsim(maxs=maxs_x, xsize=iters, plex=100, cores=cores, seed=my.seed, sig=4)

# Simulate 4S rule (imprecision) at 1/10/50/100-plex, 1 level of QC
pres.prec.4s.1 <- qcsim(maxs=maxs_x, xsize=iters, plex=1, precision=T, cores=cores, seed=my.seed, sig=4)
pres.prec.4s.10 <- qcsim(maxs=maxs_x, xsize=iters, plex=10, precision=T, cores=cores, seed=my.seed, sig=4)
pres.prec.4s.50 <- qcsim(maxs=maxs_x, xsize=iters, plex=50, precision=T, cores=cores, seed=my.seed, sig=4)
pres.prec.4s.100 <- qcsim(maxs=maxs_x, xsize=iters, plex=100, precision=T, cores=cores, seed=my.seed, sig=4)

# Simulate 4S rule (bias) at 1/10/50/100-plex, 2 levels of QC
pres.bias.4s2.1 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=1, cores=cores, seed=my.seed, sig=4)
pres.bias.4s2.10 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=10, cores=cores, seed=my.seed, sig=4)
pres.bias.4s2.50 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=50, cores=cores, seed=my.seed, sig=4)
pres.bias.4s2.100 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=100, cores=cores, seed=my.seed, sig=4)

# Simulate 4S rule (imprecision) at 1/10/50/100-plex, 2 levels of QC
pres.prec.4s2.1 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=1, precision=T, cores=cores, seed=my.seed, sig=4)
pres.prec.4s2.10 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=10, precision=T, cores=cores, seed=my.seed, sig=4)
pres.prec.4s2.50 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=50, precision=T, cores=cores, seed=my.seed, sig=4)
pres.prec.4s2.100 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=100, precision=T, cores=cores, seed=my.seed, sig=4)

# Simulate 5S rule (bias) at 1/10/50/100-plex, 1 level of QC
pres.bias.5s.1 <- qcsim(maxs=maxs_x, xsize=iters, plex=1, cores=cores, seed=my.seed, sig=5)
pres.bias.5s.10 <- qcsim(maxs=maxs_x, xsize=iters, plex=10, cores=cores, seed=my.seed, sig=5)
pres.bias.5s.50 <- qcsim(maxs=maxs_x, xsize=iters, plex=50, cores=cores, seed=my.seed, sig=5)
pres.bias.5s.100 <- qcsim(maxs=maxs_x, xsize=iters, plex=100, cores=cores, seed=my.seed, sig=5)

# Simulate 5S rule (imprecision) at 1/10/50/100-plex, 1 level of QC
pres.prec.5s.1 <- qcsim(maxs=maxs_x, xsize=iters, plex=1, precision=T, cores=cores, seed=my.seed, sig=5)
pres.prec.5s.10 <- qcsim(maxs=maxs_x, xsize=iters, plex=10, precision=T, cores=cores, seed=my.seed, sig=5)
pres.prec.5s.50 <- qcsim(maxs=maxs_x, xsize=iters, plex=50, precision=T, cores=cores, seed=my.seed, sig=5)
pres.prec.5s.100 <- qcsim(maxs=maxs_x, xsize=iters, plex=100, precision=T, cores=cores, seed=my.seed, sig=5)

# Simulate 5S rule (bias) at 1/10/50/100-plex, 2 levels of QC
pres.bias.5s2.1 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=1, cores=cores, seed=my.seed, sig=5)
pres.bias.5s2.10 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=10, cores=cores, seed=my.seed, sig=5)
pres.bias.5s2.50 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=50, cores=cores, seed=my.seed, sig=5)
pres.bias.5s2.100 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=100, cores=cores, seed=my.seed, sig=5)

# Simulate 5S rule (imprecision) at 1/10/50/100-plex, 2 levels of QC
pres.prec.5s2.1 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=1, precision=T, cores=cores, seed=my.seed, sig=5)
pres.prec.5s2.10 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=10, precision=T, cores=cores, seed=my.seed, sig=5)
pres.prec.5s2.50 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=50, precision=T, cores=cores, seed=my.seed, sig=5)
pres.prec.5s2.100 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=100, precision=T, cores=cores, seed=my.seed, sig=5)

save.image(file="Results/qcdata_interim.RData")

# T and F test data, 50-plex, 2 levels of QC (see individual function calls for cutoff values and group lenghts)
pres.bias.tf2.50_25ref_3rep_1en3 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                         levels=2, reps=3,
                                         ruleP_2S=F, rule1_sigS=F,
                                         ruleT=T, cutoff.T=1e-3, ruleF=T, cutoff.F=1e-3,
                                         FTrefLength=25, FTtestLength=3,
                                         precision=F, cores=cores, seed=my.seed)

pres.bias.tf2.50_25ref_3rep_5en4 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                          levels=2, reps=3,
                                          ruleP_2S=F, rule1_sigS=F,
                                          ruleT=T, cutoff.T=5e-4, ruleF=T, cutoff.F=5e-4,
                                          FTrefLength=25, FTtestLength=3,
                                          precision=F, cores=cores, seed=my.seed)

pres.bias.tf2.50_25ref_3rep_5en5 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                           levels=2, reps=3,
                                           ruleP_2S=F, rule1_sigS=F,
                                           ruleT=T, cutoff.T=5e-5, ruleF=T, cutoff.F=5e-5,
                                           FTrefLength=25, FTtestLength=3,
                                           precision=F, cores=cores, seed=my.seed)

pres.bias.tf2.50_25ref_5rep_1en3 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                         levels=2, reps=5,
                                         ruleP_2S=F, rule1_sigS=F,
                                         ruleT=T, cutoff.T=1e-3, ruleF=T, cutoff.F=1e-3,
                                         FTrefLength=25, FTtestLength=5,
                                         precision=F, cores=cores, seed=my.seed)

pres.bias.tf2.50_25ref_5rep_5en4 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                          levels=2, reps=5,
                                          ruleP_2S=F, rule1_sigS=F,
                                          ruleT=T, cutoff.T=5e-4, ruleF=T, cutoff.F=5e-4,
                                          FTrefLength=25, FTtestLength=5,
                                          precision=F, cores=cores, seed=my.seed)

pres.bias.tf2.50_25ref_5rep_5en5 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                           levels=2, reps=5,
                                           ruleP_2S=F, rule1_sigS=F,
                                           ruleT=T, cutoff.T=5e-5, ruleF=T, cutoff.F=5e-5,
                                           FTrefLength=25, FTtestLength=5,
                                           precision=F, cores=cores, seed=my.seed)

pres.bias.tf2.50_25ref_10rep_1en3 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                          levels=2, reps=10,
                                          ruleP_2S=F, rule1_sigS=F,
                                          ruleT=T, cutoff.T=1e-3, ruleF=T, cutoff.F=1e-3,
                                          FTrefLength=25, FTtestLength=10,
                                          precision=F, cores=cores, seed=my.seed)

pres.bias.tf2.50_25ref_10rep_5en4 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                           levels=2, reps=10,
                                           ruleP_2S=F, rule1_sigS=F,
                                           ruleT=T, cutoff.T=5e-4, ruleF=T, cutoff.F=5e-4,
                                           FTrefLength=25, FTtestLength=10,
                                           precision=F, cores=cores, seed=my.seed)

pres.bias.tf2.50_25ref_10rep_5en5 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                            levels=2, reps=10,
                                            ruleP_2S=F, rule1_sigS=F,
                                            ruleT=T, cutoff.T=5e-5, ruleF=T, cutoff.F=5e-5,
                                            FTrefLength=25, FTtestLength=10,
                                            precision=F, cores=cores, seed=my.seed)

pres.prec.tf2.50_25ref_3rep_1en3 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                          levels=2, reps=3,
                                          ruleP_2S=F, rule1_sigS=F,
                                          ruleT=T, cutoff.T=1e-3, ruleF=T, cutoff.F=1e-3,
                                          FTrefLength=25, FTtestLength=3,
                                          precision=T, cores=cores, seed=my.seed)

pres.prec.tf2.50_25ref_3rep_5en4 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                          levels=2, reps=3,
                                          ruleP_2S=F, rule1_sigS=F,
                                          ruleT=T, cutoff.T=5e-4, ruleF=T, cutoff.F=5e-4,
                                          FTrefLength=25, FTtestLength=3,
                                          precision=T, cores=cores, seed=my.seed)

pres.prec.tf2.50_25ref_3rep_5en5 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                          levels=2, reps=3,
                                          ruleP_2S=F, rule1_sigS=F,
                                          ruleT=T, cutoff.T=5e-5, ruleF=T, cutoff.F=5e-5,
                                          FTrefLength=25, FTtestLength=3,
                                          precision=T, cores=cores, seed=my.seed)

pres.prec.tf2.50_25ref_5rep_1en3 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                          levels=2, reps=5,
                                          ruleP_2S=F, rule1_sigS=F,
                                          ruleT=T, cutoff.T=1e-3, ruleF=T, cutoff.F=1e-3,
                                          FTrefLength=25, FTtestLength=5,
                                          precision=T, cores=cores, seed=my.seed)

pres.prec.tf2.50_25ref_5rep_5en4 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                          levels=2, reps=5,
                                          ruleP_2S=F, rule1_sigS=F,
                                          ruleT=T, cutoff.T=5e-4, ruleF=T, cutoff.F=5e-4,
                                          FTrefLength=25, FTtestLength=5,
                                          precision=T, cores=cores, seed=my.seed)

pres.prec.tf2.50_25ref_5rep_5en5 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                          levels=2, reps=5,
                                          ruleP_2S=F, rule1_sigS=F,
                                          ruleT=T, cutoff.T=5e-5, ruleF=T, cutoff.F=5e-5,
                                          FTrefLength=25, FTtestLength=5,
                                          precision=T, cores=cores, seed=my.seed)

pres.prec.tf2.50_25ref_10rep_1en3 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                           levels=2, reps=10,
                                           ruleP_2S=F, rule1_sigS=F,
                                           ruleT=T, cutoff.T=1e-3, ruleF=T, cutoff.F=1e-3,
                                           FTrefLength=25, FTtestLength=10,
                                           precision=T, cores=cores, seed=my.seed)

pres.prec.tf2.50_25ref_10rep_5en4 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                           levels=2, reps=10,
                                           ruleP_2S=F, rule1_sigS=F,
                                           ruleT=T, cutoff.T=5e-4, ruleF=T, cutoff.F=5e-4,
                                           FTrefLength=25, FTtestLength=10,
                                           precision=T, cores=cores, seed=my.seed)

pres.prec.tf2.50_25ref_10rep_5en5 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                           levels=2, reps=10,
                                           ruleP_2S=F, rule1_sigS=F,
                                           ruleT=T, cutoff.T=5e-5, ruleF=T, cutoff.F=5e-5,
                                           FTrefLength=25, FTtestLength=10,
                                           precision=T, cores=cores, seed=my.seed)

# same as above, with T or F test alone

pres.bias.t2.50_25ref_10rep_5en5 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                            levels=2, reps=10,
                                            ruleP_2S=F, rule1_sigS=F,
                                            ruleT=T, cutoff.T=5e-5, ruleF=F,
                                            FTrefLength=25, FTtestLength=10,
                                            precision=F, cores=cores, seed=my.seed)

pres.bias.f2.50_25ref_10rep_5en5 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                            levels=2, reps=10,
                                            ruleP_2S=F, rule1_sigS=F,
                                            ruleT=F, ruleF=T, cutoff.F=5e-5,
                                            FTrefLength=25, FTtestLength=10,
                                            precision=F, cores=cores, seed=my.seed)

pres.prec.t2.50_25ref_10rep_5en5 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                            levels=2, reps=10,
                                            ruleP_2S=F, rule1_sigS=F,
                                            ruleT=T, cutoff.T=5e-5, ruleF=F,
                                            FTrefLength=25, FTtestLength=10,
                                            precision=T, cores=cores, seed=my.seed)

pres.prec.f2.50_25ref_10rep_5en5 <- qcsim(maxs=maxs, xsize=iters, plex=50,
                                            levels=2, reps=10,
                                            ruleP_2S=F, rule1_sigS=F,
                                            ruleT=F, ruleF=T, cutoff.F=5e-5,
                                            FTrefLength=25, FTtestLength=10,
                                            precision=T, cores=cores, seed=my.seed)

save.image(file="Results/qcdata_interim.RData")

#1000-plex (1-3S, 1-2S at 1 and 2 levels of QC)

pres.bias.3s.1000 <- qcsim(maxs=maxs, xsize=iters, plex=1000, cores=cores, seed=my.seed)
pres.prec.3s.1000 <- qcsim( maxs=maxs, xsize=iters, plex=1000, precision=T, cores=cores, seed=my.seed)
pres.bias.2s.1000 <- qcsim(maxs=maxs, xsize=iters, plex=1000, cores=cores, seed=my.seed, sig=2)
pres.prec.2s.1000 <- qcsim( maxs=maxs, xsize=iters, plex=1000, precision=T, cores=cores, seed=my.seed, sig=2)
pres.bias.3s2.1000 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=1000, cores=cores, seed=my.seed)
pres.prec.3s2.1000 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=1000, precision=T, cores=cores, seed=my.seed)
pres.bias.3s2_x.1000 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=1000, cores=cores, seed=my.seed)
pres.prec.3s2_x.1000 <- qcsim(maxs=maxs_x, xsize=iters, levels=2, plex=1000, precision=T, cores=cores, seed=my.seed)
pres.bias.2s2.1000 <- qcsim(maxs=maxs, xsize=iters, levels=2, plex=1000, cores=cores, seed=my.seed, sig=2)
pres.prec.2s2.1000 <- qcsim( maxs=maxs, xsize=iters, levels=2, plex=1000, precision=T, cores=cores, seed=my.seed, sig=2)

# 1000-plex, ctd., full Westgard rules (see above)
pres.bias.wg2.1000 <- qcsim(maxs=maxs, xsize=iters, plex=1000,
                            levels=2,
                            rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                            cores=cores, seed=my.seed)
pres.prec.wg2.1000 <- qcsim(maxs=maxs, xsize=iters, plex=1000,
                            levels=2,
                            rule2_2S=T, ruleR_4S=T, rule4_1S=T, rule10X=T,
                            precision=T, cores=cores, seed=my.seed)

# 1000-plex, 1-4S and 1-5S, 2 levels of QC
pres.bias.4s2.1000 <- qcsim( maxs=maxs_x, xsize=iters, levels=2, plex=1000, precision=F, cores=cores, seed=my.seed, sig=4)
pres.bias.5s2.1000 <- qcsim( maxs=maxs_x, xsize=iters, levels=2, plex=1000, precision=F, cores=cores, seed=my.seed, sig=5)
pres.prec.4s2.1000 <- qcsim( maxs=maxs_x, xsize=iters, levels=2, plex=1000, precision=T, cores=cores, seed=my.seed, sig=4)
pres.prec.5s2.1000 <- qcsim( maxs=maxs_x, xsize=iters, levels=2, plex=1000, precision=T, cores=cores, seed=my.seed, sig=5)

# 1000-plex, T and F tests, 2 levels of QC

# (note that T and F tests with large number of analytes lead to high memory useage,
#  so it may be necessary to lower the number of cores for highly multiplex situations)

pres.bias.tf2.1000_25ref_10rep_2.7en6 <- qcsim(maxs=maxs, xsize=iters, plex=1000,
                                           levels=2, reps=10,
                                           ruleP_2S=F, rule1_sigS=F,
                                           ruleT=T, cutoff.T=2.7e-6, ruleF=T, cutoff.F=2.7e-6,
                                           FTrefLength=25, FTtestLength=10,
                                           precision=F, cores=2, seed=my.seed)

pres.prec.tf2.1000_25ref_10rep_2.7en6 <- qcsim(maxs=maxs, xsize=iters, plex=1000,
                                            levels=2, reps=10,
                                            ruleP_2S=F, rule1_sigS=F,
                                            ruleT=T, cutoff.T=2.7e-6, ruleF=T, cutoff.F=2.7e-6,
                                            FTrefLength=25, FTtestLength=10,
                                            precision=T, cores=2, seed=my.seed)

save.image(file="Results/qcdata.RData")

# sample figure generation
#
# pdf("Figure_1.pdf")
#
# par(mfrow=c(1,2))
#
# plot(pres.bias.3s.1)
# lines(pres.bias.3s.10, lty=2)
# lines(pres.bias.3s.50, lty=3)
# lines(pres.bias.3s.100, lty=4)
# legend("bottomright", legend=c("1-plex", "10-plex", "50-plex", "100-plex"),
#        lty=1:4)
#
# plot(pres.prec.3s.1)
# lines(pres.prec.3s.10, lty=2)
# lines(pres.prec.3s.50, lty=3)
# lines(pres.prec.3s.100, lty=4)
# legend("bottomright", legend=c("1-plex", "10-plex", "50-plex", "100-plex"),
#        lty=1:4)
#
# dev.off()

