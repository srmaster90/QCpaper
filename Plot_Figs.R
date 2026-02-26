# QC power figures for Unsihuay and Master (figures 1-4, supplemental figures 1-6)
# Plot figures using data generated in qcsim_simdata.R
# S Master, Children's Hospital of Philadelphia
# masters@chop.edu
# Last edit 2/6/26

library(qcsim)

load("Results/qcdata.RData")

# figure 1: 1-3s, 1-2s, full Westgard power curves

pdf("Figures/Figure1.pdf", width=9, height=9)

oldpar <- par(mfrow=c(2,2))

plot(pres.bias.3s.1)
abline(v=2, col=2, lty=3, lwd=2)

plot(pres.prec.3s.1)
abline(v=2, col=2, lty=3, lwd=2)

plot(pres.bias.3s.1)
lines(pres.bias.2s.1, col=2, lty=2)
lines(pres.bias.wg1.1, col=4, lty=4,)
legend("bottomright", legend=c(expression(1["3s"]),expression(1["2s"]),
                               expression(1["3s"]/2["2s"]/R["4s"]/4["1s"]/10["x"])),
       lty=c(1,2,4),col=c(1,2,4))

plot(pres.prec.3s.1)
lines(pres.prec.2s.1, col=2, lty=2)
lines(pres.prec.wg1.1, col=4, lty=4)
legend("bottomright", legend=c(expression(1["2s"]),expression(1["2s"]),
                               expression(1["3s"]/2["2s"]/R["4s"]/4["1s"]/10["x"])),
       lty=c(1,2,4),col=c(1,2,4))

par(oldpar)
dev.off()

# figure 2: multiplex 1-3s, full Westgard; bias

pdf("Figures/Figure2.pdf", width=9, height=4.5)
oldpar <- par(mfrow=c(1,2))
plot(pres.bias.3s2.1, main=expression(1["3s"]~"(2 levels)"))
abline(v=2, col=1, lty=3)
lines(pres.bias.3s2.10, lty=2, col=2)
lines(pres.bias.3s2.50, lty=3, col=3)
lines(pres.bias.3s2.100, lty=4, col=4)
lines(pres.bias.3s2.1000, lty=5, col=5)
legend("bottomright", legend=c("1-plex", "10-plex", "50-plex", "100-plex", "1000-plex"),
       lty=1:5,col=1:5)

plot(pres.bias.wg2.1, main=expression(1["3s"]/2["2s"]/R["4s"]/4["1s"]/10["x"]~"(2 levels)"))
abline(v=2, col=1, lty=3)
lines(pres.bias.wg2.10, lty=2, col=2)
lines(pres.bias.wg2.50, lty=3, col=3)
lines(pres.bias.wg2.100, lty=4, col=4)
lines(pres.bias.wg2.1000, lty=5, col=5)
legend("bottomright", legend=c("1-plex", "10-plex", "50-plex", "100-plex", "1000-plex"),
       lty=1:5,col=1:5)

par(oldpar)
dev.off()

# figure 3: Analytical performance (1-3s, 1-4s, 1-5s; 1-4s multiplexing), bias

# calculate SD for given power point for each curve (interpolate as necessary)
findthresh <- \(x, y, thresh) {
  ylow <- max(y[y < thresh])
  # interpolate
  idx <- which(y == ylow)
  if (idx == length(y)) return(NA)
  pctinc <- (thresh - y[idx]) / (y[idx+1] - y[idx])
  xd <- (x[idx+1] - x[idx]) * pctinc
  xout <- x[idx] + xd
  return(xout)
}

# now calc 80% power for fig3a curves
x3s <- findthresh(pres.bias.3s2.1$x, pres.bias.3s2.1$y, 0.8)
x4s <- findthresh(pres.bias.4s2.1$x, pres.bias.4s2.1$y, 0.8)
x5s <- findthresh(pres.bias.5s2.1$x, pres.bias.5s2.1$y, 0.8)

pdf("Figures/Figure3.pdf", width=9, height=4.5)
oldpar <- par(mfrow=c(1,2))
plot(pres.bias.5s2.1, col=4, lty=4, main=expression(1["3s"]~vs~1["4s"]~vs~1["5s"]~"(1-plex, 2 levels)"))
#annotation lines
lines(c(x3s,x3s), c(.8,-.05), lty=2, col=3)
lines(c(x4s,x4s), c(.8,-.05), lty=2, col=3)
lines(c(x5s,x5s), c(.8,-.05), lty=2, col=3)
lines(c(x3s,x5s), c(.8,.8), lty=2, col=3)
#...and more curves
lines(pres.bias.4s2.1, lty=2, col=2)
lines(pres.bias.3s2_x.1, lty=1, col=1)
legend("bottomright", legend=c(expression(1["3s"]), expression(1["4s"]), expression(1["5s"])),
       lty=c(1,2,4),col=c(1,2,4))

plot(pres.bias.4s2.1000, col=5, lty=5, main=expression(1["4s"]~"Multiplex Performance"))
lines(pres.bias.4s2.100, lty=4, col=4)
lines(pres.bias.4s2.50, lty=3, col=3)
lines(pres.bias.4s2.10, lty=2, col=2)
lines(pres.bias.4s2.1, lty=1, col=1)
legend("bottomright", legend=c("1-plex", "10-plex", "50-plex", "100-plex", "1000-plex"),
       lty=1:5,col=1:5)

par(oldpar)
dev.off()

# figure 4: T- and F-test QC rules (varying cutoff, varying replicates, 1000-plex comparison with single-plex Westgard)

pdf("Figures/Figure4.pdf", width=13.5, height=4.5)
oldpar <- par(mfrow=c(1,3))
plot(pres.bias.tf2.50_25ref_10rep_1en3, col=1, lty=1,main="50-plex, TF with varying P cutoff")
lines(pres.bias.tf2.50_25ref_10rep_5en4, col=2, lty=2)
lines(pres.bias.tf2.50_25ref_10rep_5en5, col=4, lty=4)
legend("bottomright", legend=c(expression("P=1x10"^-3), expression("P=5x10"^-4), expression("P=5x10"^-5)),
       lty=c(1,2,4),col=c(1,2,4))

plot(pres.bias.tf2.50_25ref_3rep_5en5, col=1, lty=1,main="50-plex, TF with varying QC replicates")
lines(pres.bias.tf2.50_25ref_5rep_5en5, col=2, lty=2)
lines(pres.bias.tf2.50_25ref_10rep_5en5, col=4, lty=4)
legend("bottomright", legend=c("3 replicates", "5 replicates", "10 replicates"),
       lty=c(1,2,4),col=c(1,2,4))

plot(pres.bias.tf2.1000_25ref_10rep_2.7en6, col=1, lty=1,main="1000-plex TF compared with traditional Westgard 1-plex")
lines(pres.bias.wg2.1, col=4, lty=4)
lines(pres.bias.3s2.1, col=3, lty=2)
legend("bottomright", legend=c(expression("1000-plex TF (P<2.7x10"^-6*")"), expression("1-plex "*1["3s"]), expression("1-plex "*1["3s"]/2["2s"]/R["4s"]/4["1s"]/10["x"])),
       lty=c(1,2,4),col=c(1,3,4))


par(oldpar)
dev.off()


# SUPPLEMENTAL FIGURES

# supplemental figure 1: bias and imprecision
pdf("Figures/Supp_Figure1.pdf", width=10, height=5)
oldpar <- par(mfrow=c(1,2))
x <- seq(from = -5, to = 5, length.out=300)
plot(x, dnorm(x), type="l", lwd=2, xaxt="n", xlab="Z", ylab="Probability", main="2 SD shift in assay bias")
lines(x, dnorm(x, mean=2), type="l", lty=2, col=2, lwd=2)
axis(side=1, at=seq(-4, 4, by=2))

x <- seq(from = -5, to = 5, length.out=300)
plot(x, dnorm(x), type="l", lwd=2, xaxt="n", xlab="Z", ylab="Probability", main="2x shift in assay imprecision")
lines(x, dnorm(x, sd=2), type="l", lty=2, col=2, lwd=2)
axis(side=1, at=seq(-4, 4, by=2))
par(oldpar)
dev.off()

# supplemental figure 2: multiplex 1-3s, full Westgard; imprecision

pdf("Figures/Supp_Figure2.pdf", width=9, height=4.5)
oldpar <- par(mfrow=c(1,2))
plot(pres.prec.3s2.1, main=expression(1["3s"]~"(2 levels)"))
lines(pres.prec.3s2.10, lty=2, col=2)
lines(pres.prec.3s2.50, lty=3, col=3)
lines(pres.prec.3s2.100, lty=4, col=4)
lines(pres.prec.3s2.1000, lty=5, col=5)
legend("bottomright", legend=c("1-plex", "10-plex", "50-plex", "100-plex", "1000-plex"),
       lty=1:5,col=1:5)

plot(pres.prec.wg2.1, main=expression(1["3s"]/2["2s"]/R["4s"]/4["1s"]/10["x"]~"(2 levels)"))
lines(pres.prec.wg2.10, lty=2, col=2)
lines(pres.prec.wg2.50, lty=3, col=3)
lines(pres.prec.wg2.100, lty=4, col=4)
lines(pres.prec.wg2.1000, lty=5, col=5)
legend("bottomright", legend=c("1-plex", "10-plex", "50-plex", "100-plex", "1000-plex"),
       lty=1:5,col=1:5)

par(oldpar)
dev.off()

# supplemental figure 3: Analytical performance (1-3s, 1-4s, 1-5s; 1-4s multiplexing), imprecision (see figure 3)
pdf("Figures/Supp_Figure3.pdf", width=9, height=4.5)
oldpar <- par(mfrow=c(1,2))
plot(pres.prec.5s2.1, col=4, lty=4, main=expression(1["3s"]~vs~1["4s"]~vs~1["5s"]~"(1-plex, 2 levels)"))
lines(pres.prec.4s2.1, lty=2, col=2)
lines(pres.prec.3s2_x.1, lty=1, col=1)
legend("bottomright", legend=c(expression(1["3s"]), expression(1["4s"]), expression(1["5s"])),
       lty=c(1,2,4),col=c(1,2,4))

plot(pres.prec.4s2.1000, col=5, lty=5, main=expression(1["4s"]~"Multiplex Performance"))
lines(pres.prec.4s2.100, lty=4, col=4)
lines(pres.prec.4s2.50, lty=3, col=3)
lines(pres.prec.4s2.10, lty=2, col=2)
lines(pres.prec.4s2.1, lty=1, col=1)
legend("bottomright", legend=c("1-plex", "10-plex", "50-plex", "100-plex", "1000-plex"),
       lty=1:5,col=1:5)

par(oldpar)
dev.off()

# supplemental figure 4: 1-5s multiplex performance
pdf("Figures/Supp_Figure4.pdf", width=10, height=5)
oldpar <- par(mfrow=c(1,2))
plot(pres.bias.5s2.1, col=1, lty=5, main=expression(1["5s"]~"Multiplex Performance (Bias)"))
lines(pres.bias.5s2.10, lty=2, col=4)
lines(pres.bias.5s2.50, lty=3, col=3)
lines(pres.bias.5s2.100, lty=4, col=2)
lines(pres.bias.5s2.1000, lty=5, col=1)
legend("bottomright", legend=c("1-plex", "10-plex", "50-plex", "100-plex", "1000-plex"),
       lty=1:5,col=1:5)

plot(pres.prec.5s2.1, col=1, lty=1, main=expression(1["5s"]~"Multiplex Performance (Imprecision)"))
lines(pres.prec.5s2.10, lty=2, col=2)
lines(pres.prec.5s2.50, lty=3, col=3)
lines(pres.prec.5s2.100, lty=4, col=4)
lines(pres.prec.5s2.1000, lty=5, col=5)
legend("bottomright", legend=c("1-plex", "10-plex", "50-plex", "100-plex", "1000-plex"),
       lty=1:5,col=1:5)

par(oldpar)
dev.off()

# supplemental figure 5: T and F alone
pdf("Figures/Supp_Figure5.pdf", width=9, height=4.5)
oldpar <- par(mfrow=c(1,2))
plot(pres.bias.tf2.50_25ref_10rep_5en5, col=1, lty=1,main=expression("T and F test Bias (50-plex, P=5x10"^-5*")"))
lines(pres.bias.t2.50_25ref_10rep_5en5, col=2, lty=2)
lines(pres.bias.f2.50_25ref_10rep_5en5, col=4, lty=4)
legend(3, 0.4, legend=c("T+F", "T only", "F only"),
       lty=c(1,2,4),col=c(1,2,4))

plot(pres.prec.tf2.50_25ref_10rep_5en5, col=1, lty=1,main=expression("T and F test Imprecision (50-plex, P=5x10"^-5*")"))
lines(pres.prec.t2.50_25ref_10rep_5en5, col=2, lty=2)
lines(pres.prec.f2.50_25ref_10rep_5en5, col=4, lty=4)
legend(4, 0.4, legend=c("T+F", "T only", "F only"),
       lty=c(1,2,4),col=c(1,2,4))

par(oldpar)
dev.off()

plot(pres.prec.t2.50_25ref_10rep_5en5, col=1, lty=1,main=expression("T-test only (50-plex, P=5x10"^-5*"): Imprecision"))
plot(pres.bias.f2.50_25ref_10rep_5en5, col=1, lty=1,main=expression("F-test only (50-plex, P=5x10"^-5*"): Bias"))
plot(pres.prec.f2.50_25ref_10rep_5en5, col=1, lty=1,main=expression("F-test only (50-plex, P=5x10"^-5*"): Imprecision"))

# supplemental figure 6: TF imprecision curves
pdf("Figures/Supp_Figure6.pdf", width=13.5, height=4.5)
oldpar <- par(mfrow=c(1,3))
plot(pres.prec.tf2.50_25ref_10rep_1en3, col=1, lty=1,main="50-plex, TF with varying P cutoff")
lines(pres.prec.tf2.50_25ref_10rep_5en4, col=2, lty=2)
lines(pres.prec.tf2.50_25ref_10rep_5en5, col=4, lty=4)
legend("bottomright", legend=c(expression("P=1x10"^-3), expression("P=5x10"^-4), expression("P=5x10"^-5)),
       lty=c(1,2,4),col=c(1,2,4))

plot(pres.prec.tf2.50_25ref_3rep_5en5, col=1, lty=1,main="50-plex, TF with varying QC replicates")
lines(pres.prec.tf2.50_25ref_5rep_5en5, col=2, lty=2)
lines(pres.prec.tf2.50_25ref_10rep_5en5, col=4, lty=4)
legend("bottomright", legend=c("3 replicates", "5 replicates", "10 replicates"),
       lty=c(1,2,4),col=c(1,2,4))

plot(pres.prec.tf2.1000_25ref_10rep_2.7en6, col=1, lty=1,main="1000-plex TF compared with traditional Westgard 1-plex")
lines(pres.prec.wg2.1, col=4, lty=4)
lines(pres.prec.3s2.1, col=3, lty=2)
legend("bottomright", legend=c(expression("1000-plex TF (P<2.7x10"^-6*")"), expression("1-plex "*1["3s"]), expression("1-plex "*1["3s"]/2["2s"]/R["4s"]/4["1s"]/10["x"])),
       lty=c(1,2,4),col=c(1,3,4))

par(oldpar)
dev.off()


