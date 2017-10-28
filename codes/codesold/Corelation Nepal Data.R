#see corelation of varibales# oir else go to line=199#
#library(corrgram)
#corrgram(myFD[,c("ld_dn",  "fw_dn",  "lh_dn",  "rd_dn",	"mb_ds",	"mu_ds",	"hw_ds",	"iw_pc")], order=TRUE, lower.panel=panel.ellipse,
upper.panel=panel.pts, text.panel=panel.txt,
diag.panel=panel.minmax, 
main="myFDin PC2/PC1 Order")
library(corrgram)
#corrgram(myFD[,c("ld_dn",  "fw_dn",  "lh_dn",  "rd_dn",	"mb_ds",	"mu_ds",	"hw_ds",	"iw_pc")], order=TRUE, lower.panel=panel.shade,
upper.panel=panel.pie, text.panel=panel.txt,
main="Car Milage Data in PC2/PC1 Order")

#corrgram(myFD[,c("fw_dn","du_dn","lh_dn","rd_dn","ir_pc","hp_dn","an_mt","an_pp","hh_dn", "ld_dn","iw_pc", "pd_pc", "fr_pc", "el_mt", "rn_fl", "rd_ds", "mb_ds", "rv_ds", "mu_ds", "hw_ds", "po_dn", "ps_cb", "ib_ds")], order=NULL, lower.panel=panel.shade,
upper.panel=NULL, text.panel=panel.txt,
main="Predictors")

#L<-cor(myFD[,c("fw_dn","du_dn","lh_dn","rd_dn","ir_pc","hp_dn","an_mt","an_pp","hh_dn", "ld_dn","iw_pc", "pd_pc", "fr_pc", "el_mt", "rn_fl", "rd_ds", "mb_ds", "rv_ds", "mu_ds", "hw_ds", "po_dn", "ps_cb", "ib_ds")], use="complete.obs", method="kendall") 

#capture.output(L,file="test.xls")

fit <- princomp(myFD, cor=TRUE)
summary(fit)
loadings(fit)
plot(fit,type="lines")
fit$scores
biplot(fit)

library(psych)
ct <- corr.test(myFD)