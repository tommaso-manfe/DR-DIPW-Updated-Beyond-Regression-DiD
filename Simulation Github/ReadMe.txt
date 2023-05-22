The name of the experiments closely follows the paper.
The estimators in the paper refer to the following functions:

IPW:std_ipw_did_rc.R (Abadie, 2005)
DIPW:dipw.R (author's work)
OR:reg_did_rc.R (Heckmann et al, 1997)
DMLDiD:DMLDiD.R (Chang, 2020)
DRDiD: drdid_rc.R (Sant'Anna and Zhao, 2020)
DR-DIPW:  drdid_rc_hayek.R (author's work)
IMP DRDiD: drdid_imp_rc.R (Sant'Anna and Zhao, 2020)
IMP DR-DIPW: drdid_imp_rc_hayek.R (author's work)
LASSO DR-DIPW: lasso_drdid_rc_hayek.R (author's work)
RF DR-DIPW: randforest_drdid_rc_hayek.R (author's work)

TWFE methods are instead run with the lm() function.
