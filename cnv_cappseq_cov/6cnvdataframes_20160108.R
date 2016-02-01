################################################################################
# Checking selected dataframe and dealing with NAs
################################################################################
colvec_titan <- cnsPT_EEx.df
colvec_titan[is.na(colvec_titan)] <- '100'
ncol(colvec_titan)
colvec_titan

colvec_titan_short <- cnsPT_EEx_short.df
colvec_titan_short[is.na(colvec_titan_short)] <- '100'
ncol(colvec_titan_short)
colvec_titan_short

colvec_seq <- cns_seqPT_EEx.df 
colvec_seq[is.na(colvec_seq)] <- '100'
ncol(colvec_seq)
colvec_seq

patientcovz.df
ncol(patientcovz.df)
colnames(patientcovz.df[2])

zLog_tumNor.df
ncol(zLog_tumNor.df)

zLogcappcov_tumNor.df
ncol(zLog_tumNor.df)

ptbafPT_AR_Ex.df
ncol(ptbafPT_AR_Ex.df)

ptbafPT_AR_EEX.df
ncol(ptbafPT_AR_EEX.df)

ptbafPT_CNS_Ex.df
ncol(ptbafPT_CNS_Ex.df)

ptbafPT_CNS_EEx.df
ncol(ptbafPT_CNS_EEx.df)




