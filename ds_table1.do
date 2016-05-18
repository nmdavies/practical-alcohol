//Neil Davies 13/05/16
//Original do file created 11/8/14 by CB
//do file for generation of table 1

cd "/Volumes/ecnmd/Post doc/Alcohol Practical Brunner"
use "working data/alcohol and prostate cancer_pheno_and_SNPs.dta", clear

//to get list of studies and no cases/controls
tab studyid caco, matcell(table1)

//country of study entered manually into excel

//to get mean age at diag and psa at diag by studyid
//Correct the PSA in EPIC:
replace psadiag=psadiag/100 if studyid=="EPIC"

tabstat agediag if caco==1, by(studyid) stat(mean sd) save nototal
tabstatmat col23
matrix table1=table1,col23

tabstat psadiag if caco==1, by(studyid) stat(median p25 p75) save nototal
tabstatmat col23
matrix table1=table1,col23

//to get proportion with family history
gen familyhistory=(famhist=="1") if famhist!="U"&famhist!=""

tabstat familyhistory if caco==1, by(studyid)  save nototal
tabstatmat col23
matrix table1=table1,col23

//to get proportion with gleason score 8-10
gen gleason8_10=(gleasonscore>=8 & gleasonscore<=10) if gleasonscore!=.
replace gleason8_10=1 if gleasonscore==. & gleasonrange=="3"
replace gleason8_10=0 if gleasonscore==. & gleasonrange!="3"

tabstat gleason8_10 if caco==1, by(studyid) nototal save
tabstatmat col23
matrix table1=table1,col23
