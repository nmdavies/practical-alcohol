// do file for generating survival results by study
//PCMUS, UTAH, IPO-Porto Poland removed as not >90% complete data and <5 deaths
//Epic norfolk, Moffitt, Protect, QLD, search, STHM1, WUGS removed as no follow up data

//remove ProMPT from meta as <5 deaths 
// remove CPCS1/2, PPF from highgrade as <5 deaths
// remove CPCS1/2, EPIC, MAYO, MEC, PPF, ULM from low grade as <5 deaths

//This repeats the analysis for high grade defined as Gleason 8 or above:

cd "/Volumes/ecnmd/Post doc/Alcohol Practical Brunner"
use "working data/alcohol and prostate cancer_pheno_and_SNPs.dta", clear

drop if studyid=="EPIC-Norfolk"|studyid=="MOFFITT"|studyid=="ProtecT"|studyid=="QLD"|studyid=="SEARCH"|studyid=="STHM1"|studyid=="WUGS"|studyid=="PCMUS"|studyid=="UTAH"|studyid=="IPO-Porto"|studyid=="Poland"

//define variables for  high grade (gleason>=7), low grade (gleason<=6), local (t1-t2 or 'localised' on seer) and advanced (t3-t4 or advanced on seer)
//eg generate highgrade=(gleason!=.)
//eg replace highgrade="1" if highgrade>...
//eg replace highgrade="0" if caco==0
//eg replace highgrade="" if gleason<...

//generate highgrade
gen highgrade=(gleasonscore)
replace highgrade=1 if gleasonscore>=8 & gleasonscore!=.
replace highgrade=0 if caco==0
replace highgrade=. if gleasonscore<=7
replace highgrade=1 if gleasonscore==. & gleasonrange=="3"


//generate lowgrade
gen lowgrade=(gleasonscore)
replace lowgrade=1 if gleasonscore<=7
replace lowgrade=0 if caco==0
replace lowgrade=. if gleasonscore>=8&gleasonscore!=.
replace lowgrade=1 if gleasonscore==. & gleasonrange=="1"

//generate a binary variable for death
gen death2=datedeath

replace death2="1" if datedeath!="" & causedeath=="1"
replace death2="0" if datedeath=="" | causedeath!="1"

destring death2 , replace

//gen lengthfu
gen lengthfudied2=(datedeath2-datediag2) if datediag2!=. & datedeath2!=. & causedeath=="1"

replace lengthfudied2=(datelastfu2-datediag2) if lengthfudied2==.
replace lengthfudied2=. if caco==0

//cox regressions

foreach k in CAPS CPCS1 CPCS2 EPIC ESTHER FHCRC MAYO MCCS MEC /// 
	PPF-UNIS ProMPT TAMPERE UKGPCS ULM {
	preserve
	
	cap{
	keep if studyid=="`k'"

	stset lengthfudied2, id(icogs_sample_id) failure(death2) 

	stcox rs10156653_CT pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8
	regsave rs10156653_CT using "results/24715survivaladjustedSNP1caco_`k'", detail(all) ci pval replace

	stcox highgrade rs10156653_CT pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8
	regsave rs10156653_CT using "results/24715survivaladjustedSNP1highgrade_`k'", detail(all) ci pval replace

	stcox lowgrade rs10156653_CT pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8
	regsave rs10156653_CT using "results/24715survivaladjustedSNP1lowgrade_`k'", detail(all) ci pval replace

	ds rs*

	foreach i in `r(varlist)' {
		foreach j in caco highgrade lowgrade{
			cap:stcox `j' `i' pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8
			cap:regsave `i' using "results/24715survivaladjustedSNP1`j'_`k'" , detail(all) ci pval append
			}
		}

	}
	restore
	}


//to append all studies back together
use "results/24715survivaladjustedSNP1cacoCAPS", clear
gen studyid="CAPS"
foreach k in CPCS1 CPCS2 EPIC ESTHER FHCRC MAYO MCCS MEC /// 
	PPF-UNIS ProMPT TAMPERE UKGPCS ULM {

	cap:append using "results/24715survivaladjustedSNP1caco`k'"
	replace studyid="`k'" if studyid==""
	duplicates drop
	save "results/24715survivaladjustedSNP1caco_allstudies" , replace
	}
	
use "results/24715survivaladjustedSNP1highgrade_CAPS", clear
gen studyid="CAPS"
foreach k in CPCS1 CPCS2 EPIC ESTHER FHCRC MAYO MCCS MEC /// 
	PPF-UNIS ProMPT TAMPERE UKGPCS ULM {

	cap:append using "results/24715survivaladjustedSNP1highgrade_`k'"
	replace studyid="`k'" if studyid==""
	duplicates drop
	save "results/24715survivaladjustedSNP1highgrade_allstudies" , replace
	}

use "results/24715survivaladjustedSNP1lowgrade_CAPS", clear 
gen studyid="CAPS"
foreach k in CPCS1 CPCS2 EPIC ESTHER FHCRC MAYO MCCS MEC /// 
	PPF-UNIS ProMPT TAMPERE UKGPCS ULM {

	cap:append using "results/24715survivaladjustedSNP1lowgrade_`k'"
	replace studyid="`k'" if studyid==""
	duplicates drop
	save "results/24715survivaladjustedSNP1lowgrade_allstudies" , replace
	}

//Clean results
use "results/24715survivaladjustedSNP1caco_allstudies" , clear

rename var snp

drop if studyid=="CPCS1"|studyid=="CPCS2"|studyid=="PPF-UNIS"|studyid=="ProMPT"

matrix A=0,0,0,0,0,0
matrix B=0,0,0,0,0,0

levels snp
foreach n in `r(levels)' {
	metan coef ci_lower ci_upper if snp=="`n'", fixed eform lcols(studyid) title(`n') ///
		sub("Random effect association with survival in highgrade prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k') nograph
	//graph export "O:\alcohol and prostate cancer\documents\metan plots\random metan plot of `n' with survival in high grade pca 24715.pdf" , as(pdf) replace
	matrix A=A\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)
	
	metan coef ci_lower ci_upper if snp=="`n'", random eform lcols(studyid) title(`n') ///
		sub("Random effect association with survival in highgrade prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k') nograph
	//graph export "O:\alcohol and prostate cancer\documents\metan plots\random metan plot of `n' with survival in high grade pca 24715.pdf" , as(pdf) replace
	matrix B=B\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)
	}
	
	
	
use "results/24715survivaladjustedSNP1highgrade_allstudies" , clear

rename var snp

drop if studyid=="CPCS1"|studyid=="CPCS2"|studyid=="PPF-UNIS"|studyid=="ProMPT"

matrix C=0,0,0,0,0,0
matrix D=0,0,0,0,0,0

levels snp
foreach n in `r(levels)' {
	metan coef ci_lower ci_upper if snp=="`n'", fixed eform lcols(studyid) title(`n') ///
		sub("Random effect association with survival in highgrade prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k') nograph
	//graph export "O:\alcohol and prostate cancer\documents\metan plots\random metan plot of `n' with survival in high grade pca 24715.pdf" , as(pdf) replace
	matrix C=C\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)
	
	metan coef ci_lower ci_upper if snp=="`n'", random eform lcols(studyid) title(`n') ///
		sub("Random effect association with survival in highgrade prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k') nograph
	//graph export "O:\alcohol and prostate cancer\documents\metan plots\random metan plot of `n' with survival in high grade pca 24715.pdf" , as(pdf) replace
	matrix D=D\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)
	}

use "results/24715survivaladjustedSNP1lowgrade_allstudies" , clear

rename var snp

drop if studyid=="CPCS1"|studyid=="CPCS2"|studyid=="EPIC"|studyid=="MAYO"|studyid=="MEC"|studyid=="PPF-UNIS"|studyid=="ProMPT"|studyid=="ULM"

matrix E=0,0,0,0,0,0
matrix F=0,0,0,0,0,0

levels snp
foreach n in `r(levels)' {
	metan coef ci_lower ci_upper if snp=="`n'", fixed eform lcols(studyid) title(`n') ///
		sub("Random effect association with survival in lowgrade prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k') nograph
	//graph export "O:\alcohol and prostate cancer\documents\metan plots\random metan plot of `n' with survival in low grade pca 24715.pdf" , as(pdf) replace

	matrix E=E\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)
	
	metan coef ci_lower ci_upper if snp=="`n'", random eform lcols(studyid) title(`n') ///
		sub("Random effect association with survival in lowgrade prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k') nograph
	//graph export "O:\alcohol and prostate cancer\documents\metan plots\random metan plot of `n' with survival in low grade pca 24715.pdf" , as(pdf) replace

	matrix F=F\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)	
	}


	
//save the matrices to excel
//to see matrix type: matrix list A
//rename columns as per matrix code
foreach i in A B C D E F{
	matrix colnames `i'= "ES" "ci_low" "ci_upp" "p_z" "i_sq" "p_het"
	matrix rownames A= "blank" `"rs10156653_CT"' `"rs1042026_TC"' `"rs10814685_AG"' `"rs10973759_AG"' `"rs10973794_GA"' `"rs10973813_GA"' `"rs10973819_GA"' `"rs10973827_CT"' `"rs1154457_CT"' `"rs1159918_AC"' `"rs11792732_TC"' `"rs1229982_TG"' `"rs1229984_TC"' `"rs12336048_GA"' `"rs12441450_AG"' `"rs12505135_AG"' `"rs12910509_CT"' `"rs12914603_CG"' `"rs14226_GA"' `"rs1424482_CT"' `"rs1441817_CG"' `"rs1550573_CT"' `"rs1582620_CT"' `"rs1662037_GA"' `"rs166892_CT"' `"rs1693439_AG"' `"rs17526590_GA"' `"rs17583753_GA"' `"rs1899355_TC"' `"rs2009181_AG"' `"rs2067062_TG"' `"rs2173201_CA"' `"rs2181139_CT"' `"rs2381872_GC"' `"rs2642630_CG"' `"rs2642632_CA"' `"rs2642637_AG"' `"rs2654849_GT"' `"rs284792_CT"' `"rs2899613_TC"' `"rs3784264_TG"' `"rs3803433_GA"' `"rs4237255_GT"' `"rs4238327_CT"' `"rs4480202_TC"' `"rs4646557_TC"' `"rs4646572_TC"' `"rs4744675_GA"' `"rs4878203_CG"' `"rs4878811_CA"' `"rs4965726_CT"' `"rs4965741_CT"' `"rs698_TC"' `"rs7169439_GA"' `"rs7497343_AG"' `"rs769153_GT"' `"rs7859877_GA"' `"rs7865905_CT"' `"rs7866277_CT"' `"rs8041644_TC"' `"rs8041922_CA"' `"rs8187876_CT"' `"rs8187890_TC"' `"rs8187950_AG"' `"rs931635_AG"' `"rs987187_TC"' `"rs994772_CT"' `"rs998226_AG"'
	di "MATRIX ******`i'"
	matrix li `i'
	svmat `i'
	}
drop *6

//Create files to create Figures 1 & 2:

keep A* B* C* D* E* F*
drop if A1==.
drop in 1
gen snp=""

//Copy in RSids

save "working data/metan_results_survival",replace

	
	
