//Neil Davies 13/05/16
//do file for generating results by individual study
//originally created 10/9/14 by CB
//note Prompt and WUGS no controls so cannot regress

cd "/Volumes/ecnmd/Post doc/Alcohol Practical Brunner"
use "working data/alcohol and prostate cancer_pheno_and_SNPs.dta", clear

//define variables for  high grade (gleason>=7), low grade (gleason<=6)
//eg generate highgrade=(gleason!=.)
//eg replace highgrade="1" if highgrade>...
//eg replace highgrade="0" if caco==0
//eg replace highgrade="" if gleason<...

//generate highgrade
gen highgrade=(gleasonscore)
replace highgrade=1 if gleasonscore>=7 & gleasonscore!=.
replace highgrade=0 if caco==0
replace highgrade=. if gleasonscore<=6
replace highgrade=1 if gleasonscore==. & gleasonrange=="3"


//generate lowgrade
gen lowgrade=(gleasonscore)
replace lowgrade=1 if gleasonscore<=6
replace lowgrade=0 if caco==0
replace lowgrade=. if gleasonscore>=7&gleasonscore!=.
replace lowgrade=1 if gleasonscore==. & gleasonrange=="1"


//logistic (outcome) (exposure)

foreach k in CAPS CPCS1 CPCS2 EPIC EPIC-Norfolk ESTHER FHCRC IPO-Porto MAYO MCCS MEC /// 
	MOFFITT PCMUS Poland PPF-UNIS ProtecT QLD SEARCH STHM1 TAMPERE UKGPCS ULM UTAH {
	
	preserve
	keep if studyid=="`k'"
	
	logistic caco rs10156653_CT pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8
	regsave rs10156653_CT using "results/prostatecancerSNP1cacoadjusted`k'", detail(all) ci pval replace


	cap:logistic highgrade rs10156653_CT pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8
	cap:regsave rs10156653_CT using "results/prostatecancerSNP1highgradeadjusted`k'", detail(all) ci pval replace


	cap:logistic lowgrade rs10156653_CT pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8
	cap:regsave rs10156653_CT using "results/prostatecancerSNP1lowgradeadjusted`k'", detail(all) ci pval replace

	ds rs* 

	foreach i in `r(varlist)' {
		foreach j in caco highgrade lowgrade{
			cap:logistic `j' `i' pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8
			cap:regsave `i' using "results/prostatecancerSNP1`j'adjusted`k'" , detail(all) ci pval append
			}
		}

	restore
	}

//open results and drop in 1 to get rid of duplicate first line

use "results/prostatecancerSNP1cacoadjustedCAPS", clear
gen studyid="CAPS"
foreach k in CPCS1 CPCS2 EPIC EPIC-Norfolk ESTHER FHCRC IPO-Porto MAYO MCCS MEC /// 
	MOFFITT PCMUS Poland PPF-UNIS ProtecT QLD SEARCH STHM1 TAMPERE UKGPCS ULM UTAH {

	append using "results/prostatecancerSNP1cacoadjusted`k'"
	replace studyid="`k'" if studyid==""
	duplicates drop
	save "results/prostatecancerSNP1cacoadjusted_allstudies" , replace
	}

use "results/prostatecancerSNP1highgradeadjustedCAPS", clear
gen studyid="CAPS"
foreach k in CPCS1 CPCS2 EPIC EPIC-Norfolk ESTHER FHCRC IPO-Porto MAYO MCCS MEC /// 
	MOFFITT PCMUS Poland PPF-UNIS ProtecT QLD SEARCH STHM1 TAMPERE UKGPCS ULM UTAH {

	append using "results/prostatecancerSNP1highgradeadjusted`k'"
	replace studyid="`k'" if studyid==""
	duplicates drop
	save "results/prostatecancerSNP1highgradeadjusted_allstudies" , replace
	}

use "results/prostatecancerSNP1lowgradeadjustedCAPS", clear 
gen studyid="CAPS"
foreach k in CPCS1 CPCS2 EPIC EPIC-Norfolk ESTHER FHCRC IPO-Porto MAYO MCCS MEC /// 
	MOFFITT PCMUS Poland PPF-UNIS ProtecT QLD SEARCH STHM1 TAMPERE UKGPCS ULM UTAH {

	append using "results/prostatecancerSNP1lowgradeadjusted`k'"
	replace studyid="`k'" if studyid==""

	duplicates drop
	save "results/prostatecancerSNP1lowgradeadjusted_allstudies" , replace
	}
	
//Meta-analyse all the results:

use "results/prostatecancerSNP1cacoadjusted_allstudies" , clear

replace var=substr(var,6,.)
rename var snp

matrix A=0,0,0,0,0,0

levels snp
foreach n in `r(levels)' {
	metan coef ci_lower ci_upper if snp=="`n'", fixed eform lcols(studyid) title(`n') sub("Fixed effect association with any prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k') nograph
	//graph export "results/metan_plots/fixed metan plot of `n' with adjust any pca.pdf" , as(pdf) replace

	matrix A=A\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)
	}
//random
matrix F=0,0,0,0,0,0

levels snp
foreach n in `r(levels)' {
	metan coef ci_lower ci_upper if snp=="`n'", random eform lcols(studyid) title(`n') sub("Random effects association with any prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k')  nograph
	//graph export "results/metan_plots/random metan plot of `n' with adjust any pca.pdf" , as(pdf) replace

	matrix F=F\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)
	}

//fixed	
use "results/prostatecancerSNP1highgradeadjusted_allstudies" , clear
replace var=substr(var,11,.)
rename var snp

matrix B=0,0,0,0,0,0

levels snp
foreach n in `r(levels)' {
	metan coef ci_lower ci_upper if snp=="`n'", fixed eform lcols(studyid) title(`n') sub("Fixed effects association with highgrade prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k') nograph 
	//graph export "results/metan_plots/fixed metan plot of `n' with adjust high grade pca.pdf" , as(pdf) replace

	matrix B=B\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)
	}
//random
matrix G=0,0,0,0,0,0

levels snp
foreach n in `r(levels)' {
	metan coef ci_lower ci_upper if snp=="`n'", random eform lcols(studyid) title(`n') sub("Random effects association with highgrade prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k')  nograph
	//graph export "results/metan_plots/random metan plot of `n' with adjust high grade pca.pdf" , as(pdf) replace

	matrix G=G\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)
	}

	
	
//fixed
use "results/prostatecancerSNP1lowgradeadjusted_allstudies" , clear
replace var=substr(var,10,.)
rename var snp

matrix C=0,0,0,0,0,0

levels snp
foreach n in `r(levels)' {
	metan coef ci_lower ci_upper if snp=="`n'", fixed eform lcols(studyid) title(`n') sub("Fixed effect association with lowgrade prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k')  nograph
	//graph export "results/metan_plots/fixed metan plot of `n' with adjust low grade pca.pdf" , as(pdf) replace

	matrix C=C\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)
	}

//random
matrix H=0,0,0,0,0,0

levels snp
foreach n in `r(levels)' {
	metan coef ci_lower ci_upper if snp=="`n'", random eform lcols(studyid) title(`n') sub("Random effects association with lowgrade prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k')  nograph
	//graph export "results/metan_plots/random metan plot of `n' with adjust low grade pca.pdf" , as(pdf) replace

	matrix H=H\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)
	}
	
	
//save the matrices to excel
//to see matrix type: matrix list A
//rename columns as per matrix code
matrix colnames A= "ES" "ci_low" "ci_upp" "p_z" "i_sq" "p_het"
matrix colnames B= "ES" "ci_low" "ci_upp" "p_z" "i_sq" "p_het"
matrix colnames C= "ES" "ci_low" "ci_upp" "p_z" "i_sq" "p_het"
matrix colnames D= "ES" "ci_low" "ci_upp" "p_z" "i_sq" "p_het"
matrix colnames E= "ES" "ci_low" "ci_upp" "p_z" "i_sq" "p_het"

//rename rows as each SNP
//`"rs1042026_TC"' `"rs1154457_CT"' `"rs1159918_AC"' `"rs1229982_TG"' `"rs1229984_TC"' `"rs12505135_AG"' `"rs1662037_GA"' `"rs166892_CT"' `"rs1693439_AG"' `"rs17526590_GA"' `"rs17583753_GA"' `"rs2009181_AG"' `"rs2173201_CA"' `"rs2654849_GT"' `"rs284792_CT"' `"rs698_TC"' `"rs931635_AG"' `"rs994772_CT"' `"rs10814685_AG"' `"rs10973759_AG"' `"rs10973794_GA"' `"rs10973813_GA"' `"rs10973819_GA"' `"rs10973827_CT"' `"rs11792732_TC"' `"rs12336048_GA"' `"rs2181139_CT"' `"rs2381872_GC"' `"rs4480202_TC"' `"rs4878203_CG"' `"rs4878811_CA"' `"rs769153_GT"' `"rs7859877_GA"' `"rs7865905_CT"' `"rs7866277_CT"' `"rs987187_TC"' `"rs998226_AG"' `"rs10156653_CT"' `"rs1424482_CT"' `"rs1582620_CT"' `"rs4237255_GT"' `"rs4744675_GA"' `"rs8187876_CT"' `"rs8187890_TC"' `"rs8187950_AG"' `"rs12910509_CT"' `"rs12914603_CG"' `"rs1441817_CG"' `"rs1550573_CT"' `"rs1899355_TC"' `"rs2067062_TG"' `"rs2642630_CG"' `"rs2642632_CA"' `"rs2642637_AG"' `"rs2899613_TC"' `"rs3784264_TG"' `"rs4238327_CT"' `"rs4646557_TC"' `"rs4646572_TC"' `"rs7169439_GA"' `"rs8041644_TC"' `"rs8041922_CA"' `"rs12441450_AG"' `"rs14226_GA"' `"rs3803433_GA"' `"rs4965726_CT"' `"rs4965741_CT"' `"rs7497343_AG"'

matrix rownames A= "blank" `"rs1042026_TC"' `"rs1154457_CT"' `"rs1159918_AC"' `"rs1229982_TG"' `"rs1229984_TC"' `"rs12505135_AG"' `"rs1662037_GA"' `"rs166892_CT"' `"rs1693439_AG"' `"rs17526590_GA"' `"rs17583753_GA"' `"rs2009181_AG"' `"rs2173201_CA"' `"rs2654849_GT"' `"rs284792_CT"' `"rs698_TC"' `"rs931635_AG"' `"rs994772_CT"' `"rs10814685_AG"' `"rs10973759_AG"' `"rs10973794_GA"' `"rs10973813_GA"' `"rs10973819_GA"' `"rs10973827_CT"' `"rs11792732_TC"' `"rs12336048_GA"' `"rs2181139_CT"' `"rs2381872_GC"' `"rs4480202_TC"' `"rs4878203_CG"' `"rs4878811_CA"' `"rs769153_GT"' `"rs7859877_GA"' `"rs7865905_CT"' `"rs7866277_CT"' `"rs987187_TC"' `"rs998226_AG"' `"rs10156653_CT"' `"rs1424482_CT"' `"rs1582620_CT"' `"rs4237255_GT"' `"rs4744675_GA"' `"rs8187876_CT"' `"rs8187890_TC"' `"rs8187950_AG"' `"rs12910509_CT"' `"rs12914603_CG"' `"rs1441817_CG"' `"rs1550573_CT"' `"rs1899355_TC"' `"rs2067062_TG"' `"rs2642630_CG"' `"rs2642632_CA"' `"rs2642637_AG"' `"rs2899613_TC"' `"rs3784264_TG"' `"rs4238327_CT"' `"rs4646557_TC"' `"rs4646572_TC"' `"rs7169439_GA"' `"rs8041644_TC"' `"rs8041922_CA"' `"rs12441450_AG"' `"rs14226_GA"' `"rs3803433_GA"' `"rs4965726_CT"' `"rs4965741_CT"' `"rs7497343_AG"'
matrix rownames B= "blank" `"rs1042026_TC"' `"rs1154457_CT"' `"rs1159918_AC"' `"rs1229982_TG"' `"rs1229984_TC"' `"rs12505135_AG"' `"rs1662037_GA"' `"rs166892_CT"' `"rs1693439_AG"' `"rs17526590_GA"' `"rs17583753_GA"' `"rs2009181_AG"' `"rs2173201_CA"' `"rs2654849_GT"' `"rs284792_CT"' `"rs698_TC"' `"rs931635_AG"' `"rs994772_CT"' `"rs10814685_AG"' `"rs10973759_AG"' `"rs10973794_GA"' `"rs10973813_GA"' `"rs10973819_GA"' `"rs10973827_CT"' `"rs11792732_TC"' `"rs12336048_GA"' `"rs2181139_CT"' `"rs2381872_GC"' `"rs4480202_TC"' `"rs4878203_CG"' `"rs4878811_CA"' `"rs769153_GT"' `"rs7859877_GA"' `"rs7865905_CT"' `"rs7866277_CT"' `"rs987187_TC"' `"rs998226_AG"' `"rs10156653_CT"' `"rs1424482_CT"' `"rs1582620_CT"' `"rs4237255_GT"' `"rs4744675_GA"' `"rs8187876_CT"' `"rs8187890_TC"' `"rs8187950_AG"' `"rs12910509_CT"' `"rs12914603_CG"' `"rs1441817_CG"' `"rs1550573_CT"' `"rs1899355_TC"' `"rs2067062_TG"' `"rs2642630_CG"' `"rs2642632_CA"' `"rs2642637_AG"' `"rs2899613_TC"' `"rs3784264_TG"' `"rs4238327_CT"' `"rs4646557_TC"' `"rs4646572_TC"' `"rs7169439_GA"' `"rs8041644_TC"' `"rs8041922_CA"' `"rs12441450_AG"' `"rs14226_GA"' `"rs3803433_GA"' `"rs4965726_CT"' `"rs4965741_CT"' `"rs7497343_AG"'
matrix rownames C= "blank" `"rs1042026_TC"' `"rs1154457_CT"' `"rs1159918_AC"' `"rs1229982_TG"' `"rs1229984_TC"' `"rs12505135_AG"' `"rs1662037_GA"' `"rs166892_CT"' `"rs1693439_AG"' `"rs17526590_GA"' `"rs17583753_GA"' `"rs2009181_AG"' `"rs2173201_CA"' `"rs2654849_GT"' `"rs284792_CT"' `"rs698_TC"' `"rs931635_AG"' `"rs994772_CT"' `"rs10814685_AG"' `"rs10973759_AG"' `"rs10973794_GA"' `"rs10973813_GA"' `"rs10973819_GA"' `"rs10973827_CT"' `"rs11792732_TC"' `"rs12336048_GA"' `"rs2181139_CT"' `"rs2381872_GC"' `"rs4480202_TC"' `"rs4878203_CG"' `"rs4878811_CA"' `"rs769153_GT"' `"rs7859877_GA"' `"rs7865905_CT"' `"rs7866277_CT"' `"rs987187_TC"' `"rs998226_AG"' `"rs10156653_CT"' `"rs1424482_CT"' `"rs1582620_CT"' `"rs4237255_GT"' `"rs4744675_GA"' `"rs8187876_CT"' `"rs8187890_TC"' `"rs8187950_AG"' `"rs12910509_CT"' `"rs12914603_CG"' `"rs1441817_CG"' `"rs1550573_CT"' `"rs1899355_TC"' `"rs2067062_TG"' `"rs2642630_CG"' `"rs2642632_CA"' `"rs2642637_AG"' `"rs2899613_TC"' `"rs3784264_TG"' `"rs4238327_CT"' `"rs4646557_TC"' `"rs4646572_TC"' `"rs7169439_GA"' `"rs8041644_TC"' `"rs8041922_CA"' `"rs12441450_AG"' `"rs14226_GA"' `"rs3803433_GA"' `"rs4965726_CT"' `"rs4965741_CT"' `"rs7497343_AG"'
matrix rownames D= "blank" `"rs1042026_TC"' `"rs1154457_CT"' `"rs1159918_AC"' `"rs1229982_TG"' `"rs1229984_TC"' `"rs12505135_AG"' `"rs1662037_GA"' `"rs166892_CT"' `"rs1693439_AG"' `"rs17526590_GA"' `"rs17583753_GA"' `"rs2009181_AG"' `"rs2173201_CA"' `"rs2654849_GT"' `"rs284792_CT"' `"rs698_TC"' `"rs931635_AG"' `"rs994772_CT"' `"rs10814685_AG"' `"rs10973759_AG"' `"rs10973794_GA"' `"rs10973813_GA"' `"rs10973819_GA"' `"rs10973827_CT"' `"rs11792732_TC"' `"rs12336048_GA"' `"rs2181139_CT"' `"rs2381872_GC"' `"rs4480202_TC"' `"rs4878203_CG"' `"rs4878811_CA"' `"rs769153_GT"' `"rs7859877_GA"' `"rs7865905_CT"' `"rs7866277_CT"' `"rs987187_TC"' `"rs998226_AG"' `"rs10156653_CT"' `"rs1424482_CT"' `"rs1582620_CT"' `"rs4237255_GT"' `"rs4744675_GA"' `"rs8187876_CT"' `"rs8187890_TC"' `"rs8187950_AG"' `"rs12910509_CT"' `"rs12914603_CG"' `"rs1441817_CG"' `"rs1550573_CT"' `"rs1899355_TC"' `"rs2067062_TG"' `"rs2642630_CG"' `"rs2642632_CA"' `"rs2642637_AG"' `"rs2899613_TC"' `"rs3784264_TG"' `"rs4238327_CT"' `"rs4646557_TC"' `"rs4646572_TC"' `"rs7169439_GA"' `"rs8041644_TC"' `"rs8041922_CA"' `"rs12441450_AG"' `"rs14226_GA"' `"rs3803433_GA"' `"rs4965726_CT"' `"rs4965741_CT"' `"rs7497343_AG"'
matrix rownames E= "blank" `"rs1042026_TC"' `"rs1154457_CT"' `"rs1159918_AC"' `"rs1229982_TG"' `"rs1229984_TC"' `"rs12505135_AG"' `"rs1662037_GA"' `"rs166892_CT"' `"rs1693439_AG"' `"rs17526590_GA"' `"rs17583753_GA"' `"rs2009181_AG"' `"rs2173201_CA"' `"rs2654849_GT"' `"rs284792_CT"' `"rs698_TC"' `"rs931635_AG"' `"rs994772_CT"' `"rs10814685_AG"' `"rs10973759_AG"' `"rs10973794_GA"' `"rs10973813_GA"' `"rs10973819_GA"' `"rs10973827_CT"' `"rs11792732_TC"' `"rs12336048_GA"' `"rs2181139_CT"' `"rs2381872_GC"' `"rs4480202_TC"' `"rs4878203_CG"' `"rs4878811_CA"' `"rs769153_GT"' `"rs7859877_GA"' `"rs7865905_CT"' `"rs7866277_CT"' `"rs987187_TC"' `"rs998226_AG"' `"rs10156653_CT"' `"rs1424482_CT"' `"rs1582620_CT"' `"rs4237255_GT"' `"rs4744675_GA"' `"rs8187876_CT"' `"rs8187890_TC"' `"rs8187950_AG"' `"rs12910509_CT"' `"rs12914603_CG"' `"rs1441817_CG"' `"rs1550573_CT"' `"rs1899355_TC"' `"rs2067062_TG"' `"rs2642630_CG"' `"rs2642632_CA"' `"rs2642637_AG"' `"rs2899613_TC"' `"rs3784264_TG"' `"rs4238327_CT"' `"rs4646557_TC"' `"rs4646572_TC"' `"rs7169439_GA"' `"rs8041644_TC"' `"rs8041922_CA"' `"rs12441450_AG"' `"rs14226_GA"' `"rs3803433_GA"' `"rs4965726_CT"' `"rs4965741_CT"' `"rs7497343_AG"'

//Repeat for the survival results:
//random
use "results/prostatecancerSNP1cacoadjusted_allstudies" , clear

replace var=substr(var,6,.)
rename var snp

matrix F=0,0,0,0,0,0

levels snp
foreach n in `r(levels)' {
	metan coef ci_lower ci_upper if snp=="`n'", random eform lcols(studyid) title(`n') sub("Random effects association with any prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k') nograph
//	graph export "results/metan_plots/random metan plot of `n' with adjust any pca.pdf" , as(pdf) replace

	matrix F=F\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)
	}


use "results/prostatecancerSNP1highgradeadjusted_allstudies" , clear
replace var=substr(var,11,.)
rename var snp

matrix G=0,0,0,0,0,0

levels snp
foreach n in `r(levels)' {
	metan coef ci_lower ci_upper if snp=="`n'", random eform lcols(studyid) title(`n') sub("Random effects association with highgrade prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k')   nograph
//	graph export "results/metan_plots/random metan plot of `n' with adjust high grade pca.pdf" , as(pdf) replace

	matrix G=G\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)
	}

use "results/prostatecancerSNP1lowgradeadjusted_allstudies" , clear
replace var=substr(var,10,.)
rename var snp

matrix H=0,0,0,0,0,0

levels snp
foreach n in `r(levels)' {
	metan coef ci_lower ci_upper if snp=="`n'", random eform lcols(studyid) title(`n') sub("Random effects association with lowgrade prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k')  nograph
	//graph export "results/metan_plots/random metan plot of `n' with adjust low grade pca.pdf" , as(pdf) replace

	matrix H=H\r(ES),r(ci_low),r(ci_upp),r(p_z), r(i_sq),r(p_het)
	}
	
	
//save the matrices to excel
//to see matrix type: matrix list A
//rename columns as per matrix code
matrix colnames F= "ES" "ci_low" "ci_upp" "p_z" "i_sq" "p_het"
matrix colnames G= "ES" "ci_low" "ci_upp" "p_z" "i_sq" "p_het"
matrix colnames H= "ES" "ci_low" "ci_upp" "p_z" "i_sq" "p_het"
matrix colnames I= "ES" "ci_low" "ci_upp" "p_z" "i_sq" "p_het"
matrix colnames J= "ES" "ci_low" "ci_upp" "p_z" "i_sq" "p_het"

//rename rows as each SNP

matrix rownames F= "blank" `"rs1042026_TC"' `"rs1154457_CT"' `"rs1159918_AC"' `"rs1229982_TG"' `"rs1229984_TC"' `"rs12505135_AG"' `"rs1662037_GA"' `"rs166892_CT"' `"rs1693439_AG"' `"rs17526590_GA"' `"rs17583753_GA"' `"rs2009181_AG"' `"rs2173201_CA"' `"rs2654849_GT"' `"rs284792_CT"' `"rs698_TC"' `"rs931635_AG"' `"rs994772_CT"' `"rs10814685_AG"' `"rs10973759_AG"' `"rs10973794_GA"' `"rs10973813_GA"' `"rs10973819_GA"' `"rs10973827_CT"' `"rs11792732_TC"' `"rs12336048_GA"' `"rs2181139_CT"' `"rs2381872_GC"' `"rs4480202_TC"' `"rs4878203_CG"' `"rs4878811_CA"' `"rs769153_GT"' `"rs7859877_GA"' `"rs7865905_CT"' `"rs7866277_CT"' `"rs987187_TC"' `"rs998226_AG"' `"rs10156653_CT"' `"rs1424482_CT"' `"rs1582620_CT"' `"rs4237255_GT"' `"rs4744675_GA"' `"rs8187876_CT"' `"rs8187890_TC"' `"rs8187950_AG"' `"rs12910509_CT"' `"rs12914603_CG"' `"rs1441817_CG"' `"rs1550573_CT"' `"rs1899355_TC"' `"rs2067062_TG"' `"rs2642630_CG"' `"rs2642632_CA"' `"rs2642637_AG"' `"rs2899613_TC"' `"rs3784264_TG"' `"rs4238327_CT"' `"rs4646557_TC"' `"rs4646572_TC"' `"rs7169439_GA"' `"rs8041644_TC"' `"rs8041922_CA"' `"rs12441450_AG"' `"rs14226_GA"' `"rs3803433_GA"' `"rs4965726_CT"' `"rs4965741_CT"' `"rs7497343_AG"'
matrix rownames G= "blank" `"rs1042026_TC"' `"rs1154457_CT"' `"rs1159918_AC"' `"rs1229982_TG"' `"rs1229984_TC"' `"rs12505135_AG"' `"rs1662037_GA"' `"rs166892_CT"' `"rs1693439_AG"' `"rs17526590_GA"' `"rs17583753_GA"' `"rs2009181_AG"' `"rs2173201_CA"' `"rs2654849_GT"' `"rs284792_CT"' `"rs698_TC"' `"rs931635_AG"' `"rs994772_CT"' `"rs10814685_AG"' `"rs10973759_AG"' `"rs10973794_GA"' `"rs10973813_GA"' `"rs10973819_GA"' `"rs10973827_CT"' `"rs11792732_TC"' `"rs12336048_GA"' `"rs2181139_CT"' `"rs2381872_GC"' `"rs4480202_TC"' `"rs4878203_CG"' `"rs4878811_CA"' `"rs769153_GT"' `"rs7859877_GA"' `"rs7865905_CT"' `"rs7866277_CT"' `"rs987187_TC"' `"rs998226_AG"' `"rs10156653_CT"' `"rs1424482_CT"' `"rs1582620_CT"' `"rs4237255_GT"' `"rs4744675_GA"' `"rs8187876_CT"' `"rs8187890_TC"' `"rs8187950_AG"' `"rs12910509_CT"' `"rs12914603_CG"' `"rs1441817_CG"' `"rs1550573_CT"' `"rs1899355_TC"' `"rs2067062_TG"' `"rs2642630_CG"' `"rs2642632_CA"' `"rs2642637_AG"' `"rs2899613_TC"' `"rs3784264_TG"' `"rs4238327_CT"' `"rs4646557_TC"' `"rs4646572_TC"' `"rs7169439_GA"' `"rs8041644_TC"' `"rs8041922_CA"' `"rs12441450_AG"' `"rs14226_GA"' `"rs3803433_GA"' `"rs4965726_CT"' `"rs4965741_CT"' `"rs7497343_AG"'
matrix rownames H= "blank" `"rs1042026_TC"' `"rs1154457_CT"' `"rs1159918_AC"' `"rs1229982_TG"' `"rs1229984_TC"' `"rs12505135_AG"' `"rs1662037_GA"' `"rs166892_CT"' `"rs1693439_AG"' `"rs17526590_GA"' `"rs17583753_GA"' `"rs2009181_AG"' `"rs2173201_CA"' `"rs2654849_GT"' `"rs284792_CT"' `"rs698_TC"' `"rs931635_AG"' `"rs994772_CT"' `"rs10814685_AG"' `"rs10973759_AG"' `"rs10973794_GA"' `"rs10973813_GA"' `"rs10973819_GA"' `"rs10973827_CT"' `"rs11792732_TC"' `"rs12336048_GA"' `"rs2181139_CT"' `"rs2381872_GC"' `"rs4480202_TC"' `"rs4878203_CG"' `"rs4878811_CA"' `"rs769153_GT"' `"rs7859877_GA"' `"rs7865905_CT"' `"rs7866277_CT"' `"rs987187_TC"' `"rs998226_AG"' `"rs10156653_CT"' `"rs1424482_CT"' `"rs1582620_CT"' `"rs4237255_GT"' `"rs4744675_GA"' `"rs8187876_CT"' `"rs8187890_TC"' `"rs8187950_AG"' `"rs12910509_CT"' `"rs12914603_CG"' `"rs1441817_CG"' `"rs1550573_CT"' `"rs1899355_TC"' `"rs2067062_TG"' `"rs2642630_CG"' `"rs2642632_CA"' `"rs2642637_AG"' `"rs2899613_TC"' `"rs3784264_TG"' `"rs4238327_CT"' `"rs4646557_TC"' `"rs4646572_TC"' `"rs7169439_GA"' `"rs8041644_TC"' `"rs8041922_CA"' `"rs12441450_AG"' `"rs14226_GA"' `"rs3803433_GA"' `"rs4965726_CT"' `"rs4965741_CT"' `"rs7497343_AG"'
matrix rownames I= "blank" `"rs1042026_TC"' `"rs1154457_CT"' `"rs1159918_AC"' `"rs1229982_TG"' `"rs1229984_TC"' `"rs12505135_AG"' `"rs1662037_GA"' `"rs166892_CT"' `"rs1693439_AG"' `"rs17526590_GA"' `"rs17583753_GA"' `"rs2009181_AG"' `"rs2173201_CA"' `"rs2654849_GT"' `"rs284792_CT"' `"rs698_TC"' `"rs931635_AG"' `"rs994772_CT"' `"rs10814685_AG"' `"rs10973759_AG"' `"rs10973794_GA"' `"rs10973813_GA"' `"rs10973819_GA"' `"rs10973827_CT"' `"rs11792732_TC"' `"rs12336048_GA"' `"rs2181139_CT"' `"rs2381872_GC"' `"rs4480202_TC"' `"rs4878203_CG"' `"rs4878811_CA"' `"rs769153_GT"' `"rs7859877_GA"' `"rs7865905_CT"' `"rs7866277_CT"' `"rs987187_TC"' `"rs998226_AG"' `"rs10156653_CT"' `"rs1424482_CT"' `"rs1582620_CT"' `"rs4237255_GT"' `"rs4744675_GA"' `"rs8187876_CT"' `"rs8187890_TC"' `"rs8187950_AG"' `"rs12910509_CT"' `"rs12914603_CG"' `"rs1441817_CG"' `"rs1550573_CT"' `"rs1899355_TC"' `"rs2067062_TG"' `"rs2642630_CG"' `"rs2642632_CA"' `"rs2642637_AG"' `"rs2899613_TC"' `"rs3784264_TG"' `"rs4238327_CT"' `"rs4646557_TC"' `"rs4646572_TC"' `"rs7169439_GA"' `"rs8041644_TC"' `"rs8041922_CA"' `"rs12441450_AG"' `"rs14226_GA"' `"rs3803433_GA"' `"rs4965726_CT"' `"rs4965741_CT"' `"rs7497343_AG"'
matrix rownames J= "blank" `"rs1042026_TC"' `"rs1154457_CT"' `"rs1159918_AC"' `"rs1229982_TG"' `"rs1229984_TC"' `"rs12505135_AG"' `"rs1662037_GA"' `"rs166892_CT"' `"rs1693439_AG"' `"rs17526590_GA"' `"rs17583753_GA"' `"rs2009181_AG"' `"rs2173201_CA"' `"rs2654849_GT"' `"rs284792_CT"' `"rs698_TC"' `"rs931635_AG"' `"rs994772_CT"' `"rs10814685_AG"' `"rs10973759_AG"' `"rs10973794_GA"' `"rs10973813_GA"' `"rs10973819_GA"' `"rs10973827_CT"' `"rs11792732_TC"' `"rs12336048_GA"' `"rs2181139_CT"' `"rs2381872_GC"' `"rs4480202_TC"' `"rs4878203_CG"' `"rs4878811_CA"' `"rs769153_GT"' `"rs7859877_GA"' `"rs7865905_CT"' `"rs7866277_CT"' `"rs987187_TC"' `"rs998226_AG"' `"rs10156653_CT"' `"rs1424482_CT"' `"rs1582620_CT"' `"rs4237255_GT"' `"rs4744675_GA"' `"rs8187876_CT"' `"rs8187890_TC"' `"rs8187950_AG"' `"rs12910509_CT"' `"rs12914603_CG"' `"rs1441817_CG"' `"rs1550573_CT"' `"rs1899355_TC"' `"rs2067062_TG"' `"rs2642630_CG"' `"rs2642632_CA"' `"rs2642637_AG"' `"rs2899613_TC"' `"rs3784264_TG"' `"rs4238327_CT"' `"rs4646557_TC"' `"rs4646572_TC"' `"rs7169439_GA"' `"rs8041644_TC"' `"rs8041922_CA"' `"rs12441450_AG"' `"rs14226_GA"' `"rs3803433_GA"' `"rs4965726_CT"' `"rs4965741_CT"' `"rs7497343_AG"'

foreach i in A F B G C H {
	di "MATRIX ******`i'"
	matrix li `i'
	svmat `i'
	}
	
