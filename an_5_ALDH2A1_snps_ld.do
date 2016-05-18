//Neil Davies 18/05/16
//This runs the joint adjustment procedure
//There are three SNPs in ALDH1A2 which are in LD:
//rs1441817_CG
//rs8041922_CA
//rs12910509_CT

//This scripts run a joint analysis to see which SNP is most strongly associated with survival:



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
	
//	cap{
		keep if studyid=="`k'"

		stset lengthfudied2, id(icogs_sample_id) failure(death2) 

		stcox highgrade rs1441817_CG rs8041922_CA rs12910509_CT pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8
		regsave  using "results/24715survivaladjustedSNP1highgrade_joint_adj_`k'", detail(all) ci pval replace

		stcox lowgrade rs1441817_CG rs8041922_CA rs12910509_CT  pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8
		regsave using "results/24715survivaladjustedSNP1lowgrade_joint_adj_`k'", detail(all) ci pval replace 



//		}
	restore
	}

use "results/24715survivaladjustedSNP1highgrade_joint_adj_EPIC",clear
foreach k in    ESTHER FHCRC MAYO MCCS MEC /// 
	PPF-UNIS ProMPT TAMPERE UKGPCS ULM {
	append using "results/24715survivaladjustedSNP1highgrade_joint_adj_`k'"
	}
