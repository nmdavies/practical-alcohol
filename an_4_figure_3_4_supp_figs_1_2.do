//Neil Davies 18/05/16
//This produces Supp Figures 1 and 2:



use "results/24715survivaladjustedSNP1caco_allstudies" , clear

rename var snp

drop if studyid=="ProMPT"|studyid=="ESTHER"|studyid=="CPCS2"|studyid=="PPF-UNIS"

local n="rs12910509_CT"
metan coef ci_lower ci_upper if snp=="`n'", fixed eform lcols(studyid) title(`n') ///
		sub("Fixed effect association with survival in any prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k')  ///
		xlabel(0.5 ,0.75, 1, 1.25 ,1.50)
graph export "results/supp_figure_1_fixed metan plot of `n' with survival any pca 24715.pdf" , as(eps) replace


local n="rs8041922_CA"
metan coef ci_lower ci_upper if snp=="`n'", fixed eform lcols(studyid) title(`n') ///
		sub("Fixed effect association with survival in any prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k') 
graph export "results/supp_figure_1_fixed metan plot of `n' with survival any pca 24715.pdf" , as(eps) replace

/*
Figure 3 – Meta-analysis of prostate cancer specific survival following a diagnosis of any prostate cancer, in association with rs1441817 in ALDH1A2.

Figure 4 – Meta-analysis of prostate cancer specific survival following a diagnosis of low grade prostate cancer, in association with rs10973794 in ALDH1B1.
*/

use "results/24715survivaladjustedSNP1caco_allstudies" , clear

drop if studyid=="ProMPT"|studyid=="CPCS2"|studyid=="PPF-UNIS"

rename var snp
local n="rs1441817_CG"
metan coef ci_lower ci_upper if snp=="`n'", fixed eform lcols(studyid)  ///
		scheme(s1mono) xtitle("Odds ratio") ylabel(`k') ///
		xlabel(0.005,0.05,0.5,1,1.5,2.5) force
graph export "results/figure_3_fixed metan plot of rs1441817_CG with survival any pca.eps" , as(eps) replace

use "results/24715survivaladjustedSNP1lowgrade_allstudies" , clear
rename var snp

drop if studyid=="ProMPT"|studyid=="CPCS1"|studyid=="CPCS2"|studyid=="PPF-UNIS"||studyid=="MEC"		

local n="rs10973794_GA"
metan coef ci_lower ci_upper if snp=="`n'", fixed eform lcols(studyid)  ///
		scheme(s1mono) xtitle("Odds ratio") ylabel(`k') ///
		xlabel(0.005,0.05,0.5,1,1.5,2.5) force
graph export "results/figure_4_fixed metan plot of rs10973794_GA with survival any pca.eps" , as(eps) replace
				
		metan coef ci_lower ci_upper if snp=="`n'", random eform lcols(studyid) title("`n'") sub("Random effect association with survival in any prostate cancer") scheme(s1mono) xtitle("Odds ratio") xlabel(0.005,0.05,0.5,1,1.5,2.5) force ytitle("Study ID") ylabel(`k')
			


local n="rs10973794_GA"
metan coef ci_lower ci_upper if snp=="`n'", fixed eform lcols(studyid) title(`n') ///
		sub("Fixed effect association with survival in any prostate cancer") scheme(s1mono) xtitle("Odds ratio") ytitle("Study ID") ylabel(`k') 
graph export "results/supp_figure_1_fixed metan plot of `n' with survival any pca 24715.pdf" , as(eps) replace
