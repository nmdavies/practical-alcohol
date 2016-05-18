//Neil Davies 15/03/16
//This reruns Clair's meta-regressions using the file she sent

//Start log
cap: log close
log using "/Volumes/filestore-2/Post doc/Alcohol Practical Brunner/meta_reg_log.txt", text replace

//Import her spreadsheet
import excel "/Volumes/filestore-2/Post doc/Alcohol Practical Brunner/Documents/sig results for metareg_nmd_160315.xlsx", ///
	sheet("Sheet1") cellrange(A10:S63) firstrow clear

	
//Clean
drop if var==""|var=="var"
destring *, replace

//Drop excluded studies
drop if substr(Q,1,3)=="exc"

//Exclude rs4965726_CT
drop if var=="rs4965726_CT"

//Setup results file
metareg coef PSA if stderr !=0 & var=="rs1441817_CG",   wsse(stderr) eform
regsave using "/Volumes/filestore-2/Post doc/Alcohol Practical Brunner/meta_reg_results", detail(all) pval ci replace


foreach i in rs1441817_CG rs8041922_CA rs12910509_CT rs10973794_GA{
	foreach j in PSA FHx USA Age{
		di "************************"
		di "`i'"
		di "************************"
		metareg coef `j'  if stderr !=0 & var=="`i'",   wsse(stderr) eform
		regsave using "/Volumes/filestore-2/Post doc/Alcohol Practical Brunner/meta_reg_results", detail(all) pval ci append
		}
	}

use "/Volumes/filestore-2/Post doc/Alcohol Practical Brunner/meta_reg_results",clear

drop if var=="_cons"
format coef ci_* pval stderr %9.2f
order depvar var coef ci_* pval

replace coef=exp(coef)
replace ci_lower=exp(ci_lower)
replace ci_upper=exp( ci_upper)

//Sort in right order
gen n=_n
gsort - n
	
log close



	
