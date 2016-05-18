//Neil Davies 18/05/16
//Updated to improve legibility of Figure 1 and 2

// do file for gen combined survival scatter plots
//created 23/9/14 modified 7/11/14 to include both metans

cd "/Volumes/ecnmd/Post doc/Alcohol Practical Brunner"

//convert adj results file so ready to merge
//pca caco adjusted
use "results/24715survivaladjustedSNP1caco_allstudies", clear
rename var snp
foreach i in _ A C G T{
	replace snp=regexr(snp, "`i'","")
	}
joinby snp using "Raw data/snp_chr_pos.dta", unmatched(both)

destring chr pos , replace

replace chr=4 if snp=="rs166892"
replace pos=100501788 if snp=="rs166892"
replace pos=38398525 if snp=="rs4878811"
replace chr=9  if snp=="rs4878811"

gen logpval=log10(1/pval)
save "results/survivalcacoadjusted_snpinfo", replace


//pca highgrade adj
use "results/24715survivaladjustedSNP1highgrade_allstudies.dta", clear

rename var snp
joinby snp using "Raw data/snp_chr_pos.dta", unmatched(both)

destring chr pos , replace

replace chr=4 if snp=="rs166892"
replace pos=100501788 if snp=="rs166892"
replace pos=38398525 if snp=="rs4878811"
replace chr=9  if snp=="rs4878811"


gen logpval=log10(1/pval)

save "results/survivalhighgradeadjuste_snpinfo" , replace

//pca lowgrade adj
use "results/24715survivaladjustedSNP1lowgrade_allstudies.dta", clear

rename var snp
joinby snp using "Raw data/snp_chr_pos.dta", unmatched(both)

destring chr pos , replace

replace chr=4 if snp=="rs166892"
replace pos=100501788 if snp=="rs166892"
replace pos=38398525 if snp=="rs4878811"
replace chr=9  if snp=="rs4878811"

gen logpval=log10(1/pval)

save "results/survivallowgradeadjuste_snpinfo" , replace

use "results/survivalcacoadjusted_snpinfo",  clear
drop _merge
joinby snp using "working data/metan_results_survival",  unmatched(both)

//First all prostate cancer

gen logpval_fix=log10(A4)*-1
gen logpval_ran=log10(B4)*-1

scatter logpval_fix logpval_ran pos if chr==4 , ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 4 snps with survival in any pca pool fix and ran.eps" , as(eps) replace

scatter logpval_fix logpval_ran pos if  chr==9 & pos<=40000000, ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 9 in gene1 snps with survival in any pca pool fix and ran.eps" , as(eps) replace

scatter logpval_fix logpval_ran pos if  chr==9 & pos>=40000000, ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 9 in gene2 snps with survival in any pca pool fix and ran.eps" , as(eps) replace

scatter logpval_fix logpval_ran pos if  chr==15 & pos<=60000000, ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 15 in gene1 snps with survival in any pca pool fix and ran.eps" , as(eps) replace

scatter logpval_fix logpval_ran pos if  chr==15 & pos>=60000000, ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 15 in gene2 snps with survival in any pca pool fix and ran.eps" , as(eps) replace

//Second high grade disease

replace logpval_fix=log10(C4)*-1
replace logpval_ran=log10(D4)*-1

scatter logpval_fix logpval_ran pos if chr==4 , ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 4 snps with survival in high grade pca pool fix and ran.eps" , as(eps) replace

scatter logpval_fix logpval_ran pos if  chr==9 & pos<=40000000, ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 9 in gene1 snps with survival in high grade pca pool fix and ran.eps" , as(eps) replace

scatter logpval_fix logpval_ran pos if  chr==9 & pos>=40000000, ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 9 in gene2 snps with survival in high grade pca pool fix and ran.eps" , as(eps) replace

scatter logpval_fix logpval_ran pos if  chr==15 & pos<=60000000, ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 15 in gene1 snps with survival in high grade pca pool fix and ran.eps" , as(eps) replace

scatter logpval_fix logpval_ran pos if  chr==15 & pos>=60000000, ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 15 in gene2 snps with survival in high grade pca pool fix and ran.eps" , as(eps) replace
	
//Second low grade disease

replace logpval_fix=log10(E4)*-1
replace logpval_ran=log10(F4)*-1

scatter logpval_fix logpval_ran pos if chr==4 , ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 4 snps with survival in low grade pca pool fix and ran.eps" , as(eps) replace

scatter logpval_fix logpval_ran pos if  chr==9 & pos<=40000000, ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 9 in gene1 snps with survival in low grade pca pool fix and ran.eps" , as(eps) replace

scatter logpval_fix logpval_ran pos if  chr==9 & pos>=40000000, ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 9 in gene2 snps with survival in low grade pca pool fix and ran.eps" , as(eps) replace

scatter logpval_fix logpval_ran pos if  chr==15 & pos<=60000000, ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 15 in gene1 snps with survival in low grade pca pool fix and ran.eps" , as(eps) replace

scatter logpval_fix logpval_ran pos if  chr==15 & pos>=60000000, ///
	scheme(s1mono) yscale(range(0 4)) ylabel(0(1)4) yline(1.30103, lpattern(dash) lcolor(black)) yline(2.3802, lpattern(dash) lcolor(gs11)) ///
	symbol(T D) msize(small small small) mfcolor(gs7 gs13 gs10) mlcolor(black black black) mlwidth(thin thin thin) xtitle("") legend(off)
graph export "results/51015manhattan plot of chr 15 in gene2 snps with survival in low grade pca pool fix and ran.eps" , as(eps) replace
