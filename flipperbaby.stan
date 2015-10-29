data {
	int N;			//total number of observations
	int n;			//number of subjects
	int m;			//number of loci
	int ii[N];		//subject number
	int jj[N];		//locus number
	int X[N];		//number of minor alleles
}

parameters {
	real<lower=0,upper=1> maf[m];	//minor allele frequencies
	real<lower=0,upper=1> f[n];	//inbreeding coeficient
}

model	{

	for (i in 1:N) {
		real pf0;
		real pf1;

		pf0 <- log(1-f[ii[i]]) + binomial_log(X[i],2,maf[jj[i]]);
		pf1 <- log(f[ii[i]]) + 
			if_else(X[i]!=1,binomial_log(X[i]/2,1,maf[jj[i]]),log(0));
		increment_log_prob(log_sum_exp(pf0,pf1));
	}
	
}
