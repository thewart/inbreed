data {
	int N;			//total number of observations
	int n;			//number of subjects
	int m;			//number of loci
	int ii[N];		//subject number
	int jj[N];		//locus number
	int X[N];		//number of minor alleles

	real mu0;		//prior scale for mu
	real s0;		//prior scale for v
}

parameters {
	real f_mu;			//mean inbreeding
	real<lower=0> f_s;		//inbreeding scale
	real<lower=0,upper=1> maf[m];	//minor allele frequencies
	real<lower=0,upper=1> f_raw[n];	//inbreeding coeficient
}

transformed parameters{
	real<lower=0> f[n];
	
	for (i in 1:n) {
		f[i] <- inv_logit(f_mu + f_raw[i]*f_s); 
	}
	
}

model {

	for (i in 1:N) {
		real pf0;
		real pf1;

		pf0 <- log1m(f[ii[i]]) + binomial_log(X[i],2,maf[jj[i]]);
		pf1 <- log(f[ii[i]]) + 
			if_else(X[i]!=1,binomial_log(X[i]/2,1,maf[jj[i]]),log(0));
		increment_log_prob(log_sum_exp(pf0,pf1));
	}
	
	f_raw ~ normal(0,1);
	f_s ~ cauchy(0,s0);
	f_mu ~ normal(0,mu0);
	
}

