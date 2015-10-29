data {
	int N;			//total number of observations
	int n;			//number of subjects
	int m;			//number of loci
	int ii[N];		//subject number
	int jj[N];		//locus number
	int X[N];		//number of minor alleles

	real v0;		//prior scale for v
}

parameters {
	real<lower=0,upper=1> f_mu;	//mean inbreeding
	real<lower=0,upper=pi()/2> f_v_raw;		//f inverse "sample size"
	real<lower=0,upper=1> maf[m];	//minor allele frequencies
	real<lower=0,upper=1> f[n];	//inbreeding coeficient
}

transformed parameters{
	real<lower=0> f_v;
	real<lower=0> a;
	real<lower=0> b;
	
	f_v <- v0 * tan(f_v_raw);

	a <- f_mu / f_v;
	b <- (1-f_mu) / f_v;
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
	
	f ~ beta(a,b);
	//f_v ~ cauchy(0,v0);
	//f_mu ~ uniform(0,1);
	
}
