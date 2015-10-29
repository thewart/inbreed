function flipperbaby(geno,burn,thin,iter)

fprior = [.5,5];
pprior = [1,1];

m,n = size(geno)
ft = zeros(Float64,n);
pt = Array(Float64,m);
for j in 1:m pt[j] = mean(geno[j,geno[j,:].!=-1])/2; end
z = Array(Bool,(m,n));

nn = sum(geno.!=-1,1)';

nz = zeros(Int64,n);
nma = zeros(Int64,m);
na = zeros(Int64,m);

#initialize gibbs chain
saveiter = (burn+1):thin:iter;
nsave = length(saveiter);
niter = maximum(saveiter);
f = Array(Float64,(n,nsave));
p = Array(Float64,(m,nsave));

for t in 1:iter

  for i in 1:n
    for j in 1:m

      if geno[j,i] == -1                       #skip if missing data
        continue
      elseif geno[j,i] == 1                    #can't be IBD if heterozygous
        z[j,i] = false;
      else
        pz1 = pt[j]^(geno[j,i]==2) * (1-pt[j])^(geno[j,i]==0) * f[i];
        pz0 = pt[j]^geno[j,i] * (1-pt[j])^(2-geno[j,i]) * (1-f[i]);
        z[j,i] = rand() < (pz1/(pz1+pz0));
      end

      nz[i] += z[j,i];                         #count IBD alleles
      na[j] += z[j,i] ? 1 : 2;                 #count only one total allele if IBD
      nma[j] += geno[j,i] / (z[j,i] ? 2 : 1);  #count only one minor allele if IBD

    end
  end

  for i in 1:n
    ft[i] = rand( Beta(fprior[1]+nz[i], fprior[2]+nn[i]-nz[i]) )
    nz[i] = 0;
  end

  for j in 1:m
    pt[j] = rand( Beta(pprior[1]+nma[j], pprior[2]+na[j]-nma[j]) );
    na[j] = 0;
    nma[j] = 0;
  end

  nsamp = findin(saveiter,t);
  if !isempty(nsamp)
    f[:,nsamp] = ft;
    p[:,nsamp] = pt;
  end

end

return f,p

end
