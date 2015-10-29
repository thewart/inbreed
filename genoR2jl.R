source('~/Dropbox/monkeybris/rscript/setup.R')

X <- as.matrix(fv$X)
#X[is.na(X)] <- -1

if (!exists("Dthresh"))
  Dthresh <- 3e5

shitlist <- c()
for (i in 1:length(unique(fv$SNP$chrom)))
{
  ic <- unique(fv$SNP$chrom)[i]
  loci <- which(fv$SNP$chrom == ic)
  loc <- fv$SNP$loc[loci]
  nl <- length(loc)
  
  D <- abs(kronecker(t(rep(1,times=nl)),loc) - kronecker(rep(1,times=nl),t(loc)))
  foo <- which(D<Dthresh & D>0,arr=T)
  dump <- vector("numeric")
  
  moo <- foo[!(foo[,1] %in% dump | foo[,2] %in% dump),]
  while (length(moo) > 0)
  {
    dump <- c(dump, as.numeric(names(which.max(table(moo[,1])))))
    moo <- foo[!(foo[,1] %in% dump | foo[,2] %in% dump),]
  }
  
  shitlist <- c(shitlist,loci[dump])
}
X <- X[,-shitlist]

source('~/Dropbox/monkeybris/rscript/pedigree.preproc.batch.R')
reped <- ped.matchnames(as.character(fv$MON$Animal.ID),pedigree$id)
pedigree <- ped.replace(pedigree,reped$oldID,reped$ID)
redped <- ped.trace(fv$MON$Animal.ID,pedigree)
A <-with(redped,kinship2::kinship(id,dadid=sire,momid=dam))
A <- A[rownames(A) %in% rownames(X),rownames(A) %in% rownames(X)]
fp = diag(A)[match(rownames(X),rownames(A))] - 0.5
