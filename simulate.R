## qR.rb -j sim -r -p skylake-himem -t "01:00:00" -c 10 -y 2-999 simulate.R
# j=A,B,C,D,E for imputed data sets
# i for the random permutation

library(data.table)

library(magrittr)

task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if(task_id_string=="")
  task_id_string="833"
task_id <- as.numeric(task_id_string)
i <- as.numeric(task_id_string)
outfile=paste0("data/output-final/out",i,".csv.gz")
if(file.exists(outfile))
  stop("cowardly refusing to overwrite output file ",outfile)

set.seed(i)

pheno=lapply(toupper(letters[1:5]), function(x)
  fread(file.path("data",paste0("merged",x,".fam")), header = F))
pheno=do.call("cbind", lapply(pheno, function(x) x$V6))
message("original data")
print(head(pheno))

## add latest ACR phenotype
merged_fam <- read.table("~/share2/People/ANNA/T1D-GWAS/1-GenotypeQCImp/Merged/merged.fam", header = F)
matchIDs <- readRDS("~/share2/People/ANNA/T1D-GWAS/1-GenotypeQCImp/sample-keys/genotypedsamples-IDs.RDS")
merged_fam$altID <- matchIDs$pheno_ID[match(merged_fam$V2, matchIDs$geno_ID)]
phenotypes.acr <- readRDS("~/share2/People/ANNA/T1D-GWAS/2-PhenotypeDataPrep/c-MakePhenotypes/Secondary-phenotypes/secondary-phenotypes.RDS")
phenotypes.acr = phenotypes.acr[match(merged_fam$altID, phenotypes.acr$StudyID),]
pheno=cbind(pheno, phenotypes.acr[,grep("latestACR",names(phenotypes.acr))])

## add average residual ACR phenotype
pheno=cbind(pheno, phenotypes.acr[,grep("averes",names(phenotypes.acr))])

## permute
## NB - some missing values - these need to stay put
idx=which(!is.na(pheno[,1]))
sidx=sample(idx)
pheno[idx,]=pheno[sidx,]
message("permuted data")
print(head(pheno))

## add a normal random deviate
x=numeric(nrow(pheno))
x[idx]=rnorm(length(idx))
x[-idx]=NA
pheno=cbind(pheno, rand=x)

file_pheno=file.path("data",paste0("merged-tmp",i,".pheno"))
write.table(pheno,
            file=file_pheno,
            row.names = F, col.names = F, quote = F, sep="\t")

npheno=ncol(pheno)

library(parallel); options(mc.cores=3)
LETTERS=toupper(letters[1:npheno])
mclapply(1:npheno, function(j) {
  f=paste0("data/output-new/out",i,"-",LETTERS[j],".assoc.txt")
  fzip=paste0("data/output-new/out",i,"-",LETTERS[j],".assoc.txt.gz")
  if(!file.exists(f) && !file.exists(fzip))
    system(
      paste0("data-orig/../../3-GWAS/gemma-0.98.5-linux-static-AMD64 ",
             "-bfile data/mergedA -lm 1 ",
             "-p data/merged-tmp",i,".pheno -n ",j,
             " -outdir data/output-new -o out",i,"-",LETTERS[j])
    )
})

list.files("data/output-new",pattern=paste0("out",i,"-.*.assoc.txt"))
message("gemma complete")
## if(!interactive())
##   q("no")

## read in GWAS results
gwas_out=mclapply(LETTERS, function(j) {
  f=paste0("data/output-new/out",i,"-",j,".assoc.txt")
  fzip=paste0("data/output-new/out",i,"-",j,".assoc.txt.gz")
  if(file.exists(fzip)) {
    fread(cmd=paste("zcat",fzip), select = c(9:10))
  } else if(file.exists(f)) {
    fread(f, select = c(9:10))
    system(paste("gzip",f))
  } else {
    NULL
  }
})

message("read ",length(gwas_out)," copies of gwas_out")
print(str(gwas_out))

beta=lapply(gwas_out, function(x) x$beta) %>% do.call("cbind", .)
varbeta=lapply(gwas_out, function(x) x$se^2) %>% do.call("cbind", .)

# Pool using Rubin's rules

pool=function(beta, varbeta) {
library(matrixStats)
  pooledMean <- rowMeans(beta)
  betweenVar <- rowMeans(varbeta) # mean of variances - think this should be called within
withinVar <- rowVars(beta) # variance of variances - think this should be called between
m=ncol(beta)
  dfCorrection <- (m+1)/m # dfCorrection

  totVar <- betweenVar + withinVar*dfCorrection # total variance
  pooledSE <- sqrt(totVar) # standard error

  lambda <- (withinVar + (withinVar/ncol(beta)))/totVar # fraction of missing information
  n <- sum(!is.na(pheno[,1]))
  k <- 2
  nu_old <- (m-1)/lambda^2
  nu_com <- n-k
  nu_obs <- (nu_com+1)/(nu_com+3)*nu_com*(1-lambda)
  nu_BR <- (nu_old*nu_obs)/(nu_old+nu_obs)
summary(nu_BR)

pooledP <- pt(q = abs(pooledMean / pooledSE), df = nu_BR, lower.tail = FALSE) * 2
}

pool.traj=pool(beta[,1:5], varbeta[,1:5])
pool.latest=pool(beta[,6:10], varbeta[,6:10])
pool.averes=pool(beta[,11:15], varbeta[,11:15])
p.norm=pnorm(-abs(beta[,16])/sqrt(varbeta[,16])) * 2

result=data.table(pool.traj=pool.traj,pool.latest=pool.latest,pool.averes=pool.averes, p.norm=p.norm)
fwrite(result, file = outfile)

message("saved to ",outfile)
unlink(file_pheno)
