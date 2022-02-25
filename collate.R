library(data.table)
files=list.files("data/output-final",full=TRUE)
message("files found: ",length(files))

get_quantiles=function(f) {
  message("reading ",f)
  x=fread(f)
  c(traj=min(x$pool.traj), latest=min(x$pool.latest), averes=min(x$pool.averes), norm=min(x$p.norm))
}

library(parallel)
if(!interactive()) {
  options(mc.cores=6)
  q=mclapply(files, get_quantiles)
  qmat=do.call(rbind,q)
  save(qmat, file="qmat.RData")
} else {
  (load("qmat.RData"))
}
summary(qmat[,1])
dim(qmat)

null_traj=qmat[,1]
null_latest=qmat[,2]
null_averes=qmat[,3]
null_norm=qmat[,4]
traj.p=8.9e-11
latest.p=6.9e-9
message("empirical traj p=",mean(null_traj <= traj.p))
message("empirical latest p=",mean(null_latest <= latest.p))
message("empirical gwas sig 5e-8 p=",mean(null_norm <= 5e-8))
message("empirical gwas sig 1e-8 p=",mean(null_norm <= 1e-8))


par(mfrow=c(1,3))
qqplot(-log10(null_norm),-log10(null_traj),main="trajectory"); abline(0,1); abline(v=-log10(5e-8)); abline(h=-log10(traj.p))
qqplot(-log10(null_norm),-log10(null_latest),main="latest"); abline(0,1); abline(v=-log10(5e-8)); abline(h=-log10(acr.p))
qqplot(-log10(null_norm),-log10(null_averes),main="averes"); abline(0,1); abline(v=-log10(5e-8)); abline(h=-log10(acr.p))

## quantile normalise
traj=approx(sort(null_traj), sort(null_norm), xout=8.9e-11)
message("traj raw p=",traj$x,"\tadjusted p=",traj$y)
acr=approx(sort(null_latest), sort(null_norm), xout=6e-9)
message("latest raw p=",acr$x,"\tadjusted p=",acr$y)

## full sumstats
traj=fread("~/share2/People/ANNA/T1D-GWAS/Data/GWAS-sumstats/trajectory-GWAS-sumstats.tsv")
latest=fread("~/share2/People/ANNA/T1D-GWAS/Data/GWAS-sumstats/latestACR-GWAS-sumstats.tsv")
averes=fread("~/share2/People/ANNA/T1D-GWAS/Data/GWAS-sumstats/averageACR-GWAS-sumstats.tsv")

## find hits
x=copy(traj)
collect_hits=function(x) {
  x[,is_hit:=FALSE]
  while(any(x$p_value < 1e-6 & x$is_hit==FALSE)) {
    wh=which.max(-log10(x$p_value) * (x$is_hit==FALSE))
    x$is_hit[wh]=TRUE
    print(x[wh,])
    x[,in_region:=chromosome == x$chromosome[wh] & abs(base_pair_location - x$base_pair_location[wh]) < 1e+6]
    x=x[in_region==FALSE | is_hit==TRUE]
  }
  x[is_hit==TRUE]
}

hits_traj=collect_hits(traj)
hits_latest=collect_hits(latest)
hits_averes=collect_hits(averes)

hits_all=rbind(hits_traj[,.(chromosome, base_pair_location,effect_allele, other_allele,hit_traj=is_hit)],
               hits_latest[,.(chromosome, base_pair_location,effect_allele, other_allele,hit_latest=is_hit)],
               hits_averes[,.(chromosome, base_pair_location,effect_allele, other_allele,hit_averes=is_hit)],fill=TRUE)
## collapse any repeats
hits_all=hits_all[,.(hit_traj=any(hit_traj,na.rm=TRUE),
                     hit_latest=any(hit_latest,na.rm=TRUE),
                     hit_averes=any(hit_averes,na.rm=TRUE)),
                  by=c("chromosome","base_pair_location","effect_allele","other_allele")]
## add estimates
hits_all=merge(hits_all,traj[,.(chromosome,base_pair_location,effect_allele,other_allele,
                                beta_traj=beta,se_traj=standard_error,p_traj=p_value)])
hits_all=merge(hits_all,latest[,.(chromosome,base_pair_location,effect_allele,other_allele,
                                beta_latest=beta,se_latest=standard_error,p_latest=p_value)])
hits_all=merge(hits_all,averes[,.(chromosome,base_pair_location,effect_allele,other_allele,
                                beta_averes=beta,se_averes=standard_error,p_averes=p_value)])
## estimate adjusted p values
hits_all$adjp_traj=approx(sort(null_traj), sort(null_norm), xout=hits_all$p_traj)$y
hits_all$adjp_latest=approx(sort(null_latest), sort(null_norm), xout=hits_all$p_latest)$y
hits_all$adjp_averes=approx(sort(null_averes), sort(null_norm), xout=hits_all$p_averes)$y

## final results
hits_all[pmin(adjp_traj,adjp_latest,adjp_averes,na.rm=TRUE) < 1e-6]



