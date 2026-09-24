options(stringsAsFactors=FALSE)
library(ggplot2)
args_full <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_full, value = TRUE)
if (!length(script_arg)) stop("Run with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
root <- normalizePath(Sys.getenv("SEER_PROJECT_DIR", unset = dirname(dirname(script_path))), winslash = "/", mustWork = TRUE)

out <- Sys.getenv("SEER_REVISION_OUTPUT_DIR", unset = file.path(root, "04_results", "revision_20260921"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
tab <- function(n) read.csv(file.path(root,'04_results/tables',n),check.names=FALSE)
write <- function(d,n) write.csv(d,file.path(out,n),row.names=FALSE)
s <- tab('table_28_age_by_first_cancer_site.csv')
mapping <- data.frame(Original=sort(unique(s$First_Cancer_Site_Group)))
mapping$Shared <- ifelse(mapping$Original %in% c('Myeloma','Mesothelioma','Kaposi Sarcoma','Miscellaneous'),'Other first cancers',mapping$Original)
write(mapping,'classification_map.csv')
s$Shared <- mapping$Shared[match(s$First_Cancer_Site_Group,mapping$Original)]
ag <- function(d,cols) aggregate(d[c('Observed','Expected','Persons','Person Years at Risk')],d[cols],sum)
g <- ag(s,c('Age_Group','Shared'))
g$Age_Group <- relevel(factor(g$Age_Group),'30-39')
stopifnot(sum(g$Persons)==276844,all(g$Expected>0),sum(g$Observed)==345)
fit <- function(d) glm(Observed~Shared+Age_Group,offset=log(Expected),family=poisson,data=d)
m <- fit(g); m0 <- glm(Observed~Shared,offset=log(Expected),family=poisson,data=g)
effects <- function(m) { k <- grep('Age_Group',names(coef(m))); b <- coef(m)[k]; se <- sqrt(diag(vcov(m)))[k]; data.frame(Age_Group=sub('Age_Group','',names(b)),Ratio=exp(b),Lower=exp(b-1.96*se),Upper=exp(b+1.96*se)) }
e <- effects(m); write(e,'pooled_effects.csv'); write(g,'shared_age_cells.csv')
tests <- data.frame(Test='Age adjusted for 15 shared groups',P=anova(m0,m,test='Chisq')[2,'Pr(>Chi)'])
l <- tab('table_04_age_by_latency.csv')
null <- glm(Observed~factor(Age_Group)+factor(Latency),offset=log(Expected),family=poisson,data=l)
# The full age by latency model has one parameter per cell, hence is saturated.
observed_dev <- deviance(null)
set.seed(20260921); B <- 1999L
boot <- replicate(B,{dd<-l;dd$Observed<-rpois(nrow(dd),fitted(null)); suppressWarnings(deviance(glm(Observed~factor(Age_Group)+factor(Latency),offset=log(Expected),family=poisson,data=dd)))})
stopifnot(all(is.finite(boot)))
bp <- (1+sum(boot>=observed_dev))/(B+1)
tests <- rbind(tests,data.frame(Test='Age by latency parametric bootstrap 1999 replicates',P=bp))
# Rounding sensitivity uses raw fine-age cells; no claim of recovering unrounded E.
raw <- read.delim(file.path(root,'02_raw_exports/05_first_cancer_diag_age_latency_mpsir.txt'),check.names=FALSE)
raw <- raw[raw[['Selected Events']]=='Brain and Other Nervous System' & raw$Latency=='Total' & raw[[4]] %in% mapping$Original & raw$Persons>0,]
a <- as.numeric(sub('^([0-9]+).*','\\1',raw[[3]]))
raw$Age_Group <- cut(a,c(-1,14,19,29,39),labels=c('0-14','15-19','20-29','30-39'))
raw$Shared <- mapping$Shared[match(raw[[4]],mapping$Original)]
set.seed(20260922)
pert <- replicate(200,{d<-raw;d$Expected<-runif(nrow(d),pmax(0,d$Expected-0.005),d$Expected+0.005);h<-ag(d,c('Age_Group','Shared'));h$Age_Group<-relevel(factor(h$Age_Group),'30-39');effects(fit(h))$Ratio})
write(data.frame(Age_Group=e$Age_Group,Min=apply(pert,1,min),Max=apply(pert,1,max)),'rounding_sensitivity.csv')
write(tests,'new_tests.csv')
s$Age_Band <- ifelse(s$Age_Group=='0-14','0-14 years','15-39 years')
h <- ag(s,c('Age_Band','Shared'))
h$SIR <- h$Observed/h$Expected
h$Lower <- ifelse(h$Observed==0,0,qchisq(.025,2*h$Observed)/2/h$Expected)
h$Upper <- qchisq(.975,2*(h$Observed+1))/2/h$Expected
h$EAR <- (h$Observed-h$Expected)/h[['Person Years at Risk']]*100000
write(h,'shared_figure_data.csv')
# Contract: shared descriptive categories permit visual comparison; no causal treatment claim.
# Paired quantitative forest panels, exact Poisson intervals, 183 x 160 mm, 600 dpi TIFF.
h$Shared <- factor(h$Shared,levels=rev(sort(unique(h$Shared))))
p <- ggplot(h,aes(y=Shared,x=SIR))+geom_vline(xintercept=1,linetype=2,colour='grey55')+
 geom_segment(aes(x=pmax(.02,Lower),xend=Upper,yend=Shared),linewidth=.4,colour='#346A85')+
 geom_point(data=h[h$Observed>0,],size=1.7,colour='#346A85')+
 geom_text(aes(x=480,label=paste0('n=',Observed)),size=2.5,hjust=1)+
 scale_x_log10(limits=c(.02,550),breaks=c(.1,1,10,100))+facet_grid(.~Age_Band)+
 labs(x='Standardized incidence ratio (95% CI)',y=NULL)+theme_bw(base_size=8,base_family='Arial')+
 theme(panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),strip.background=element_rect(fill='#EDF2F4'),strip.text=element_text(face='bold'),plot.margin=margin(8,10,8,8))
ggsave(file.path(out,'shared_first_cancer.tiff'),p,width=183,height=160,units='mm',dpi=600,compression='lzw')
ggsave(file.path(out,'shared_first_cancer.png'),p,width=183,height=160,units='mm',dpi=160)
ggsave(file.path(out,'shared_first_cancer.svg'),p,width=183,height=160,units='mm',device=svglite::svglite)
print(e);print(tests);print(read.csv(file.path(out,'rounding_sensitivity.csv')))
cat('Pooled model:',nrow(g),'cells;',sum(g$Persons),'persons;',sum(g$Observed),'events; converged',m$converged,'\n')
