### mutant lines
#this may not be necessary
rm(list = ls(all = TRUE)) 
.libPaths("C:/Program Files/R/R-3.5.0/library")
library(data.table)
library(behavr)
library(scopr)
library(sleepr)
library(ggetho)
library(plotly)
library(survival)  
library(zeitgebr)
library(bit64)
library(stringr)
library(readr)
library(EnvStats)
source("Starting_material.R")

# This is the placement of the data in this computer
DATA_DIR <- "/mnt/ethoscope_results"
#This is the placement of the real data - this is accessible only through imperial VPN
REMOTE_DATA_SOURCE <- "ftp://turing.lab.gilest.ro/auto_generated_data/ethoscope_results/"
#this is the place were a cache version of the data is store. This make the loading faster on the second time.
#CACHE = "etho_cache"

#This is the path where your query/METADATA of the experiment in txt or csv format is saved
METADATA = "mutant SD metadata_new_2.csv"
fread(METADATA)


#if you have the data locally and you are disconected form the internet you should run this instead:
query <- link_ethoscope_metadata(fread(METADATA),
                                 result_dir = DATA_DIR)


#This is the magic step, it loads the data to R and applies a function at the same time, in this case, the sleep annotation.
dt <- load_ethoscope(query,
                     reference_hour = 9.0, 
                     FUN = sleep_annotation)
#cache = CACHE)

dt<-dt[xmv(status)=="ok"]

dt[,t:=t+days(xmv(baseline_days))]

#dt <- dt[,curateDeadAnimals(.SD, max_immobile_live=hours(6)), by=id]

dt<- curate_dead_animals(dt,
                         moving_var = moving,
                         time_window = hours(4),
                         prop_immobile = 0.01,
                         resolution = 3)
#remove dead animals
valid_animals <- dt[,
                    .(is_valid = max(t) > days(0)),
                    by=key(dt)]
unique(valid_animals)[,.N,by="is_valid"]
valid_animals <- valid_animals[is_valid == T]
valid_animals[, is_valid := NULL]
key(valid_animals)
key(dt)
dt <- dt[valid_animals]
unique(dt[meta=T])[,.N,by="species"]
gc()

#  To add day number, and light phase
#  To add day number, and light phase
dt [,day:=floor(t/days(1))]
dt [,phase:=ifelse(t %% hours(24)>hours(12),"Dark","Light")]
dt [,phase:=factor(phase, levels= c("Dark","Light"))]
gc()

#sleep architecture

pl <- ggetho(dt[xmv(exp)=="1"],
             aes(X=t,y=interaction(id, genotype, sleep_deprived,machine_name, sep = " : "),z=asleep)
) 
pl_overview <- pl + 
  #facet_grid(rep ~ .) +
  ggetho::stat_tile_etho() + 
  stat_ld_annotations() 

pl_overview

pl <- ggetho(dt,
             aes(X=t,y=interaction(id, genotype, sleep_deprived, sep = " : "),z=asleep)
) 
pl_overview <- pl + 
  #facet_grid(rep ~ .) +
  ggetho::stat_tile_etho() + 
  stat_ld_annotations() 

pl_overview

# asleep fraction population graph

pdf("male_male_SD_reps1.pdf") 

pl <- ggetho(dt[t>days(2)&t<days(5) & xmv(status)=="ok"],
             aes(x=t, y=asleep, colour=sleep_deprived))


pl_1 <- pl +
  stat_pop_etho() +
  stat_ld_annotations()+
  ethogram_theme+
  facet_grid(SD_method~genotype)+
  labs(title = "Asleep Fraction Population Graph") +
  xlab("time")+ ylab("Asleep Fraction")+
  scale_fill_manual(values=c("gray","seagreen","green"))+
  scale_colour_manual(values=c("gray","seagreen","green"))+
  scale_y_continuous(limits=c(-0.10,1),expand = c(0.01, 0.01))

pl_1
dev.off()


###boxplots

####calculate mean sleep. You can change "asleep" to other behaviours form line 128-134. You can add "phase" to by= c("id", "phase")

dt_rebound<- dt[t>days(4)&t<days(4.125) &xmv(status)=="ok"]

dt_summary_rebound  = dt_rebound[,
                                 .(sleep_fraction=mean(asleep),
                                   moving_fraction=mean(moving))
                                 ,by=c("id")]#"V1","V2", "phase"

###this is to add columns from metdata file to summary tables

setkeyv(dt_summary_rebound,"id")
metadata<-dt[meta=T]
dt_summary_rebound_2 <- metadata[dt_summary_rebound]

# # # # to compare the difference after sd and control for each genotype and both methods

# Create and save box_plots. This will be saved in working direction 
pdf("boxplot_rebound.pdf") 

box_plot <- ggplot(dt_summary_rebound_2,
                   aes(sleep_deprived, sleep_fraction, fill=sleep_deprived)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(size=2,alpha=.3, height=0)+
  ethogram_theme+ 
  scale_y_continuous(limits=c(-0.3,1))+
  labs(title = "Mean Sleep Fraction During 3h Rebound") +
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette ="Paired")+
  facet_grid(SD_method~genotype)+
  stat_n_text()
box_plot
dev.off()
#mean sleep fraction for all

mean_sleep_fraction<-tapply(X=dt_summary_rebound_2$sleep_fraction,
                            INDEX = list(dt_summary_rebound_2$sleep_deprived,dt_summary_rebound_2$genotype,dt_summary_rebound_2$SD_method),
                            FUN = mean)
print(mean_sleep_fraction)



###stats
###stats for mean sleep fraction

### 2 way anova

### see differences in rebound sleep for each genotype after dynamic SD


dt_summary_rebound_dynamic<-dt_summary_rebound_2[SD_method=='Dynamic'& sleep_deprived=='TRUE']

mean_sleep_fraction_model_dynamic<- aov(sleep_fraction ~ genotype, dt_summary_rebound_dynamic)

summary(mean_sleep_fraction_model_dynamic)



### see differences in rebound sleep for each genotype with or without dynamic SD
#CS
dt_summary_rebound_dynamic_2<-dt_summary_rebound_2[SD_method=='Dynamic'&genotype=='CS']
mean_sleep_fraction_model_dynamic_2<- aov(sleep_fraction ~ sleep_deprived, dt_summary_rebound_dynamic_2)
summary(mean_sleep_fraction_model_dynamic_2) 

#TAR2-TAR3 Δ124/tm6b
dt_summary_rebound_dynamic_2<-dt_summary_rebound_2[SD_method=='Dynamic'&genotype=='TAR2-TAR3 Δ124/tm6b']
mean_sleep_fraction_model_dynamic_2<- aov(sleep_fraction ~ sleep_deprived, dt_summary_rebound_dynamic_2)
summary(mean_sleep_fraction_model_dynamic_2) 

#w+;;, TAR3/tm3,sb
dt_summary_rebound_dynamic_2<-dt_summary_rebound_2[SD_method=='Dynamic'&genotype=='w+;;, TAR3/tm3,sb']
mean_sleep_fraction_model_dynamic_2<- aov(sleep_fraction ~ sleep_deprived, dt_summary_rebound_dynamic_2)
summary(mean_sleep_fraction_model_dynamic_2) 




### see differences in rebound sleep for each genotype after MM and after control treatment

dt_summary_rebound_MM<-dt_summary_rebound_2[SD_method=='MM'& sleep_deprived=='TRUE']

mean_sleep_fraction_model_MM<- aov(sleep_fraction ~ genotype, dt_summary_rebound_MM)

summary(mean_sleep_fraction_model_MM)

mean_sleep_fraction_model_MM1<-lm(sleep_fraction ~ genotype, dt_summary_rebound_MM)
summary(mean_sleep_fraction_model_MM1)

###mean_sleep_fraction_model_MM<-lm(dt_summary_rebound_MM$sleep_fraction~dt_summary_rebound_MM$sleep_deprived+dt_summary_rebound_MM$genotype)

par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=c(2, 0.8,0))
plot(mean_sleep_fraction_model_MM)

summary(mean_sleep_fraction_model_MM)

setwd("/stats")

#name of the file that will be generated
sink(file="rebound_wildtype_mm.txt") 

### see differences in rebound sleep for each genotype with or without MM
#CS
dt_summary_rebound_MM_2<-dt_summary_rebound_2[SD_method=='MM'&genotype=='CS']

mean_sleep_fraction_model_MM_2<- aov(sleep_fraction ~ sleep_deprived, dt_summary_rebound_MM_2)

summary(mean_sleep_fraction_model_MM_2)

#TAR2-TAR3 Δ124/tm6b
dt_summary_rebound_MM_2<-dt_summary_rebound_2[SD_method=='MM'&genotype=='TAR2-TAR3 Δ124/tm6b']

mean_sleep_fraction_model_MM_2<- aov(sleep_fraction ~ sleep_deprived, dt_summary_rebound_MM_2)

summary(mean_sleep_fraction_model_MM_2)

#w+;;, TAR3/tm3,sb

dt_summary_rebound_MM_2<-dt_summary_rebound_2[SD_method=='MM'&genotype=='w+;;, TAR3/tm3,sb']

mean_sleep_fraction_model_MM_2<- aov(sleep_fraction ~ sleep_deprived, dt_summary_rebound_MM_2)

summary(mean_sleep_fraction_model_MM_2)



#####################
##  bout analysis  ##
#####################

# bout analysis using boutAnalysis function

bout_dt_light<- bout_analysis(asleep,
                              dt[is_interpolated==F & t > days(4) & t < days(4.125) & phase=="Light"])

#bout_dt_dark <- bout_analysis(asleep,
#dt[is_interpolated==F & t > days(4) & t < days(4.125) & phase=="Dark"])


## select only sleep bouts (asleep = FALSE defines activity bouts)
sleep_bout_light <- bout_dt_light[asleep==T]
#sleep_bout_dark <- bout_dt_dark[asleep==T]

summary_bout_light <- sleep_bout_light[,
                                       .(  n=.N, 
                                           median_length=median(duration)),
                                       by=id]
summary_bout_light<-rejoin(summary_bout_light)
#
#summary_bout_dark <- sleep_bout_dark[,
.(  n=.N, 
    median_length=median(duration)),
by=id]
#summary_bout_dark<-rejoin(summary_bout_dark)

median_bout_length_during_sleep_rebound <- ggplot(summary_bout_light,
                                                  aes(sleep_deprived, median_length/60, fill=genotype)) + 
  geom_boxplot() +
  geom_jitter(size=2,alpha=.3, height=0)+
  generic_theme +
  scale_y_continuous(limits=c(-0.3,15))+
  facet_grid(SD_method~genotype)+
  labs(title = "Median Sleep bout length during sleep rebound") +
  ylab("Median Sleep Bout Length (min)")+
  stat_n_text()

#scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
median_bout_length_during_sleep_rebound

#### median bout length value
median_bout_length<-tapply(X=summary_bout_light$median_length/60,
                           INDEX=list(summary_bout_light$sleep_deprived,summary_bout_light$genotype,summary_bout_light$SD_method),
                           FUN = median)
print(median_bout_length)

##### stats for median bout length

############### compare rebound sleep after MM SD for each genotype

median_bout_length_MM_genotype<-subset(summary_bout_light,SD_method=='MM'&sleep_deprived=='TRUE')

kruskal.test(median_length~genotype,data=median_bout_length_MM_genotype)

print(pairwise.wilcox.test(median_bout_length_MM_genotype$median_length, 
                           median_bout_length_MM_genotype$genotype,
                           p.adjust.method ="fdr"))



############### compare rebound sleep with or without MM SD for each genotype
############### compare rebound sleep bout length in each genotype with or without MM, use wilcox test as there
#               only 2 groups
# CS
median_bout_length_MM_genotype<-subset(summary_bout_light,SD_method=='MM'&genotype=='CS')
print(pairwise.wilcox.test(median_bout_length_MM_genotype$median_length, 
                           median_bout_length_MM_genotype$sleep_deprived,
                           p.adjust.method ="fdr"))

# TAR2-TAR3 Δ124/tm6b
median_bout_length_MM_genotype<-subset(summary_bout_light,SD_method=='MM'&genotype=='TAR2-TAR3 Δ124/tm6b')
print(pairwise.wilcox.test(median_bout_length_MM_genotype$median_length, 
                           median_bout_length_MM_genotype$sleep_deprived,
                           p.adjust.method ="fdr"))
# w+;;, TAR3/tm3,sb
median_bout_length_MM_genotype<-subset(summary_bout_light,SD_method=='MM'&genotype=='w+;;, TAR3/tm3,sb')
print(pairwise.wilcox.test(median_bout_length_MM_genotype$median_length, 
                           median_bout_length_MM_genotype$sleep_deprived,
                           p.adjust.method ="fdr"))



############### compare rebound sleep after dynamic SD for each genotype

median_bout_length_dynamic_genotype<-subset(summary_bout_light,SD_method=='Dynamic'&sleep_deprived=='TRUE')

kruskal.test(median_length~genotype,data=median_bout_length_dynamic_genotype)

print(pairwise.wilcox.test(median_bout_length_dynamic_genotype$median_length, 
                           median_bout_length_dynamic_genotype$genotype,
                           p.adjust.method ="fdr"))

############### compare rebound sleep with or without dynamic SD for each genotype
############### compare rebound sleep bout length in each genotype with or without dynamic, use wilcox test as there
#               only 2 groups

# CS
median_bout_length_dynamic_genotype<-subset(summary_bout_light,SD_method=='Dynamic'&genotype=='CS')
print(pairwise.wilcox.test(median_bout_length_dynamic_genotype$median_length, 
                           median_bout_length_dynamic_genotype$sleep_deprived,
                           p.adjust.method ="fdr"))



# TAR2-TAR3 Δ124/tm6b
median_bout_length_dynamic_genotype<-subset(summary_bout_light,SD_method=='Dynamic'&genotype=='TAR2-TAR3 Δ124/tm6b')
print(pairwise.wilcox.test(median_bout_length_dynamic_genotype$median_length, 
                           median_bout_length_dynamic_genotype$sleep_deprived,
                           p.adjust.method ="fdr"))

# w+;;, TAR3/tm3,sb
median_bout_length_dynamic_genotype<-subset(summary_bout_light,SD_method=='Dynamic'&genotype=='w+;;, TAR3/tm3,sb')
print(pairwise.wilcox.test(median_bout_length_dynamic_genotype$median_length, 
                           median_bout_length_dynamic_genotype$sleep_deprived,
                           p.adjust.method ="fdr"))



#### median bout number

summary_rebound_bout_number_light <- sleep_bout_light[,
                                                      .(  n=.N),
                                                      by=id]

summary_rebound_bout_number_light<-rejoin(summary_rebound_bout_number_light)

median_bout_number_during_sleep_rebound <- ggplot(summary_rebound_bout_number_light,
                                                  aes(sleep_deprived, n, fill=genotype)) + 
  geom_boxplot() +
  geom_jitter(size=2,alpha=.3, height=0)+
  generic_theme +
  scale_y_continuous(limits=c(-3,15))+
  facet_grid(SD_method~genotype)+
  labs(title = "Median Sleep Bout Number During Sleep Rebound") +
  ylab("Median Sleep Bout Number") +
  stat_n_text()


median_bout_number_during_sleep_rebound

#median_bout_number_during_sleep_rebound
median_bout_number<-tapply(X=summary_rebound_bout_number_light$n,
                           INDEX=list(summary_rebound_bout_number_light$sleep_deprived,
                                      summary_rebound_bout_number_light$genotype,
                                      summary_rebound_bout_number_light$SD_method),
                           FUN = median)
print(median_bout_number)

##### stats for median bout number

##### test difference in rebound bout number in after MM in each genotype

median_bout_number_MM_genotype<-subset(summary_rebound_bout_number_light,SD_method=='MM'&sleep_deprived=='TRUE')

kruskal.test(n~genotype,data=median_bout_number_MM_genotype)

print(pairwise.wilcox.test(median_bout_number_MM_genotype$n, 
                           median_bout_number_MM_genotype$genotype,
                           p.adjust.method ="fdr"))

############### compare rebound sleep with or without MM SD for each genotype
############### compare rebound sleep bout length in each genotype with or without MM, use wilcox test as there
#               only 2 groups

#CS
median_bout_number_MM_genotype<-subset(summary_bout_light,SD_method=='MM'&genotype=='CS')
print(pairwise.wilcox.test(median_bout_number_MM_genotype$n, 
                           median_bout_number_MM_genotype$sleep_deprived,
                           p.adjust.method ="fdr"))
# TAR2-TAR3 Δ124/tm6b

median_bout_number_MM_genotype<-subset(summary_bout_light,SD_method=='MM'&genotype=='TAR2-TAR3 Δ124/tm6b')
print(pairwise.wilcox.test(median_bout_number_MM_genotype$n, 
                           median_bout_number_MM_genotype$sleep_deprived,
                           p.adjust.method ="fdr"))


# w+;;, TAR3/tm3,sb
median_bout_number_MM_genotype<-subset(summary_bout_light,SD_method=='MM'&genotype=='w+;;, TAR3/tm3,sb')
print(pairwise.wilcox.test(median_bout_number_MM_genotype$n, 
                           median_bout_number_MM_genotype$sleep_deprived,
                           p.adjust.method ="fdr"))


############### compare rebound sleep after dynamic SD for each genotype

median_bout_number_dynamic_genotype<-subset(summary_rebound_bout_number_light,SD_method=='Dynamic'&sleep_deprived=='TRUE')

kruskal.test(n~genotype,data=median_bout_number_dynamic_genotype)

print(pairwise.wilcox.test(median_bout_number_dynamic_genotype$n, 
                           median_bout_number_dynamic_genotype$genotype,
                           p.adjust.method ="fdr"))
############### compare rebound sleep with or without dynamic SD for each genotype
#CS
median_bout_number_dynamic_genotype<-subset(summary_bout_light,SD_method=='Dynamic'&genotype=='CS')
print(pairwise.wilcox.test(median_bout_number_dynamic_genotype$n, 
                           median_bout_number_dynamic_genotype$sleep_deprived,
                           p.adjust.method ="fdr"))
#TAR2-TAR3 Δ124/tm6b

median_bout_number_dynamic_genotype<-subset(summary_bout_light,SD_method=='Dynamic'&genotype=='TAR2-TAR3 Δ124/tm6b')
print(pairwise.wilcox.test(median_bout_number_dynamic_genotype$n, 
                           median_bout_number_dynamic_genotype$sleep_deprived,
                           p.adjust.method ="fdr"))
# w+;;, TAR3/tm3,sb
median_bout_number_dynamic_genotype<-subset(summary_bout_light,SD_method=='Dynamic'&genotype=='w+;;, TAR3/tm3,sb')
print(pairwise.wilcox.test(median_bout_number_dynamic_genotype$n, 
                           median_bout_number_dynamic_genotype$sleep_deprived,
                           p.adjust.method ="fdr"))


