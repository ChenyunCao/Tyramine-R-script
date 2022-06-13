####aggression assay
library(EnvStats)
library(gplots)
library(ggplot2)
dt<-read.csv('aggression assay.csv')
#### mean fighting frquency
mean_fighting_frequency<-tapply(X=dt$number.of.fighting,
                            INDEX = dt$genotype,
                            FUN = mean)
mean_fighting_frequency

number_model<-lm(dt$number.of.fighting~dt$genotype)
summary(frequency_model)

pl_mean_fighting_frequency<-barplot(mean_fighting_frequency,
                                    xlab='Genotype',
                                    ylab = 'Mean Fighting Frequency')


#### boxplot median fighting frequency
#median_fighting_frequency<-tapply(X=dt$number.of.fighting,
                                INDEX = dt$genotype,
                                FUN = median)

median_fighting_frequency
median_fighting_boxplot<-boxplot(dt$number.of.fighting~dt$genotype,
                                 xlab = 'Genotype',
                                 ylab = 'Median Frequency of Fighting',
                                 main='Median Frequency of Fighting for Each Genotype')
                            
median_fighting_boxplot

kruskal.test(dt$number.of.fighting~dt$genotype)

#### mean fighting index

mean_fighting_index<-tapply(X=dt$fighting.index,
                                INDEX = dt$genotype,
                                FUN = mean)
mean_fighting_index

fighting_index_model<-lm(dt$fighting.index~dt$genotype)
summary(fighting_index_model)

pl_mean_fighting_index<-barplot(mean_fighting_index,
                                xlab = 'Genotype',
                                ylab = 'Mean Fighting Index',
                                ylim = c(0,0.4))

#### boxplot median fighting index

#median_fighting_index<-tapply(X=dt$fighting.index,
                                  INDEX = dt$genotype,
                                  FUN = median)
median_fighting_index

median_fighting_boxplpt<-boxplot(dt$fighting.index~dt$genotype,
                                 xlab = 'Genotype',
                                 ylab = 'Median Fighting Index',
                                 main='Median Fighting Index for each Genotype')
                                

kruskal.test(dt$fighting.index~dt$genotype)


#### fighting length
dt_length<-read.csv('Aggression Assay Fighting Length.csv')

#### mean fighting length

mean_fighting_length<-tapply(X=dt_length$Length.s.,
                                INDEX = dt_length$Genotype,
                                FUN = mean)
mean_fighting_length

fighting_length_model<-lm(dt_length$Length.s.~dt_length$Genotype)
summary(fighting_length_model)

pdf('mean fight length.pdf')

pl<-barplot(mean_fighting_length,
            xlab = 'Genotype',
            ylab = 'Mean Fighting Length(s)',
            ylim = c(0,25))
pl

dev.off()

#### median fighting lenghth

#median_fighting_length<-tapply(X=dt_length$Length.s.,
                             INDEX = dt_length$Genotype,
                             FUN = median)
median_fighting_length

median_fighting_boxplpt<-boxplot(dt_length$Length.s.~dt_length$Genotype,
                                 xlab = 'Genotype',
                                 ylab = 'Median Length of Fighting (s)',
                                 main='Median Length of Fighting for Each Genotype',
                                 ylim=c(0,25))

kruskal.test(dt_length$Length.s.~dt_length$Genotype)
pairwise.wilcox.test(dt_length$Length.s.,
                     dt_length$Genotype,
                     p.adjust.method ="fdr")
