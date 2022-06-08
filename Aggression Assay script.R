####aggression assay
library(EnvStats)
dt<-read.csv('aggression_assay.csv')
#### mean fighting frquency
mean_fighting_frequency<-tapply(X=dt$frequency.of.fighting,
                            INDEX = dt$genotype,
                            FUN = mean)
mean_fighting_frequency



#### boxplot median fighting frequency
median_fighting_frequency<-tapply(X=dt$frequency.of.fighting,
                                INDEX = dt$genotype,
                                FUN = median)

median_fighting_frequency
median_fighting_boxplot<-boxplot(dt$frequency.of.fighting~dt$genotype,
                                 xlab = 'Genotype',
                                 ylab = 'Median Frequency of Fighting',
                                 main='Median Frequency of Fighting for Each Genotype',
                                 col=c('pink','light green','light blue'))

median_fighting_boxplot

kruskal.test(dt$frequency.of.fighting~dt$genotype)

#### mean fighting index

mean_fighting_index<-tapply(X=dt$fighting.index,
                                INDEX = dt$genotype,
                                FUN = mean)
mean_fighting_index

#### boxplot median fighting index

median_fighting_index<-tapply(X=dt$fighting.index,
                                  INDEX = dt$genotype,
                                  FUN = median)
median_fighting_index

median_fighting_boxplpt<-boxplot(dt$fighting.index~dt$genotype,
                                 xlab = 'Genotype',
                                 ylab = 'Median Fighting Index',
                                 main='Median Fighting Index for each Genotype',
                                 col=c('pink','light green','light blue'))

kruskal.test(dt$fighting.index~dt$genotype)


#### fighting length
dt_length<-read.csv('Aggression Assay Fighting Length.csv')

#### mean fighting length

mean_fighting_length<-tapply(X=dt_length$Length.s.,
                                INDEX = dt_length$Genotype,
                                FUN = mean)
mean_fighting_length

#### median fighting lenghth

median_fighting_length<-tapply(X=dt_length$Length.s.,
                             INDEX = dt_length$Genotype,
                             FUN = median)
median_fighting_length

median_fighting_boxplpt<-boxplot(dt_length$Length.s.~dt_length$Genotype,
                                 xlab = 'Genotype',
                                 ylab = 'Median Length of Fighting (s)',
                                 main='Median Length of Fighting for Each Genotype',
                                 col=c('pink','light green','light blue'),
                                 ylim=c(0,30))
kruskal.test(dt_length$Length.s.~dt_length$Genotype)
