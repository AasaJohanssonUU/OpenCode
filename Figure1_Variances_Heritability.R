## Script for Figure 1 in the article.

## The frequencies have been estimated using plink 2.0.
## Eigen and CADD values has been added for each SNV:
## Starting from the frequenices for the varints in variable: freq with with 8 columns:
## Chr  RSid    Position    MAF   N_Alleles   Chr:Pos CADD    Eigen


##########   Figure 1A ##############

temp<-hist(freq[,4], breaks=50, plot=F)$density
temp<-temp/sum(temp)*100
barplot(temp, ylab="Proportion (%)", xlab="MAF bins (0-50%)", space=0.4, col="white")

##########   Figure 1 D, E, F ##############


# Allele counts per MAF-bin recalculated to variances per MAF-bin
#N=1,484, B) N=14,844, and C) N=148,435).

maf1<-floor(freq[,4]*100)+0.99999999999999999  ## Make bins for the MAF

#maf1<-maf
freq1<-freq ## Just to fit with the script below

## Estimate the variance for the frequencies
## Variance = 𝑠2=Σ(𝑥𝑖−𝑥m)^2/(𝑛−1)
## We assumre the genotyep is 1
# xm=the mean = the frequency
# Xi = 1 eller 0


N1<-round(freq1[,4]*freq1[,5],0) ## Number of rare alleles
N2<-freq1[,5]-N1 ## Nubmer of common allalaes

## Genotype variance per SNV:
## Variance = 𝑠2=Σ(𝑥𝑖−𝑥m)^2/(𝑛−1)
freq1$var<-(N1*((freq1[,4]-1)^2)+ N2*((freq1[,4]-0)^2))/(freq1[,5]-1)
# Where (N1*((freq1[,4]-1)^2) N1, the number of alternative alleles. And  abs(freq1[,4]-1) is 𝑥𝑖−𝑥m



## Phenotype variance alternative 2
freq1$varPHC<- (freq1[,7]^2*N1*((freq1[,4]-1)^2)+ freq1[,7]^2*N2*((freq1[,4]-0)^2))/(freq1[,5]-1)
freq1$varPHE<- (freq1[,8]^2*N1*((freq1[,4]-1)^2)+ freq1[,8]^2*N2*((freq1[,4]-0)^2))/(freq1[,5]-1)

V1<-{}
V1PHC<-{}
V1PHE<-{}
for(i in 1:50){
    ind1<-which(maf1==i) ## SNPs withink the maf-bin
    V1<-c(V1, sum(freq1$var[ind1])) ## Summ variances genotypes
    V1PHC<-c(V1PHC, sum(freq1$varPHC[ind1], na.rm=T)) ## Summ variacnes Phentoyps (CADD)
    V1PHE<-c(V1PHE, sum(freq1$varPHE[ind1], na.rm=T)) ## Summ variacnes Phentoyps (Eigen)
    
}
    
# Figure 1D:
barplot(V1/sum(V1), ylab="Variance (%)", xlab="MAF bins (0-50%)", space=0.4, col="white")

#Figure 1E:
barplot(V1PHC/sum(V1PHC), ylab="Variance (%)", xlab="MAF bins (0-50%)", space=0.4, col="white")

## Figure 1F:
barplot(V1PHE/sum(V1PHE), ylab="Variance (%)", xlab="MAF bins (0-50%)", space=0.4, col="white")




########### Figure 1C ########################


#Split data into rare and common - CADD
rare<-freq[,4]<0.024
rare[which(rare==TRUE)]<-"rare"
rare[which(rare==FALSE)]<-"common"
cadd_bin<-floor(freq$cadd.ind..3./10)
cadd_bin[which(cadd_bin>4)]<-4
tab<-table(rare, cadd_bin)
print(chisq.test(tab))
for(i in 1:5){
tab[,i]<-tab[,i]/sum(tab[,i])
}

colnames(tab)<-c("<10", "10-20", "20-30", "30-40", ">40")
cadd_tab<-tab


#Split data into rare and common - Eigen
eigen_bin<-floor(freq$eigen.ind..4./10)
eigen_bin[which(eigen_bin>4)]<-4

tab<-table(rare, eigen_bin)
print(chisq.test(tab))
for(i in 1:5){
tab[,i]<-tab[,i]/sum(tab[,i])
}
colnames(tab)<-c("<10", "10-20", "20-30", "30-40", ">40")
eigen_tab<-tab


barplot(cbind(cadd_tab, eigen_tab)*100, angle = 45, ylab="Proportion (%)", xlab="Cadd bins              Eigen bins", legend.text=rownames(tab))




############## Figure 1B #######################

## one frequency file was created (using plink) for each individual as well as one for whole cohort.
## For each individuals the number of rare alleles withink each of 50 frequency bins 0-1%; 1-2%, ...  was counted for each individual and saved as a matririx
## With one columnt per frequency bin, and one row per indvidual. Variable name: freqmat



weight_mat<-freqmat

## Weight the columns for each individual
for(i in 1:NCOL(weight_mat)){
    weight_mat[,i]<-freqmat[,i]/sum(freqmat[,i])
}
M<-{}
Sd<-{}

## Loop across the freq-bins and create mean and SDs across alla individuals
for(i in 1:50){
    M<-c(M,mean(as.numeric(as.matrix(weight_mat[i,1:NCOL(weight_mat)])), na.rm=T))
    Sd<-c(Sd,sd(as.numeric(as.matrix(weight_mat[i,1:NCOL(weight_mat)])), na.rm=T))
}
## Change to %
M<-M*100
Sd<-Sd*100

mp <- barplot(M, ylab="Burden (%)", xlab="MAF bins (0-50%)", ylim=c(0,3.5), space=0.4, col="white")
CI1<-M-(Sd*1.95)
CI2<-M+(Sd*1.95)
arrows(x0 = mp, y0=CI1,x1=mp,y1=CI2, length = 0.02, angle=90, code=3)  ## Change length to fit the figure.


