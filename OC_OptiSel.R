#### This is the Optimum contribution estimation
### File1 PedKelp file:
### All fndrs, Sire/Dam  0, sex NA
### GPs, Sire is its fndrs, Dam NA, sex according to itself
### ProgSP, Sire is parental MG, Dam is parental FG, sex NA
### ProgSP_GP, Sire being NA, Dam being ProgSP, sex NA

### File 2 phen:
### Merged, to get all the pedigree info (equiGen, etc)
### File 3 kinship:
### I used our outCovComb output estimated at diploid level

rm(list=ls())
setwd("/Users/maohuang/Desktop/Kelp/SugarKelpBreeding/")
#install.packages("here")

here::i_am("TraitAnalyses201003/OptiSel/ReplicationScript.R")
setwd("TraitAnalyses201003/OptiSel/Examples")


### Input the outCovComb_dip
load(here("TraitAnalyses201003/ReorderPedigree","outCovComb_dip_0116_2021.Rdata"))
  ls()
  outCovComb4_dipOrder[1:4,1:4]

### Input GEBVs
phenoBV<-read.csv(here("allBLUPs_PlotsOnly_withSGP_866_AddfndrsMrkData_0116_2021_dip.csv"),sep=",",header=T,row.names=1)

library("optiSel")
library("data.table")
library(ggplot2)

### These were manipulated using the codes at the end to save files

###############################JL codes

load("Prepared_PedKelp_Pedig_phen_0202_2021.RData")

library(stringr)
phenGP <- phen[phen$isCandidate,]

### RM SNE ones
phenGP$loc<-str_split_fixed(rownames(phenGP),"-", 4)[,2] 

GOM<-c("CB","CC","JS","LD","LL","NC","NL","OD","OI","SF","UCONN")
SNE<-c("BL","FI","FW","PI","TI","LR","DB")

phenGP$Region<-ifelse(phenGP$loc%in%GOM,"GOM","SNE")
  head(phenGP)
  phenGP[phenGP$Region=="SNE",]$loc%in%SNE

phenGP<-phenGP[phenGP$Region=="GOM",]  # 405 individuals
phenGP<-phenGP[,-c((ncol(phenGP)-1),ncol(phenGP))] # rm loc and Region cols

grmGP <- fPED[rownames(fPED)%in%rownames(phenGP), colnames(fPED)%in%rownames(phenGP)]
  ### Order phenGP and grmGP?
phenGP<-phenGP[match(rownames(phenGP),rownames(grmGP)),]
  
forOptiSel <- data.frame(Indiv=rownames(phenGP), Trait=phenGP$DWpM,Sex=phenGP$Sex)  ### Trait here!!!
#invisible(capture.output(cand2 <- optiSel::candes(forOptiSel, grm=grmGP, quiet=T)))
cand <- optiSel::candes(forOptiSel, grm=grmGP)
cand$mean
allNe <- c(6, 30, 60, 300, 600, 1200)
#allNe<-c(60,600)
ocList<-NULL
ocGP <- NULL
for (Ne in allNe){
  con <- list(
    ub.grm = max(cand$mean$grm+(1-cand$mean$grm)/(2*Ne), 0)
  )
  fit <- opticont("max.Trait", cand=cand, con=con, trace=F)
  oc<-fit$parent
  print(Ne)
  print(fit$info)
  print(fit$mean)
  oc$Ne<-paste0("Ne",Ne)
  ocList<-rbind(ocList,oc)
  ocGP <- cbind(ocGP, oc[,"oc"])
}

sum(ocList[ocList$Ne=="Ne60",]$oc>0.005)
sum(ocList[ocList$Ne=="Ne600",]$oc>0.005)
data2<-ocList[ocList$Ne=="Ne60",]
data3<-ocList[ocList$Ne=="Ne600",]  #ocList from JL codes estimation

ggplot(data=NULL,aes(x=data2$oc,y=data3$oc))+
  geom_point()

#default was using "cccp2"
################### Run twice for the trait !!!

AshFDWpM_oc<-data2  ## !!!

DWpM_oc<-data2      ## !!!
 
TwoTraits_oc<-AshFDWpM_oc
library(expss)
TwoTraits_oc$oc_using_AshFDWpM<-AshFDWpM_oc$oc
TwoTraits_oc$AshFDWpM<-AshFDWpM_oc$Trait
TwoTraits_oc$oc_using_DWpM<-vlookup(TwoTraits_oc$Indiv,dict=DWpM_oc,result_column="oc",lookup_column="Indiv")
TwoTraits_oc$DWpM<-vlookup(TwoTraits_oc$Indiv,dict=DWpM_oc,result_column="Trait",lookup_column="Indiv")

TwoTraits_oc$loc<-str_split_fixed(rownames(TwoTraits_oc),"-", 4)[,2] 
TwoTraits_oc$Region<-ifelse(TwoTraits_oc$loc%in%GOM,"GOM","SNE")
  TwoTraits_oc[TwoTraits_oc$Region=="SNE",]$loc%in%SNE
  
write.csv(TwoTraits_oc,"Candidates_TwoTraits_oc_0204_2021.csv")

TwoTraits_oc<-TwoTraits_oc[order(TwoTraits_oc$Sex,-TwoTraits_oc$oc_using_DWpM),]
    TwoTraits_oc$DWpM[50]  # 0.05
    
Select_Sex<-TwoTraits_oc[TwoTraits_oc$Sex=="female",]  

Select_DwPM<-Select_Sex[Select_Sex$DWpM>0,]  # DWpM >0
Select_DwPM_oc<-Select_DwPM[order(-Select_DwPM$oc_using_DWpM),]  

Select_DwPM_oc$FemaleRank_DwPM_oc<-c(1:nrow(Select_DwPM_oc))

####Plot the maximum contribution as a function of Ne
colnames(ocGP)<-paste0("Ne",allNe)

plotA<-ggplot(data=NULL,aes(log10(allNe),apply(ocGP, 2, max)))+
  geom_point()+
  ylab("Maximum OC value")+
  xlab("log10(Ne)")+
  theme_bw()+
  theme(panel.grid.major = element_blank())

#scale_x_continuous(trans="log10",breaks=c(),label=c())

# Plot how many GP should have contribution of at least 10% of max contribution
closeToMax <- apply(ocGP, 2, function(vec) sum(vec > max(vec)*0.05))

plotB<-ggplot(data=NULL,aes(log10(allNe),closeToMax))+
  geom_point()+
  xlab("log10(Ne)")+
  ylab("Number of GPs > 5% of maximum OC")+
  theme_bw()+
  theme(panel.grid.major = element_blank())

library(ggpubr)
ggarrange(plotA, plotB,
          ncol = 2, nrow = 1,
          common.legend = TRUE)

##################################
# 
data2<-read.csv("Candidates_TwoTraits_oc_0204_2021.csv",sep=",",header=T)
data2$Trait<-data2$DWpM
data2$oc<-data2$oc_using_DWpM
#### Plot out Trait vs oc separating The 50th value and <=0

AboveZero<-function(oc){
  color<-ifelse(oc>sort(oc,decreasing = TRUE)[50],">Top50oc","<=Top50oc") 
  return(color)
}

### Trait !!!
plot<-function(data){
  plot<-ggplot(data,aes(Trait,oc))+
    geom_point(aes(color=AboveZero(oc)))+
    theme_bw()
  return(plot)
}

# Parent60$Trait<-Parent60$DWpM
# plot1<-plot(Parent60)
# #plot2<-plot(Parent240)
# #plot3<-plot(Parent600)

dataF<-droplevels(data2[data2$Sex=="female",])  # subset sex
dataM<-droplevels(data2[data2$Sex=="male",])  # subset sex


plot2<-plot(dataF)
plot3<-plot(dataM)

### Trait !!!
pdf(file="Candidates_Ne60_FG_MG_DwPM_withOnlyub_GOMonly_0204_2021.pdf",
    width=8.5,height =11)

library(ggpubr)
ggarrange(plot2, plot3,
          labels = c("Ne60_FG", "Ne60_MG"),
          ncol = 1, nrow = 2)
dev.off()

data2<-data2[order(data2$Sex,-data2$Trait,-data2$oc),]
data3<-data3[order(data3$Sex,-data3$Trait,-data3$oc),]

#### Plotting
oc2F<-data2[order(data2$Sex,-data2$oc),][1:100,]
oc3F<-data3[order(data3$Sex,-data3$oc),][1:100,]
sum(oc2F$Indiv%in%oc3F$Indiv)

plot<-ggplot(data=NULL,aes(oc2F$oc,oc3F$oc))+
  geom_point()
print(plot)
#### Done plotting

#### Comparing the list of individuals selected as top 100, for F and M  
F_M_list<-function(data=data){
  Top50F<-data[data$Sex=="female",]$Indiv[1:50]
  Top50M<-data[data$Sex=="male",]$Indiv[1:50]
  Top50F<-sort(Top50F)
  Top50M<-sort(Top50M)
  return(list(Top50F=Top50F,Top50M=Top50M))
}

list<-list(data2,data3) #
tmpF<-NULL
tmpM<-NULL

for (i in 1:length(list)){
  L1F<-F_M_list(list[[i]])$Top50F
  L1M<-F_M_list(list[[i]])$Top50M
  tmpF<-cbind(tmpF,L1F)
  tmpM<-cbind(tmpM,L1M)
}

head(tmpF)
head(tmpM)
identical(tmpF[,1],tmpF[,2])
identical(tmpF[,2],tmpF[,3])




### FILE 1: 
# ### Input Pedi
# pedir<-here("TraitAnalyses201003/ReorderPedigree","Ped_in_Order_866_Individuals_Fndr_New_Order_0116_2021.csv")
# biphasicPedNH<-read.csv(file=pedir,sep=",",header=TRUE,row.names=1)

# ### Input raw Pheno
# load(here("TraitAnalyses201003","dataNHpi_withChk_3_sets_PhotoScore23.rdata"))
# dataNHpi<-dataNHpiBoth_C  

# ### Reformat pedigree matrix to match the package format ###
# cols<-c("Indiv","Sire","Dam","Sex","Born")
# 
# PedKelp<-matrix(nrow=nrow(biphasicPedNH),ncol=length(cols))
# rownames(PedKelp)<-rownames(biphasicPedNH)
# 
# colnames(PedKelp)<-cols
# PedKelp<-as.data.frame(PedKelp)
# 
# PedKelp[,1:3]<-biphasicPedNH[,c(1,3,2)]  ### !!!! biphasicPedNH: fndRow, F(Dame), M(Sire)
#   head(PedKelp)

# ### Add the generation number to PedKelp$Born !!!!!!!!!!
# PedKelp$Born[fndrRow]<-1 # fndrs
# PedKelp$Born[GPRow]<-2
# PedKelp$Born[SPRow]<-3
# PedKelp$Born[SP_GPRow]<-4


# ### Prepare/Pruning the Pedigree file "dataf" (PedKelp)
# PedKelp$Indiv2<-rownames(PedKelp)
# Pedig<-prePed(PedKelp)
# rownames(Pedig)<-Pedig$Indiv2
# 
# PedKelp<-PedKelp[,!colnames(PedKelp)=="Indiv2"]
# Pedig<-Pedig[,!colnames(Pedig)=="Indiv2"]
#   
# ### The sexes turned out to be wrong, after checking their source code!!!
# library(expss)
# Pedig$Sex<-vlookup(rownames(Pedig),dict=PedKelp,result_column = "Sex",lookup_column = "row.names")
# 

### FILE 2:  
### Merge phenotypic file to the pedigree info
# phen<-merge(phenoBV,PedKelp,all.x=TRUE,by.x="row.names",by.y="row.names")
#   head(phen)
#phen<-phen[order(match(phen$Row.names,rownames(phenoBV))),]

### Compute the number of equivalent complete generations in the pedigree of each individual ###
### equiGen: Number of equivalent complete generations, defined as the sum of the proportion of known ancestors over all generations traced,
# Sy    <- summary(Pedig) 
#   head(Sy)
# phen  <- merge(phen, Sy[, c("Indiv", "equiGen")], on="Indiv")

# ##### Checking on completeness
# compl <- completeness(Pedig, keep=phen$Indiv, by="Indiv")
#   head(compl)
#   tail(compl)
# compl <- completeness(Pedig, keep=phen$Indiv, by="Sex")
# library("ggplot2")
# plot<-ggplot(compl, aes(x=Generation, y=Completeness, col=Sex)) + 
#   geom_line()
# print(plot)
# #######
     
#   diag(fPED2)<-NA 
#   mean(c(fPED2),na.rm=TRUE)
#   hist(c(fPED2))
# #### Assessing relationship among Top 20 FGs
#   TwoTraits_oc<-TwoTraits_oc[order(-TwoTraits_oc$oc_using_DWpM),]
#   Top20<-as.character(TwoTraits_oc$Indiv[1:20])
#     str(Top20)
#   fPED3<-fPED[rownames(fPED)%in%Top20,colnames(fPED)%in%Top20]  
#     dim(fPED3)
#   diag(fPED3)<-NA 
#   mean(c(fPED3),na.rm=TRUE)
#     hist(c(fPED3)) 
###########

# phen2 Not keeping the equiGen column   
# cand1 <- candes(phen=phen2, fPED=fPED2, cont=NULL,quiet=T)
#   cand1$mean

### Traditional OCS with upper bounds for not just female, no care sex, contributions ###
  # n<-5 # number of progeny per individual
  # N0<-400 # evaluation population size
  # threshold<-n/(2*N0)  
  #females <- cand$phen$Sex=="female" & cand$phen$isCandidate
  #ub  <- setNames(rep(threshold, sum(females)), cand$phen$Indiv[females])
  
### Define upper limits for kinships and native kinships ###
# Compare_Ne<-function(Ne=60,cand=cand1){
#   library(optiSel)
# ub.fPED  <- cand1$mean$fPED  + (1-cand1$mean$fPED)/(2*Ne*L)
# con1 <- list(
#             ub.fPED    = ub.fPED)
# return(con1)
# }
# #Can set up other constrains in con
# #lb.equiGen = cand$mean$equiGen,
# #ub         = ub, 
# 
# ### Calculate the optimum contribution
# Compare_Ne(60,cand1)
# #Compare_Ne(240,cand)
# #Compare_Ne(600,cand)

# ### Change trait here!!!
# fit60<-opticont("max.DWpM",cand=cand1,con=Compare_Ne(Ne=60,cand=cand1),solver="cccp2",trace=FALSE)  ### choose a trait!  "cccp"
# 
# #fit240<-opticont("max.AshFDwPM",cand=cand,con=Compare_Ne(Ne=240,cand=cand),solver="slsqp",trace=FALSE)  ### choose a trait!  "cccp"
# #fit600<-opticont("max.AshFDwPM",cand=cand,con=Compare_Ne(Ne=600,cand=cand),solver="slsqp",trace=FALSE)  ### choose a trait!  "cccp"
#   fit60$info
#   fit60$mean
#   
# OC<-function(fit){
#   fit$parent$oc
# }
# 
# sum(OC(fit60)>0.01) #319 >0; 20 >0.01   # 306  #515
# #sum(OC(fit240)>0.00) #468 >0; 20 >0.01  # 497  #314
# #sum(OC(fit600)>0.00) #517 >0; 20 >0.01  #280   #255
# 
# Parent60<-fit60$parent
# #Parent240<-fit240$parent
# #Parent600<-fit600$parent


# 
# ### Distinguish each individual to fndr, GP, SP, SP_GP in order to define generation # and sexes of Indiv
# find_Individual<-function(pedNameList){
#   fndrRow<-which(!grepl("FG",pedNameList) & !grepl("MG",pedNameList))
#   GPRow<-which(!grepl("x",pedNameList) & !grepl("UCONN-S",pedNameList) & c(grepl("FG",pedNameList) | grepl("MG",pedNameList)))
#   SPRow<-which(grepl("x",pedNameList))
#   SP_GPRow<-which(grepl("UCONN-S",pedNameList))
#   return(list(fndrRow=fndrRow,GPRow=GPRow,SPRow=SPRow,SP_GPRow=SP_GPRow))
# }
# 
# dataf<-phenoBV #### This is when I do not need 
# #dataf<-PedKelp
# find_Individual<-find_Individual(rownames(dataf))
# 
# fndrRow<-find_Individual$fndrRow
# GPRow<-find_Individual$GPRow
# SPRow<-find_Individual$SPRow
# SP_GPRow<-find_Individual$SP_GPRow
# 
# 
# ### Finding the FG and MG sexes
# find_FGMG<-function(GPNameList){
#   FGP<-GPNameList[grep("FG",GPNameList)]  # FGP list
#   MGP<-GPNameList[grep("MG",GPNameList)]  # MGP list
#   return(list(FGP=FGP,MGP=MGP))
# }
# 
# 
# GPNameList<-rownames(dataf)[c(GPRow,SP_GPRow)]  ### !!!!
# 
# FGP<-find_FGMG(GPNameList)$FGP
# MGP<-find_FGMG(GPNameList)$MGP
# 
# dataf$Sex<-NA
# for (i in 1:nrow(dataf)){
#   if (rownames(dataf)[i]%in%FGP==TRUE){
#     dataf$Sex[i]<-2                           # female
#   }else if(rownames(dataf)[i]%in%MGP==TRUE){
#     dataf$Sex[i]<-1                           # male
#   }else{
#     is.na<-dataf$Sex[i]
#   }
# }
# 
# phen<-dataf
# 
# phen$Sex[phen$Sex==1]<-"male"
# phen$Sex[phen$Sex==2]<-"female"
# phen$isCandidate<-rownames(phen)%in%rownames(phen)[c(GPRow,SP_GPRow)]
# 
# # ######## Keep only the GPs---- this works the same as the phen of all lines,
# phen$Indiv<-rownames(phen)
# phen2<-phen[phen$isCandidate,]
# dim(phen2)
# 
# #### Assessing relationship among 517 GPs
# fPED2<-fPED[rownames(fPED)%in%phen2$Indiv,colnames(fPED)%in%phen2$Indiv]  
# dim(fPED2)



#Parent60$oc240<-vlookup(Parent60$Indiv,Parent240,result_column = "oc",lookup_column = "Indiv")
#Parent60$oc600<-vlookup(Parent60$Indiv,Parent600,result_column = "oc",lookup_column = "Indiv")

### Trait !!!
#Parent60<-Parent60[order(Parent60$Sex,-Parent60$AshFDwPM,-Parent60$oc,-Parent60$oc240,-Parent60$oc600),]  

#write.csv(Parent60,"Candidates_Ne60_240_600_AshFDwPM_withOnly_ub_866Indiv.csv")

# Provide the data.frame, and order by Sex and the column of oc
# data1<-Parent60[order(Parent60$Sex,-Parent60$Trait,-Parent60$oc),]  ### Trait !!!
#data2<-Parent60[order(Parent60$Sex,-Parent60$AshFDwPM,-Parent60$oc240),] ### Trait !!!
#data3<-Parent60[order(Parent60$Sex,-Parent60$AshFDwPM,-Parent60$oc600),] ### Trait !!!


# ## Compare the list for different Nes, with different trait and different lb ub settings
# Model7F<-tmpF
# Model7M<-tmpM
#   identical(Model7M[,1],Model1M[,1])
#   identical(Model7F[,1],Model1F[,1])
  
# AshFDwPM_oc<-read.csv("Candidates_Ne60_240_600_AshFDwPM_withOnly_ub_866Indiv.csv",sep=",",header=T)
# DWpM_oc<-read.csv("Candidates_Ne60_240_600_DWpM_withOnly_ub_866Indiv.csv",sep=",",header=T)    
  
# colkeep<-c("Indiv","AshFDwPM","DWpM","Sex","oc","AshOnly","WWP","PDW")
# TwoTraits_oc<-AshFDwPM_oc[,colkeep]  
# TwoTraits_oc$oc_using_AshFDwPM<-TwoTraits_oc$oc
# library(expss)
# TwoTraits_oc$oc_using_DWpM<-vlookup(TwoTraits_oc$Indiv,dict=DWpM_oc,result_column="oc",lookup_column="Indiv")
#   head(TwoTraits_oc)
# Order data frame


# library(tidyr)
# Longdf<-gather(Parent60,NeValue,OCs,c(oc,oc240:oc600))  # (data, new variable for Cols, values for Cols, Cols)
# Longdf$NeValue[Longdf$NeValue=="oc"]<-"Ne60"
# Longdf$NeValue[Longdf$NeValue=="oc240"]<-"Ne240"
# Longdf$NeValue[Longdf$NeValue=="oc600"]<-"Ne600"
# Longdf$NeValue<-as.factor(Longdf$NeValue)
  # plot<-ggplot(data=Longdf,aes(AshFDwPM,OCs))+
  #   geom_point(aes(color=AboveZero(OCs)))+ 
  #   theme_bw()+
  #   facet_grid(rows=vars(NeValue))
# facet_grid was not able to find "NeValue" ??

# print(plot)
#  
# 
#   plot(OC(fit60))
#   plot(OC(fit240))
#   plot(OC(fit600))
#   plot(OC(fit60),OC(fit600))
# 
#   
  # ### L is the generation interval in years
  # ### approximately set to be 2
  #   L<-2
  #  
  # #### FILE 3:
  # fPED<-outCovComb4_dipOrder
  #   #fPED<-pedIBD(PedKelp)      # This one is identical to using Pedig
  #   #fPED2  <- pedIBD(Pedig)    # only r=0.737 with the cor(c(outCovComb4_dipOrder),c(fPED))
  # 
  # save(phen,fPED,file="Prepared_PedKelp_Pedig_phen_0202_2021.Rdata") 
  
  ### Convert the Indiv from numeric back to character, so that it matches those in fPED
  # phen$IndivNum<-as.factor(phen$Indiv)    
  # phen$Indiv<-phen$Row.names    # has to be "Indiv", matching fPED names
  #   head(phen)
  #   str(phen)
  
# #head(Candidate[Candidate$oc>0.01,c("Sex","AshFDwPM","oc")])   ### choose a trait!
# Candidate_rank<-Candidate60[order(-Candidate60$oc),]
# head(Candidate_rank)
#   
### Making plots, GEBV vs OC, two figure given the Nes, then color code those with larger than 0

### Making plots, OC in different color given the two Nes?

# ####!!!!!!!!!!! Question, why did this not work???
# ### Estimate population means ###
# cont <- agecont(Pedig, use=!is.na(Pedig$Born),maxAge = NA)  ### Did not work ?????
# head(cont) 
# 
# L  <- 1/(4*cont$male[1]) + 1/(4*cont$female[1])  ### Cannot get L due to NAs in the cont step
# 
# cont<-genecont(Pedig,from=NULL,to=NULL)  ### Did not work??
# str(cont)
# 
# 





# Ignore these below:
# ###########OLDER way of estimating the PedKelp
# 
# rownames(PedKelp)<-biphasicPedNH$fndRow
# PedKelp<-as.data.frame(PedKelp)
# PedKelp$Indiv<-c(rownames(biphasicPedNH))
# 
# Name_Order<-c(rownames(biphasicPedNH))
# names(Name_Order)<-c(biphasicPedNH[,1])
# head(Name_Order)
# ### Find the "sire" "dam" name according to the biphasicPedNH for PedKelp
# PedName<-function(row){
#   if(row==0 | is.na(row)){
#     value<-row
#   }else{
#     value<-Name_Order[row]
#   }
#   return(value)
# }
# 
# for (i in 1:nrow(PedKelp)){
#   PedKelp$Sire[i]<-PedName(biphasicPedNH[i,2])
#   PedKelp$Dam[i]<-PedName(biphasicPedNH[i,3])
# }
# 
# # if without FG, it is fndr, 2018 (generation 1)
# # if with FG, MG but no "x" or UCONN, it is GP, 2018 (generation 2)
# # if with "x", it is crosses made in 2018 and 2019 -> need to match to which year (generation 3)
# # if with UCONN, it is UCONN_GP, 2019 (generation 4)
# # 
# find_Individual<-function(pedNameList){
#   fndrRow<-which(!grepl("FG",pedNameList) & !grepl("MG",pedNameList))
#   GPRow<-which(!grepl("x",pedNameList) & !grepl("UCONN-S",pedNameList) & c(grepl("FG",pedNameList) | grepl("MG",pedNameList)))
#   SPRow<-which(grepl("x",pedNameList))
#   SP_GPRow<-which(grepl("UCONN-S",pedNameList))
#   return(list(fndrRow=fndrRow,GPRow=GPRow,SPRow=SPRow,SP_GPRow=SP_GPRow))
# }
# # 
# ### Add the generation number to PedKelp$Born
# find_Individual<-find_Individual(PedKelp$Indiv)
# fndrRow<-find_Individual$fndrRow
# GPRow<-find_Individual$GPRow
# SPRow<-find_Individual$SPRow
# SP_GPRow<-find_Individual$SP_GPRow
# 
# 
# 
# PedKelp$Born[fndrRow]<-1 # fndrs
# PedKelp$Born[GPRow]<-2
# PedKelp$Born[SPRow]<-3
# PedKelp$Born[SP_GPRow]<-4
# 
# is.na<-PedKelp$Sex
# str(PedKelp)
# 
# ### Finding the FG and MG rows
# find_FGMG<-function(GPNameList){
#   FGP<-GPNameList[grep("FG",GPNameList)]  # FGP list
#   MGP<-GPNameList[grep("MG",GPNameList)]  # MGP list
#   return(list(FGP=FGP,MGP=MGP))
# }
# 
# FGP<-find_FGMG(PedKelp$Indiv[c(GPRow,SP_GPRow)])$FGP
# MGP<-find_FGMG(PedKelp$Indiv[c(GPRow,SP_GPRow)])$MGP
# 
# 
# for (i in 1:nrow(PedKelp)){
#   if (rownames(PedKelp)[i]%in%FGP==TRUE){
#     PedKelp$Sex[i]<-2
#   }else if(rownames(PedKelp)[i]%in%MGP==TRUE){
#     PedKelp$Sex[i]<-1
#   }else{
#     is.na<-PedKelp$Sex[i]
#   }
# }
# 
# #### Sorry, Decided to convert the PedKelp back to numeric form
# rownames(PedKelp)<-rownames(biphasicPedNH)
# PedKelp$Indiv<-biphasicPedNH$fndRow
# PedKelp$Sire<-biphasicPedNH$X.1
# PedKelp$Dam<-biphasicPedNH$X.2
# PedKelp[103:105,]
# 
# head(PedKelp)  
# PedKelp[103:105,]





#### Other codes might be useful

# ### Estimate native contributions from pedigrees ###
# Pedig2 <- prePed(Ped, lastNative=1970, thisBreed="Angler", keep=animals)
# BC     <- pedBreedComp(Pedig2, thisBreed="Angler")
# phen   <- merge(phen, BC[, c("Indiv", "unknown", "native")], on="Indiv")
# setnames(phen, old="native", new="pedNC")
# BC[I, 1:4, on="Indiv"]

# ### Estimate native kinships from pedigrees ###
# fPEDN  <- pedIBDatN(Pedig2, thisBreed="Angler", keep.only=phen$Indiv)
# natKin <- fPEDN$Q1/fPEDN$Q2
# natKin[I, I]
# 

#con<-list(uniform="female",ub.sKin=0.010)
# This is assuming female contributions uniform. derive the ub.sKin/ub.fPED. But it needs L !!!! Did not quite work
# Ne <- 60
# con <- list(
#   uniform = "female",
#   ub.fPED = 1-(1-cand$mean$fPED)*(1-1/(2*Ne))^(1/L),
#   lb.equiGen = cand$mean$equiGen
# )

# ### Mate allocation
# 
# con <- list(ub=ub, ub.fSEG=ub.fSEG)
# fit <- opticont("max.EBV", cand, con, solver="cccp")
# 
# Candidate <- fit$parent
# Candidate$n <- noffspring(Candidate, N=200, random=FALSE)$nOff
# Mating      <- matings(Candidate, fSEG)
# head(Mating)
# 
# attributes(Mating)$objval

