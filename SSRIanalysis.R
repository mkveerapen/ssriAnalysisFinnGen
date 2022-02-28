#Exploring longitudinal data for SSRI, antibiotics, and painkillers
#analyst: Kumar Veerapen, PhD
#Date: February 19, 2020
#Logs: Evernote/InterventionAnalysis/Searching longitudinal phenotypes

#working dir on vm 
#gcloud beta compute --project "finngen-refinery-dev" ssh --zone "europe-west1-b" "veerapen"

setwd("/mnt/disks/workDir/r7")

#### ADHERENCE
#totalquantity/totaldays
#1.1 and above is overbuying


###################################
######### Loading libs ############
###################################
library(data.table)
library(ggplot2)
library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)
library(RNOmni)
library(VennDiagram)
library(ggfortify)
library(survival)
library(gridExtra)
library(grid)
library(tidyverse)
library(sunburstR)


#In various European locales, as the comma character serves as decimal point, the read.csv2 function should be used instead.
#endpoints <- read.csv2("FINNGEN_ENDPOINTS_DF5_V2_2020-02-11.csv", na.strings="-9", header=T)

data <- fread("zcat src/finngen_R7_detailed_longitudinal.txt.gz", data.table=F, tmpdir=".")
#data <- fread("zcat finngen_R5_v2_detailed_longitudinal.gz", data.table=F)
#data <- fread("zcat src/finngen_R6_v2_detailed_longitudinal.gz", data.table=F)
#phenos.cov.r5 <- fread("zcat ../r5/20200713_R5_COV_PHENO_KV.txt.gz", data.table=F)

nrow(data)
#[1] 84477978
# R6 115740625 
# R7 148237657
# 53031209???

table(data$SOURCE)
#CANC    DEATH    INPAT  OPER_IN OPER_OUT   OUTPAT PRIM_OUT    PURCH
#46876    35145  2856550  1898128  4155526 10902854 21496596 42776151
#REIMB
#310152

### R6 
#    CANC    DEATH    INPAT  OPER_IN OPER_OUT   OUTPAT PRIM_OUT    PURCH
#67045    52681  3826868  2383206  5403973 14928598 28977406 59605493
#REIMB
#495355

data.purchases <- data[data$SOURCE=="PURCH",]
data.purchases$APPROX_EVENT_DAY <- as.Date(data.purchases$APPROX_EVENT_DAY)
#(unique(data.purchases[nchar(data.purchases$CODE3)!=6 & !is.na(data.purchases$CODE3),]$CODE3) )
#the 5 characters are either accidental movement of the columns to fix later or that they need to be padded.
#Ari mentioned that the VNR codes are 6 digits and no 5 digit ones exist. The 7 digit ones are highly secured prescriptions
#https://wiki.vnr.fi/?page_id=36
#7 digit codes and what https://www.kela.fi/documents/10180/3612716/Korvattavat%20potilaskohtaiset%20erityislupavalmisteet%20%28pdf%29/02594f8d-ded2-4d71-92cc-3f3e03402394

data.purchases$CODE3 <- str_pad(data.purchases$CODE3, width=6, side="left", pad="0")

### Kela drug purchase register
#Field |Register        |Register long name         |
 # ------|----------------|---------------------------|
#  SOURCE|PURCH           |Kela drug purchase register|
  
 # Field |Code in register|Code description.          |
#  ------|----------------|---------------------------|
#  CODE1 | ATC_CODE       | ATC code                  |
#  CODE2 | SAIR 	        | Kela reimbursement code   |
#  CODE3 | VNRO	          | Product number		        |
#  CODE4 | PLKM           | Number of packages        |

#Kela reimbursement code = there are only certain severe disease. Special reimbursement code. What was drug prescribed for?

#vnr codes which contain the drug dosage and packaging
vnr.codes <- read.table("~/workDir/drugIdentities/vnr_vahvuus_pkoko_Arille_210220.txt", sep="\t", header=T, colClasses="character")
vnr.codes$vnr <- str_pad(vnr.codes$vnr, width=6, side="left", pad="0")


###Switch around the ATC Codes for the following with the following commented lines
#N06AB SSRI
#N06A is antidepressants 
#J01C for antibiotics (beta-lactams); J01 is antibiotics general use -- on Veerapen
#M01A painkillers, Non-steroidal -- on tyo high mem
#use ^<code>
#N05BA

#### CHECK ATC code here! Some of the ATC codes are MISSING! That caused the next step to lose out on a lot of data
vnr.codes.drug <- vnr.codes[grepl("^N06A", vnr.codes$ATC) ,]
#| grepl("^N05AN", vnr.codes$ATC)

#### CHECK ATC code here! 
data.purchases.drug <- data.purchases[grepl("^N06A", data.purchases$CODE1) ,]
#| grepl("^N05AN", data.purchases$CODE1)
#TO INCLUDE Lithium, see N05AN - Lithium


## For testing only
#vnr.codes.drug.test <- vnr.codes[grepl("^N05BA", vnr.codes$ATC),]
#data.purchases.drug.test <- data.purchases[grepl("^N05BA", data.purchases$CODE1),]

data.purchases.drug <- data.purchases.drug[order(data.purchases.drug$FINNGENID, data.purchases.drug$APPROX_EVENT_DAY),]
data.purchases.drug <- data.purchases.drug[order(data.purchases.drug$CODE3),]
vnr.codes <- vnr.codes[order(vnr.codes$vnr),]
data.purchases.vnr <- 
  #merge(by.x="CODE3", by.y="vnr", x=data.purchases.drug, y=vnr.codes, all.x=T) #left join -- has all and VNR codes which are not curated
   #merge(by.x="CODE3", by.y="vnr", x=data.purchases.drug, y=vnr.codes, all.y=T) #right join
  #merge(by.x="CODE3", by.y="vnr", x=data.purchases.drug, y=vnr.codes, all=T) #-- has most and VNR codes which are not curated
  merge(by.x="CODE3", by.y="vnr", x=data.purchases.drug, y=vnr.codes) 

#data.purchases.vnr <- data.purchases.vnr[!is.na(data.purchases.vnr$ATC),] ##this will remove the ATC codes that I initially filtered for!

nrow(data.purchases.drug)
#[1] 1849497
#R6 2452329
nrow(data.purchases.vnr)
#[1] 1849422
#only 75 missing
#R6 2452228
# missing 101 only 

#birth year
#APROX_EVENT_DATE - EVENT_AGE
#test <- data.purchases.vnr %>% 
#  arrange(APPROX_EVENT_DAY) %>%
#  distinct(FINNGENID, .keep_all = TRUE) %>% 
#  mutate(bday = lubridate::as_date(APPROX_EVENT_DAY - lubridate::dyears(EVENT_AGE)))

#####################
##### AGE Filter ####
#####################
### NOT DONE ON 10/12/2020
phenos.cov <- fread("zcat src/R6_COV_PHENO_V0.txt.gz", data.table=F)
#nrow(phenos.cov[phenos.cov$AGE_AT_DEATH_OR_NOW>35,])
#[1] 226950 from 256404 (29454 removed)

phenos.cov <- phenos.cov[phenos.cov$AGE_AT_DEATH_OR_NOW>35,]
id.list <- phenos.cov$FINNGENID
length(id.list)
data.purchases.vnr <- data.purchases.vnr[order(data.purchases.vnr$FINNGENID),]

#nrow(data.purchases.vnr[data.purchases.vnr$FINNGENID %in% id.list,])
#data.purchases.vnr[data.purchases.vnr$FINNGENID %in% id.list,]
#Adding the >35 year old age threshold from Juha
data.purchases.vnr <- data.purchases.vnr[data.purchases.vnr$FINNGENID %in% id.list,]
#103619 to 88428 individuals
#####################
#####################
#####################

#What's missing?
intersect.vnr.diff.vnr.drug <- intersect((unique(data.purchases.vnr$CODE3)), (unique(data.purchases.drug$CODE3) ))
missing.vnr.diff.vnr <- setdiff((unique(data.purchases.vnr$CODE3)), (unique(data.purchases.drug$CODE3) ) )
missing.vnr.diff.drug <- setdiff((unique(data.purchases.drug$CODE3)), (unique(data.purchases.vnr$CODE3) ) )
#10 (9 in R5) VNR codes are not there

#write.table(missing.vnr.diff.drug, "R5_missingVNRssriAnalysis.txt", quote=F, sep = "\t", row.names = F, col.names = F)
write.table(missing.vnr.diff.drug, "R6_missingVNRssriAnalysis.txt", quote=F, sep = "\t", row.names = F, col.names = F)

length(unique(data.purchases.drug$CODE3) )
#[1] 736
# R6 758
length(unique(data.purchases.vnr$CODE3) )
#[1] 727
# R6 748

## what are the drug families?
a <- c("N06AA", "N06AB", "N06AF", "N06AG", "N06AX")
#a <- c("N06AA", "N06AB", "N06AF", "N06AG", "N06AX", "N05AN01")

b <- c("MRI", "SSRI", "MAO.I", "MAO.A.I", "Others")
#b <- c("MRI", "SSRI", "MAO.I", "MAO.A.I", "Others", "Lithium")

ATC.codes <- as.data.frame(cbind(a,b))
colnames(ATC.codes) <- c("ATC.codes", "drugFamily") 
data.purchases.vnr$drug.Family.code <- 0
data.purchases.vnr$drug.Family <- 0

for(j in 1:nrow(ATC.codes)){
  tempName <- paste("^", as.character(ATC.codes[j,1]), sep="") 
  tempName.fam <- as.character(ATC.codes[j,2])
  print(tempName)
  data.purchases.vnr$drug.Family.code <-
    ifelse(grepl(tempName, data.purchases.vnr$CODE1), as.character(ATC.codes[j,1]), data.purchases.vnr$drug.Family.code)
  data.purchases.vnr$drug.Family <-
    ifelse(grepl(tempName, data.purchases.vnr$CODE1), tempName.fam, data.purchases.vnr$drug.Family)
}

nrow(data)
#[1] 84477978
# R6 115740625
nrow(data.purchases)
#[1] 42776151
# R6 59605493
nrow(data.purchases.drug)
#[1] 1849497
# R6 2452329
nrow(data.purchases.drug[!duplicated(data.purchases.drug$FINNGENID),])
#[1] 81096
# R6 103280
nrow(data.purchases.vnr)
#1849422
# R6 2452228 
nrow(data.purchases.vnr[!duplicated(data.purchases.vnr$FINNGENID),])
#81096 (no one is missing!)
# R6 103280 (no one is missing! ) 

##############################
###### with drug dosage ######
##############################
todaysDate <- format.Date(Sys.Date(), "20%y%m%d")
rm(data, data.purchases, phenos.cov)
save.image(paste(todaysDate, ".RData", sep=""))

for(i in 1:nrow(ATC.codes)) {
  print(as.character(ATC.codes[i,2]))
  data.purchases.vnr.drug <- data.purchases.vnr[data.purchases.vnr$drug.Family.code == as.character(ATC.codes[i,1]),]
  #catching empties
  if(dim(data.purchases.vnr.drug)[1]==0){
    next
  }else{
    drug.purchase.freq <- as.data.frame(table(data.purchases.vnr.drug$FINNGENID) )
    drug.purchase.freq$label <- as.character(ATC.codes[i,1])
    
    unique.drug <- data.purchases.vnr.drug[!duplicated(data.purchases.vnr.drug$FINNGENID),]
    unique.drug.drug <- data.purchases.vnr.drug[!duplicated(data.purchases.vnr.drug$CODE3),]
    
    ### Bursts of purchases
    data.purchases.vnr.drug <- data.purchases.vnr.drug[order(data.purchases.vnr.drug$FINNGENID, 
                                                             data.purchases.vnr.drug$APPROX_EVENT_DAY),]
    
    data.purchases.vnr.bursts <- as.data.frame(data.purchases.vnr.drug %>% 
                                                 group_by(FINNGENID) %>%
                                                 mutate(daysPassed.diff=c(NA,diff(APPROX_EVENT_DAY))) )
    
    data.purchases.vnr.bursts <- data.purchases.vnr.bursts[!is.na(data.purchases.vnr.bursts$daysPassed.diff),]
    
    summary(data.purchases.vnr.bursts$daysPassed.diff)
    #Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    #0.0    35.0    84.0   193.2   112.0  8724.0
    
    quantile(data.purchases.vnr.bursts$daysPassed.diff, c(.05, .1, .2, .8, .9, .95))
    #5%    10%    20%    80%    90%    95%
    #  14.00  17.00  30.00 125.00 224.00 652.75
    
    median.plotting <- 
      as.data.frame(data.purchases.vnr.bursts %>%
                      group_by(FINNGENID) %>%
                      summarise(M=mean(daysPassed.diff), Med=median(daysPassed.diff), Q1=quantile (daysPassed.diff, probs=0.25), 
                                Q2=quantile(daysPassed.diff, probs=0.50), Q3=quantile(daysPassed.diff, probs=0.75), 
                                Q4=quantile(daysPassed.diff, probs=1.00)) ) 
    
    #pdf(paste("R5.finngenFreq.", as.character(ATC.codes[i,2]),".vnr.pdf", sep=""))
    pdf(paste(todaysDate, "_R6.finngenFreq.", as.character(ATC.codes[i,2]),".vnr.pdf", sep=""))
    
    a <- ggplot(drug.purchase.freq, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_col() +
      labs(title=paste(as.character(ATC.codes[i,2]), "purchase"))  +
      theme_classic()
    print("purchase plotting")
    
    b <- ggplot(drug.purchase.freq, aes(x=Freq)) + geom_bar() +
      labs(title=paste(as.character(ATC.codes[i,2]), "purchase Freq dist") ) +
      theme_classic()
    print("purchase freq dist plot")
    
    c <- ggplot(drug.purchase.freq, aes(x=label, y=Freq)) +
      geom_boxplot() +
      stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
      geom_jitter(shape=16, position=position_jitter(0.02)) +
      labs(title=paste(as.character(ATC.codes[i,2]), "purchase Freq dist") ) +
      theme_classic()
    print("purchase freq dist 2")
    
    d <- ggplot(median.plotting, aes(x=reorder(FINNGENID, Med), y=Med)) +
      geom_point(size=.2, shape=21, fill="white") +
      geom_errorbar(aes(ymin=Q1, ymax=Q3), width=.01, position = position_dodge2(width = 0.05, padding = 0.05), linetype="4C88C488") +
      theme_classic()
    print("median plotting")
    
    e <- ggplot(data.purchases.vnr.bursts[!is.na(data.purchases.vnr.bursts$daysPassed.diff),], 
                aes(x=reorder(FINNGENID, daysPassed.diff, FUN=median), y=daysPassed.diff)) +
      geom_boxplot(outlier.size = 0) +
      # stat_summary(fun.y=median, geom="point", shape=23, size=4) +
      labs(title=paste(as.character(ATC.codes[i,2]), "purchase Freq dist") ) +
      theme_classic()
    print("purchase freq dist 3 plotting")
    
    print(a)
    print(b)
    print(c)
    print(d)
    print(e)
    
    dev.off()
  }
 
  print("done")
  print(as.character(ATC.codes[i,1]))
}

# pkoko_num number of pills in packages
# vahvuus_num amount in mg
data.purchases.vnr$pkoko_vahvuus <- as.numeric(data.purchases.vnr$pkoko_num) * as.numeric(data.purchases.vnr$vahvuus_num)
data.purchases.vnr <- data.purchases.vnr[order(data.purchases.vnr$FINNGENID, data.purchases.vnr$APPROX_EVENT_DAY),]


### Bursts of purchases
data.purchases.vnr.bursts <- as.data.frame(data.purchases.vnr %>% 
                                             group_by(FINNGENID) %>%
                                             mutate(daysPassed.diff=c(NA,diff(APPROX_EVENT_DAY))) )

data.purchases.vnr.bursts <- data.purchases.vnr.bursts[!is.na(data.purchases.vnr.bursts$daysPassed.diff),]

### How about what types of drugs?
drugType.freq <- as.data.frame(table(data.purchases.vnr$CODE3))
colnames(drugType.freq) <- c("VNRO", "Freq")
drugType.freq <- merge(by.x="VNRO", by.y="CODE3", x=drugType.freq, y=data.purchases.vnr.bursts) 
drugType.freq <- drugType.freq[order(-drugType.freq$Freq),]

#would this normalise the data?
drugType.freq$dosePerDay <- as.numeric(drugType.freq$pkoko_vahvuus) / as.numeric(drugType.freq$daysPassed.diff)

# number of unique drug names
table(as.data.frame(table(drugType.freq$valmiste))$Freq!=0)
#TRUE
#142
# R6 143 
drugName <- as.data.frame(table(drugType.freq$valmiste))
drugName <- drugName[drugName$Freq!=0,]

drugList <- drugName$Var1
drugName.ATC <- drugType.freq[!is.na(drugType.freq$valmiste) & !duplicated(drugType.freq$valmiste),][,c("valmiste", "ATC")]


columnClasses.drugSummary <- c("numeric", 
                              "factor",
                              "numeric",
                              "factor",
                              "numeric",
                              "character",
                              "numeric",
                              "character",
                              "numeric",
                              "numeric", 
                              "character",
                              "character",
                              "character"
                              )
columnNames.drugSummary<- c("numberOfVNR",
            "VNR.code", 
            "mean.EVENT_AGE",
            "drug.Name",
            "mean.dosePerDay",
            "quartiles.dosePerDay",
            "mean.daysPassed",
            "quartiles.daysPassed",
            "numberOfPeoplePurchase",
            "numberOfPurchases",
            "ATC.code",
            "drug.Family.code",
            "drug.Family"
            )
drugSummary <- read.table(text = "",
                          colClasses = columnClasses.drugSummary,
                          col.names = columnNames.drugSummary)

for(i in 1:length(drugList)) {
  temp.drug <- drugType.freq[drugType.freq$valmiste %in% as.character(drugList[i]),]
  temp.drug <- temp.drug[order(-temp.drug$dosePerDay),]
  print(as.character(drugList[i]) )
  print(nrow(temp.drug))
  
  #what VNR codes are present?
  numberOfVNR <- length(unique(temp.drug$VNRO)) #number of VNRO
  VNR.code <- paste(as.character(unique(temp.drug$VNRO)), collapse=",") #what VNRO separated by ","
  #or drop factors unused
  temp.drug$VNRO <- temp.drug$VNRO[, drop=TRUE]
  drug.purchase.freq.temp <- as.data.frame(table(temp.drug$VNRO))
  
  #mean Age of event
  mean.EVENT_AGE<- mean(as.numeric(temp.drug$EVENT_AGE))
  
  #drugName
  drug.Name<- as.character(unique(temp.drug$valmiste) ) 
  
  #dosePerDay mean and sd
  mean.dosePerDay <- paste(mean(as.numeric(temp.drug$dosePerDay)), sd(as.numeric(temp.drug$dosePerDay)), sep=",")
  quartiles.dosePerDay <- paste(t(as.data.frame(quantile(temp.drug$dosePerDay, c(0,.25, .5, .75, 1.0), na.rm = TRUE))), collapse=",")
  
  #daysPassed
  mean.daysPassed <- paste(mean(as.numeric(temp.drug$daysPassed.diff)), sd(as.numeric(temp.drug$daysPassed.diff)), sep=",")
  quartiles.daysPassed <- paste(t(as.data.frame(quantile(temp.drug$daysPassed.diff, c(0,.25, .5, .75, 1.0), na.rm = TRUE))), collapse=",")
  
  #number of Purchasers (people)
  numberOfPeoplePurchase <- length(unique(temp.drug$FINNGENID))
  
  #datapoints available
  numberOfPurchases <- nrow(temp.drug)
  
  #ATC code
  ATC.code <- as.character(unique(temp.drug$CODE1))
  drug.Family.code <- as.character(unique(temp.drug$drug.Family.code))
  drug.Family <- as.character(unique(temp.drug$drug.Family))
  
  drugSummary <- rbind(drugSummary, cbind(numberOfVNR,
                             VNR.code, 
                             mean.EVENT_AGE,
                             drug.Name,
                             mean.dosePerDay,
                             quartiles.dosePerDay,
                             mean.daysPassed,
                             quartiles.daysPassed,
                             numberOfPeoplePurchase,
                             numberOfPurchases,
                             ATC.code, 
                             drug.Family.code,
                             drug.Family ) )
  
  pdf(paste(todaysDate, "_R6.vnr.",as.character(drugList[i]), "_", nrow(temp.drug), "_ATC", ATC.code, "_", drug.Family, ".pdf", sep="") )
  
  plot1 <- ggplot(drug.purchase.freq.temp, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_col() +
    labs(title=paste(drug.Name, "purchase", sep=" "))  +
    theme_classic()
  print(plot1)
  
  plot2 <- ggplot(temp.drug[!is.na(temp.drug$dosePerDay),], 
         aes(x=reorder(FINNGENID, dosePerDay, FUN=median), y=dosePerDay)) +
    geom_boxplot(outlier.size = 0) +
    # stat_summary(fun.y=median, geom="point", shape=23, size=4) +
    labs(title=paste("R6", drug.Name, "purchase Freq dist", sep = " ") ) +
    theme_classic()
  print(plot2)
  
  dev.off()
  
  
  print("DONE")
  #what's package size, dosage? 
  #title needs to be stitched with number of times
}

###package size and dose?
drugSummary <- within(drugSummary, 
                      dosePerDay.quartile <- do.call(rbind, (strsplit(as.character(drugSummary$quartiles.dosePerDay), 
                                                             ',', fixed=T))))
drugSummary <- (within(drugSummary, 
                       daysPassed.quartile <- do.call(rbind, (strsplit(as.character(drugSummary$quartiles.daysPassed), 
                                                              ',', fixed=T)))))

drugSummary$numberOfPurchases <- as.numeric(as.character(drugSummary$numberOfPurchases))
drugSummary$numberOfPeoplePurchase <- as.numeric(as.character(drugSummary$numberOfPeoplePurchase))

#plotting of drugs by drug families
for(i in 1:nrow(ATC.codes)){
  tempName.fam <- as.character(ATC.codes[i,2])
  temp.drugSummary <- drugSummary[grepl(paste("^", as.character(ATC.codes[i,1]), sep=""), drugSummary$ATC.code),]
  print(as.character(ATC.codes[i,1]))
  print(nrow(temp.drugSummary))
  
  pdf(paste(todaysDate, "_R6.SSRI.vnr.", tempName.fam, ".summary.pdf", sep=""))
  ##dose per day
  #quartiles split by 0,.25, .5, .75, .1
  a <- ggplot(temp.drugSummary, aes(x=reorder(drug.Name, as.numeric(numberOfPurchases)), y=as.numeric(dosePerDay.quartile[,3])) ) +
    geom_point(size=3, shape=21) +
    geom_errorbar(aes(ymin=as.numeric(dosePerDay.quartile[,2]), ymax=as.numeric(dosePerDay.quartile[,4]) ), width=.01, 
                  position = position_dodge2(width = 0.05, padding = 0.05)) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title=paste("drug summary for ", tempName.fam , ", dosePerDay", sep=""))
  
  ## days passed
  b <- ggplot(temp.drugSummary, aes(x=reorder(drug.Name, numberOfPurchases), y=as.numeric(daysPassed.quartile[,3])) ) +
    geom_point(size=3, shape=21) +
    geom_errorbar(aes(ymin=as.numeric(daysPassed.quartile[,2]), ymax=as.numeric(daysPassed.quartile[,4]) ), width=.01, 
                  position = position_dodge2(width = 0.05, padding = 0.05) ) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title=paste("drug summary for ",tempName.fam , "days passed", sep=""))
  
  print(a)
  print(b)
  dev.off()
}

write.table(drugSummary, paste(todaysDate, "R6_drugSummary.txt", sep="_"), sep="\t", quote=FALSE, row.names=F)

############################## 
###### Switching drugs? ###### 
############################## 

drug.switch <- as.data.frame(drugType.freq %>% group_by(FINNGENID) %>% count(drug.Family) )
temp.switch <- as.data.frame(table(drug.switch$FINNGENID))
nrow(temp.switch)
#[1] 40468
#R6 86004

#How many people have switched?
nrow(temp.switch[temp.switch$Freq>1,])
#[1] 14259
# R6 42897 

drugType.freq$switchStatus <- ifelse(as.character(drugType.freq$FINNGENID) %in% as.character(temp.switch[temp.switch$Freq>1,]$Var1), 
                                       1, 0)
drugType.freq$SSRIorNot <- ifelse(drugType.freq$drug.Family=="SSRI", 1, 0)
table(drugType.freq$SSRIorNot)

#0       1
#1188018 1160930

table(drugType.freq[drugType.freq$switchStatus==0,]$drug.Family)
#Lithium MAO.A.I     MRI  Others    SSRI
#  10189    3428   25605  134458  273107
#R6 
# MAO.A.I     MRI  Others    SSRI
#    4361   33003  185703  354132

table(drugType.freq[drugType.freq$switchStatus==1,]$drug.Family)
#Lithium MAO.A.I     MRI  Others    SSRI
#  46442   22688  116276  555954  580179
#R6  after age filter
# MAO.A.I     MRI  Others    SSRI
#   31322  158775  774854  806798

#who switched specifically from SSRI?
drugType.freq.ssri <- drugType.freq[drugType.freq$drug.Family=="SSRI",]
ssri.ids <- as.data.frame(table(drugType.freq.ssri$FINNGENID))
ssri.ids <- ssri.ids[order(ssri.ids$Var1),]
nrow(ssri.ids)
#[1] 49677
# R6 63310 

drugType.freq.notssri <- drugType.freq[drugType.freq$drug.Family!="SSRI",]
notssri.ids <- as.data.frame(table(drugType.freq.notssri$FINNGENID))
notssri.ids <- notssri.ids[order(notssri.ids$Var1),]
nrow(notssri.ids)
#[1] 48770, this value also includes people who are solely the others! 
#R6 62782

ssriAlone <- setdiff(ssri.ids$Var1, notssri.ids$Var1)
ssriSwitch <- intersect(ssri.ids$Var1, notssri.ids$Var1)
othersAlone <- setdiff(notssri.ids$Var1, ssri.ids$Var1)
length(ssriAlone)
#[1] 18691
#R6 23222 
length(ssriSwitch)
#[1] 30986
#R6 40088 
length(othersAlone)
#[1] 17784
#R6 22694 

18691 + 30986
#[1] 49677
23222 + 40088
#R6 63310

venn.diagram(
  x = list(ssri.ids$Var1, notssri.ids$Var1),
  category.names = c("SSRI" , "notSSRI"),
  filename = 'ssriIDsoverLap_SSRInoSSRI.png',
  output=TRUE
)


venn.diagram(
  x = list(na.omit(phenos.cov[phenos.cov$ANTIDEPRESSANTS==1,]$FINNGENID), na.omit(phenos.cov[phenos.cov$F5_DEPRESSIO==1,]$FINNGENID), na.omit(phenos.cov[phenos.cov$SSRI==1,]$FINNGENID) ),
  category.names = c("ANTIDEPRESSANTS","MDD" , "SSRI"),
  filename = 'MDD_SSRI_overlap.png',
  output=TRUE
)

venn.diagram(
  x = list(ssriAlone, ssriSwitch, othersAlone),
  category.names = c("SSRIalone" , "SSRIswitch", "othersAlone"),
  filename = 'ssriIDsoverLap_SSRIandOthers.png',
  output=TRUE
)

drugType.freq$SSRIswitchStatus <- ifelse(as.character(drugType.freq$FINNGENID) %in% as.character(ssriSwitch), 1, 
                                         ifelse(as.character(drugType.freq$FINNGENID) %in% as.character(ssriAlone), 0, 
                                                -9) )
# 0 = SSRI alone
# 1 = SSRI switched to something else
# -9 = Others alone
table(drugType.freq$SSRIswitchStatus)
#-9       0       1
#241652  273107 1253567
#R6
#     -9       0       1
# 299157  354132 1695659

#make drug families numerical
#drugType.freq$drug.Family.num <- ifelse(as.character(drugType.freq$drug.Family)=="SSRI", 1, 
#                                        ifelse(as.character(drugType.freq$drug.Family)=="Lithium", 2, 
 #                                              ifelse(as.character(drugType.freq$drug.Family)=="Others", 3, 
  #                                                    ifelse(as.character(drugType.freq$drug.Family)=="MAO.A.I", 4, 
   #                                                          ifelse(as.character(drugType.freq$drug.Family)=="MRI", 5,
    #                                                                -9 ) ) ) ) )

drugType.freq$drug.Family.num <- ifelse(as.character(drugType.freq$drug.Family)=="SSRI", 1, 
                                               ifelse(as.character(drugType.freq$drug.Family)=="Others", 2, 
                                                      ifelse(as.character(drugType.freq$drug.Family)=="MAO.A.I", 3, 
                                                             ifelse(as.character(drugType.freq$drug.Family)=="MRI", 4,
                                                                    -9 ) ) ) ) 

# 1 = SSRI
# 2 = Lithium
# 3 = Others
# 4 = Monoamine oxidase A inhibitors (MAO.A.I)
# 5 = Monoamine oxidase A inhibitors (MRI)
table(drugType.freq$drug.Family.num)

#1      2      3      4      5
#853286  56631 690412  26116 141881
#1       2       3       4       5
#1161182   70677  960663   35698  191853
#1       2       3       4
#1160930  960557   35683  191778



pdf(paste(todaysDate,"_R6_purchaseDist.SSRI.pdf", sep="") )

 ggplot(drugType.freq, aes(x=drug.Family)) + 
  geom_histogram(color="black", fill="white", stat="count") + 
  #geom_vline(aes(xintercept=median(SSRIfreq)), color="blue", linetype="dashed", size=1) + 
  theme_classic() + 
  labs(title="Antidepressant purchases")
dev.off()


#phenos.cov <- fread("zcat R5_COV_PHENO_V1.txt.gz", data.table=F)
#phenos.cov.drugs <- fread("zcat R5_COV_PHENO_V1_DRUGS_V1_ALL.txt.gz", data.table=F)
phenos.cov <- fread("zcat src/R6_COV_PHENO_V0.txt.gz", data.table=F)
#how about if I just obtain who switched or didn't
#based on ssriIDsoverLap_SSRInoSSRI.png
phenos.cov$SSRIswitchStatus <- ifelse(order(as.character(phenos.cov$FINNGENID) ) %in% order(as.character(ssriAlone) ), 0,
                                      ifelse(order(as.character(phenos.cov$FINNGENID) ) %in% order(as.character(ssriSwitch) ), 1, NA ) )

table(phenos.cov$SSRIswitchStatus)
#0     1
#18691 12295
#R6
#0     1
#23222 16866
#R6 after age filter
#    0     1
#18660 16438

#phenos.cov$SSRIswitchStatus <- ifelse( phenos.cov$AGE_AT_DEATH_OR_NOW > 35,
#                                       phenos.cov$SSRIswitchStatus, NA ) 
#table(phenos.cov$SSRIswitchStatus)
#    0     1
#16589 14574

#what was in R5?
phenos.cov.r5 <- fread("zcat ../r5/20200713_R5_COV_PHENO_KV.txt.gz", data.table=F)
r5.switchStatus <- as.data.frame(cbind(phenos.cov.r5$FINNGENID, phenos.cov.r5$SSRIswitchStatus) )
colnames(r5.switchStatus) <- c("FINNGENID", "SSRIswitchStatus.R5" )
phenos.cov <- phenos.cov[order(phenos.cov$FINNGENID),]
r5.switchStatus <- r5.switchStatus[order(r5.switchStatus$FINNGENID),]

phenos.cov <- merge(by.x="FINNGENID", by.y="FINNGENID", x=phenos.cov, y=r5.switchStatus, all=T)

table(phenos.cov$SSRIswitchStatus.R5)

#0     1
#18691 12295


#what were the na's from SSRI determined from Aki?
na.ids <- phenos.cov[is.na(phenos.cov$SSRI),]$FINNGENID
length(na.ids)
#[1] 22268


pdf(paste(todaysDate,"_R6_switch.SSRI.pdf", sep="") )

ggplot(phenos.cov[!is.na(phenos.cov$SSRIswitchStatus),], aes(x=SSRIswitchStatus)) + 
  geom_histogram(color="black", fill="white", stat="count") + 
  #geom_vline(aes(xintercept=median(SSRIfreq)), color="blue", linetype="dashed", size=1) + 
  theme_classic() + 
  labs(title="Antidepressant purchases")

ggplot(phenos.cov[!is.na(phenos.cov$SSRIswitchStatus.R5),], aes(x=SSRIswitchStatus.R5)) + 
  geom_histogram(color="black", fill="white", stat="count") + 
  #geom_vline(aes(xintercept=median(SSRIfreq)), color="blue", linetype="dashed", size=1) + 
  theme_classic() + 
  labs(title="Antidepressant purchases in R5")


dev.off()

### SSRI vs no SSRI
phenos.cov$SSRIstatus <- ifelse(order(as.character(phenos.cov$FINNGENID) ) %in% order(as.character(ssri.ids[ssri.ids$Freq >= 1,]$Var1) ), 1, 0)
table(phenos.cov$SSRIstatus)
#0      1
#197060  63345
#R6 
#     0      1
#198159  63310

#Aki's definition
table(phenos.cov$SSRI)
#0      1
#182776  51360

#repeat R5
r5.SSRI <- as.data.frame(cbind(phenos.cov.r5$FINNGENID, phenos.cov.r5$SSRI) )
colnames(r5.SSRI) <- c("FINNGENID", "SSRI.R5" )
phenos.cov <- phenos.cov[order(phenos.cov$FINNGENID),]
r5.SSRI <- r5.SSRI[order(r5.SSRI$FINNGENID),]

phenos.cov <- merge(by.x="FINNGENID", by.y="FINNGENID", x=phenos.cov, y=r5.SSRI, all=T)
table(phenos.cov$SSRI.R5)

#0      1
#159950  41147
why.na <- phenos.cov[is.na(phenos.cov$SSRI),]
why.na.ids <- why.na$FINNGENID
summary(why.na$AGE_AT_DEATH_OR_NOW)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#16.50   54.00   66.00   63.83   75.60  103.50    5065

drugType.freq <- drugType.freq[order(drugType.freq$FINNGENID),]
long.drug.na <- drugType.freq[as.character(drugType.freq$FINNGENID) %in% as.character(why.na.ids),]
nrow(long.drug.na)
#[1] 362703
freq.na <- as.data.frame(table(long.drug.na$FINNGENID) )
summary(freq.na$Freq)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1       3       7      17      19     314

pdf(paste(todaysDate,"_R6_SSRIstatus.pdf", sep="") )

ggplot(phenos.cov[!is.na(phenos.cov$SSRIstatus),], aes(x=SSRIstatus)) + 
  geom_histogram(color="black", fill="white", stat="count") + 
  #geom_vline(aes(xintercept=median(SSRIfreq)), color="blue", linetype="dashed", size=1) + 
  theme_classic() + 
  labs(title="SSRI vs no SSRI")

ggplot(phenos.cov[!is.na(phenos.cov$SSRI),], aes(x=SSRI)) + 
  geom_histogram(color="black", fill="white", stat="count") + 
  #geom_vline(aes(xintercept=median(SSRIfreq)), color="blue", linetype="dashed", size=1) + 
  theme_classic() + 
  labs(title="SSRI vs no SSRI (Aki's definition)")

ggplot(phenos.cov[!is.na(phenos.cov$SSRI.R5),], aes(x=SSRI.R5)) + 
  geom_histogram(color="black", fill="white", stat="count") + 
  #geom_vline(aes(xintercept=median(SSRIfreq)), color="blue", linetype="dashed", size=1) + 
  theme_classic() + 
  labs(title="SSRI vs no SSRI R5 (Aki's definition)")


dev.off()


todaysDate <- format.Date(Sys.Date(), "20%y%m%d")
#write.table(phenos.cov, paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="\t", quote = F, 
#            col.names = T, row.names = F)
write.table(phenos.cov, paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="\t", quote = F, 
            col.names = T, row.names = F)
#command<- paste("gzip -9 -v -f ", paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="")
command<- paste("gzip -9 -v -f ", paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="")
system(command)
save.image(file=paste(todaysDate,"RData", sep=".") )

#obtaining phenotype to test list for SAIGE
#phenos.list <- colnames(phenos.cov[,grep("^F5",colnames(phenos.cov) )] )
switchID <- as.data.frame(colnames(phenos.cov[grep("SSRIswitchStatus", colnames(phenos.cov))] ) )
colnames(switchID) <- "test"
SSRIstatus <- as.data.frame(colnames(phenos.cov[grep("SSRIstatus", colnames(phenos.cov))] ) )
colnames(SSRIstatus) <- "test"
SSRIstatus.Aki <- as.data.frame(colnames(phenos.cov[grep("SSRI", colnames(phenos.cov))] ) )
colnames(SSRIstatus.Aki) <- "test"
binaryStat <- as.data.frame(rbind(switchID, SSRIstatus, SSRIstatus.Aki))
#phenos.list <- c(phenos.list, as.character(switchID))  
write.table(binaryStat, paste(todaysDate, "phenoList.txt", sep="_"), sep="\t", quote=F, col.names = F, row.names = F)

############################## 
###### Purchase number? ###### 
############################## 

#total number of purchases
nrow(drugType.freq[drugType.freq$drug.Family=="SSRI",])
#853286
#R6 1161182 aafter age filter 1037840
drugType.freq.SSRI <- drugType.freq[drugType.freq$drug.Family=="SSRI",]
SSRI.num <- as.data.frame(table(drugType.freq$FINNGENID))

nrow(SSRI.num)
#[1] 67461
#R6 86004 
colnames(SSRI.num) <- c("FINNGENID", "SSRIfreq")
SSRI.num <- SSRI.num[order(as.character(SSRI.num$FINNGENID)),]
phenos.cov <- phenos.cov[order(as.character(phenos.cov$FINNGENID)),]
SSRI.num$FINNGENID <- as.character(SSRI.num$FINNGENID)
phenos.cov$FINNGENID <- as.character(phenos.cov$FINNGENID)

phenos.cov2 <- merge(by.x="FINNGENID", by.y="FINNGENID", 
                     x=phenos.cov[order(phenos.cov$FINNGENID),], y=SSRI.num[order(SSRI.num$FINNGENID),], 
                     all.x=T, sort=T)
table(is.na(phenos.cov2$SSRIfreq))
#FALSE   TRUE
#82835 178634

#I know there's an overlap!
#intersect(order(as.character(phenos.cov$FINNGENID)), order(as.character(SSRI.num$FINNGENID)))
test <- phenos.cov2[!is.na(phenos.cov2$SSRIfreq),]$FINNGENID
#Seem to be missing some?
nrow(SSRI.num)
#[1] 67461
#R6 86004 after age filter 74359
 length(test)
#[1] 65833 (1628 is missing?!?)
 #R6 83228 (3176 missing)
 missing <- SSRI.num[!(SSRI.num$FINNGENID %in% test),]$FINNGENID
 length(missing)
 #[1] 1628
 #R6 3169 
 
 #I know that the IDs exist!
 table(order(as.character(phenos.cov$FINNGENID)) %in% order(as.character(SSRI.num$FINNGENID)))
 # FALSE   TRUE
 #151331  67461
 
 #R6
 # FALSE   TRUE
 #175465  86004
 
 
 data.purchases.vnr.ssri <- data.purchases.vnr[grepl("SSRI", data.purchases.vnr$drug.Family),]
 ssriFreq <- as.data.frame(table(data.purchases.vnr.ssri$FINNGENID))
 nrow(ssriFreq)
 #[1] 56519 
 #R6 71860 
 
 colnames(ssriFreq) <- c("FINNGENID", "SSRIfreq")
 ssriFreq$FINNGENID <- as.character(ssriFreq$FINNGENID)
 ssriFreq <- ssriFreq[order(ssriFreq$FINNGENID),]
 phenos.cov <- phenos.cov[order(phenos.cov$FINNGENID),]
 phenos.cov <- merge(by.x="FINNGENID", by.y="FINNGENID", 
                      x=phenos.cov[order(phenos.cov$FINNGENID),], y=ssriFreq[order(ssriFreq$FINNGENID),], 
                      all.x=T, sort=T)
 table(is.na(phenos.cov$SSRIfreq))
 #R6
 # FALSE   TRUE
 # 69190 192279
 
 #found out that I need to save the NAs as 0 or 0.5 for the IRN (inverse rank normalized) transformation
 
phenos.cov$SSRIfreq <- ifelse(is.na(phenos.cov$SSRIfreq), 0.5, phenos.cov$SSRIfreq)
phenos.cov$SSRIfreq <- ifelse(phenos.cov$SSRIfreq==0, 0.5, phenos.cov$SSRIfreq)
phenos.cov$SSRIfreq <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35, phenos.cov$SSRIfreq, NA)
phenos.cov$SSRIfreq_Konrad <- ifelse(phenos.cov$SSRIfreq==0.5, NA, phenos.cov$SSRIfreq)

phenos.cov$SSRIfreq_IRN <- NA
phenos.cov[!is.na(phenos.cov$SSRIfreq),]$SSRIfreq_IRN <- 
  rankNorm(phenos.cov[!is.na(phenos.cov$SSRIfreq),]$SSRIfreq)

phenos.cov$SSRIfreq_Konrad_IRN <- NA
phenos.cov[!is.na(phenos.cov$SSRIfreq_Konrad),]$SSRIfreq_Konrad_IRN <-
  rankNorm(phenos.cov[!is.na(phenos.cov$SSRIfreq_Konrad),]$SSRIfreq_Konrad)


### What did Juha's dataset look like?
phenos.cov.drug.r5 <- fread("zcat ../r5/R5_COV_PHENO_V1_DRUGS_V1_ALL.txt.gz", data.table=F)
quantile(phenos.cov.drug.r5$ATC_N06AB, prob=c(0.1,0.25,0.5,0.75,0.9,1.0))
#10%   25%   50%   75%   90%  100%
#  0.5   0.5   0.5   1.0  12.0 302.0
r5.SSRIfreq <- as.data.frame(cbind(phenos.cov.drug.r5$FINNGENID, phenos.cov.drug.r5$ATC_N06AB, phenos.cov.drug.r5$ATC_N06AB_IRN )  )
colnames(r5.SSRIfreq) <- c("FINNGENID", "ATC_N06AB.R5", "ATC_N06AB_IRN.R5")
r5.SSRIfreq <- r5.SSRIfreq[order(r5.SSRIfreq$FINNGENID), ]
phenos.cov <- phenos.cov[order(phenos.cov$FINNGENID), ]
phenos.cov <- merge(by.x="FINNGENID", by.y="FINNGENID", x=phenos.cov, y=r5.SSRIfreq, all=T)

phenos.cov$SSRIfreq_IRN <- NA
phenos.cov$ATC_N06AB.R5 <- as.numeric(as.character(phenos.cov$ATC_N06AB.R5) )
phenos.cov$ATC_N06AB_IRN.R5 <- as.numeric(as.character(phenos.cov$ATC_N06AB_IRN.R5) )

phenos.cov$ATC_N06AB.R5_Konrad <-  ifelse(phenos.cov$ATC_N06AB.R5==0.5, NA, phenos.cov$ATC_N06AB.R5)
phenos.cov$ATC_N06AB.R5_Konrad_IRN <- NA
phenos.cov[!is.na(phenos.cov$ATC_N06AB.R5_Konrad),]$ATC_N06AB.R5_Konrad_IRN <-
  rankNorm(phenos.cov[!is.na(phenos.cov$ATC_N06AB.R5_Konrad),]$ATC_N06AB.R5_Konrad)


phenos.cov[!is.na(phenos.cov$ATC_N06AB.R5),]$SSRIfreq_IRN <- 
  rankNorm(phenos.cov[!is.na(phenos.cov$ATC_N06AB.R5),]$ATC_N06AB.R5)

quantile(phenos.cov$SSRIfreq_IRN, prob=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), na.rm = T)
quantile(phenos.cov$ATC_N06AB_IRN.R5, prob=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), na.rm = T)


 pdf(paste(todaysDate,"_purchaseDist.SSRI.pdf", sep="") )
 ggplot(phenos.cov[!is.na(phenos.cov$SSRIfreq),], aes(x=SSRIfreq)) + 
   geom_histogram(color="black", fill="white") + 
   geom_vline(aes(xintercept=median(SSRIfreq)), color="blue", linetype="dashed", size=1) + 
   theme_classic() + 
   labs(title="SSRI Purchase without IRN")
 
 ggplot(phenos.cov[!is.na(phenos.cov$SSRIfreq_IRN),], aes(x=SSRIfreq_IRN)) + 
   geom_histogram(color="black", fill="white") + 
   geom_vline(aes(xintercept=median(SSRIfreq_IRN)), color="blue", linetype="dashed", size=1) + 
   theme_classic() + 
   labs(title="SSRI Purchase with IRN")

 ggplot(phenos.cov[!is.na(phenos.cov$ATC_N06AB.R5),], aes(x=ATC_N06AB.R5)) + 
   geom_histogram(color="black", fill="white") + 
   geom_vline(aes(xintercept=median(ATC_N06AB.R5)), color="blue", linetype="dashed", size=1) + 
   theme_classic() + 
   labs(title="SSRI Purchase without IRN (Juha's)")
 
 ggplot(phenos.cov[!is.na(phenos.cov$ATC_N06AB_IRN.R5),], aes(x=ATC_N06AB_IRN.R5)) + 
   geom_histogram(color="black", fill="white") + 
   geom_vline(aes(xintercept=median(ATC_N06AB_IRN.R5)), color="blue", linetype="dashed", size=1) + 
   theme_classic() + 
   labs(title="SSRI Purchase with IRN (Juha's)")
 
 dev.off()
 
 
 todaysDate <- format.Date(Sys.Date(), "20%y%m%d")
 #write.table(phenos.cov, paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="\t", quote = F, 
  #           col.names = T, row.names = F)
 write.table(phenos.cov, paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="\t", quote = F, 
             col.names = T, row.names = F)
 
 command<- paste("gzip -9 -vf ", paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="")
 system(command)
 phenos.list <- colnames(phenos.cov[grep("ATC_N06AB_IRN.R5|SSRIfreq_IRN", colnames(phenos.cov))] )
 #phenos.list <- c(phenos.list, as.character(switchID))  
 write.table(phenos.list, paste(todaysDate, "phenoList_quantTrait.txt", sep="_"), sep="\t", quote=F, col.names = F, row.names = F)
 save.image(file=paste(todaysDate,"RData", sep=".") )
 
 
 ### FIGURE THIS OUT LATER
 missing <- ssriFreq[!(ssriFreq$FINNGENID %in% test),]$FINNGENID
 length(missing)
 
 
 ############################## 
 ##### Type of SSRI used ###### 
 ##############################
 
 table(data.purchases.vnr$drug.Family)
#Lithium MAO.A.I     MRI  Others    SSRI
#57370   27588  151900  714719  897845
 #R6 
 #Lithium MAO.A.I     MRI  Others    SSRI
 #71464   37529  204533  992450 1217716
 #R6 after age
 #  64788   35910  189269  891095 1085404
 
 #vnr subset for SSRI
 #data.purchases.vnr.ssri
 
 #VNR bursts and how many days between
 #data.purchases.vnr.bursts
 #additional information from data.purchases.vnr.bursts
 #drugType.freq
 
 #created my own map file to generic name binning for each commercial SSRI to create smaller bins source was Google
 ssriMap <- read.table("/home/kumar/workDir/20200528_SSRIgenericNamesMap.txt", header=T, sep="\t")
 #there are 72 drug.Family in  data.purchases.vnr.ssri
 #there are 72 in ssriMap
 ssriMap$DDDmg <- 0
 ssriMap$DDDmg <- ifelse(ssriMap$generic.Name=="Citalopram", 
                         20, ssriMap$DDD)
 ssriMap$DDDmg <- ifelse(ssriMap$generic.Name=="Escitalopram", 
                         10, ssriMap$DDD)
 ssriMap$DDDmg <- ifelse(ssriMap$generic.Name=="Fluoxetine", 
                         20, ssriMap$DDD)
 ssriMap$DDDmg <- ifelse(ssriMap$generic.Name=="Fluvoxamine", 
                         100, ssriMap$DDD)
 ssriMap$DDDmg <- ifelse(ssriMap$generic.Name=="Paroxetine", 
                         20, ssriMap$DDD)
 ssriMap$DDDmg <- ifelse(ssriMap$generic.Name=="Sertraline", 
                         50, ssriMap$DDD)
 ssriMap <- ssriMap[!grepl("^X", colnames(ssriMap))]
 
 ##################################
### using data.purchases.vnr.ssri ###
 ##################################
data.purchases.vnr.ssri <- merge(by.x="valmiste", by.y="drug.Name", x=data.purchases.vnr.ssri, y=ssriMap, all.x=T)
table(data.purchases.vnr.ssri$generic.Name)
#Citalopram Escitalopram   Fluoxetine  Fluvoxamine   Paroxetine   Sertraline
#330138       235856       134748        16216        66154       114733
#R6
#  Citalopram Escitalopram   Fluoxetine  Fluvoxamine   Paroxetine   Sertraline
#447764       308305       185617        23376        94187       158467
#after age filter
#      411221       262891       165932        21883        87017       136460
#DUH!! Could have also used this https://www.whocc.no/atc_ddd_index/?code=N06AB&showdescription=yes
#could have searched for full ATC code

#accidental inclusion of empty columns from excel
data.purchases.vnr.ssri <- data.purchases.vnr.ssri[!grepl("^X", colnames(data.purchases.vnr.ssri))]
data.purchases.vnr.ssri$vahvuus_num <- as.numeric(as.character(data.purchases.vnr.ssri$vahvuus_num))
summary(data.purchases.vnr.ssri[data.purchases.vnr.ssri$generic.Name=="Citalopram",]$vahvuus_num)
summary(data.purchases.vnr.ssri[data.purchases.vnr.ssri$generic.Name=="Escitalopram",]$vahvuus_num)
summary(data.purchases.vnr.ssri[data.purchases.vnr.ssri$generic.Name=="Fluoxetine",]$vahvuus_num)
summary(data.purchases.vnr.ssri[data.purchases.vnr.ssri$generic.Name=="Fluvoxamine",]$vahvuus_num)
summary(data.purchases.vnr.ssri[data.purchases.vnr.ssri$generic.Name=="Paroxetine",]$vahvuus_num)
summary(data.purchases.vnr.ssri[data.purchases.vnr.ssri$generic.Name=="Sertraline",]$vahvuus_num)


################################
### using drugType.freq.ssri ###
################################

drugType.freq.ssri <- merge(by.x="valmiste", by.y="drug.Name", x=drugType.freq.ssri, y=ssriMap, all.x=T)
table(drugType.freq.ssri$generic.Name)
#Citalopram Escitalopram   Fluoxetine  Fluvoxamine   Paroxetine   Sertraline
#313244       225342       125453        15236        63688       110323
#R6
# Citalopram Escitalopram   Fluoxetine  Fluvoxamine   Paroxetine   Sertraline
#426370       295302       173735        22128        90865       152782
### From without the first one
#Citalopram Escitalopram   Fluoxetine  Fluvoxamine   Paroxetine   Sertraline
#330138       235856       134748        16216        66154       114733
#after age filter
#      392448       253092       155636        20720        84041       131903
drugType.freq.ssri$dosePerDay <- as.numeric(as.character(drugType.freq.ssri$dosePerDay))


summary(drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Citalopram" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.02   14.20   20.62   37.47   37.04 5000.00
summary(drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Escitalopram" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0084    7.1770   10.7143   19.8105   19.3548 2000.0000
summary(drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Fluoxetine" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.043   12.903   20.000   33.466   32.787 6000.000
summary(drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Fluvoxamine" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.271   50.000   90.000  118.089  132.353 9000.000
summary(drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Paroxetine" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0628   14.7059   20.2020   32.7693   33.3333 2000.0000
summary(drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Sertraline" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.154    43.860    64.935   123.176   113.636 12500.000

pdf(paste(todaysDate,"_purchaseDistPerGeneric.SSRI.pdf", sep="") )
ggplot(drugType.freq.ssri[drugType.freq.ssri$daysPassed.diff!=0 & drugType.freq.ssri$dosePerDay < 500 ,], 
         aes(x=generic.Name, y=dosePerDay)) + 
  geom_boxplot(varwidth = TRUE, outlier.color = "red", outlier.shape = 1) +
  theme_classic() + 
  labs(title="SSRI Generic Drug per day dose", y = "mg/day", x="Generic Name") + 
  annotate("text", x = 3, y = 300, label = "Daily Drug Dosage\nCitalopram=20mg\nEscitalopram=10mg\nFluoxetine=20mg\nFluvoxamine=100mg\nParoxetine=20mg\nSertraline=50mg")

ggplot(drugType.freq.ssri[drugType.freq.ssri$daysPassed.diff!=0 & drugType.freq.ssri$dosePerDay < 500,], 
       aes(x=generic.Name, y=dosePerDay)) + 
  geom_boxplot(varwidth = TRUE, outlier.shape = NA) + 
  scale_y_continuous(limits = quantile(drugType.freq.ssri[drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay,
                                       c(0.1, 0.9))) +
  theme_classic() + 
  labs(title="SSRI Generic Drug per day dose", y = "mg/day", x="Generic Name") + 
  annotate("text", x = 3, y = 75, label = "Daily Drug Dosage\nCitalopram=20mg\nEscitalopram=10mg\nFluoxetine=20mg\nFluvoxamine=100mg\nParoxetine=20mg\nSertraline=50mg")

dev.off()

#MADddd = median(| DDD - median(x)|)
drugType.freq.ssri<-drugType.freq.ssri[order(drugType.freq.ssri$FINNGENID),]
toFilter <- drugType.freq.ssri[drugType.freq.ssri$FINNGENID %in% id.list,]

#hard statistical removal of outliers, removed lower 10th and upper 90th
quantile(toFilter[toFilter$daysPassed.diff!=0 & 
                    toFilter$generic.Name=="Sertraline",]$dosePerDay,  
         probs = c(0.1, 0.5, 0.75, 0.9))
##### SEPTEMBER 25th, 2020!!!!   !!!!!!!!VERIFY RANGE
##### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#daysPassed.diff> 15, daysPassed.diff<167
#Citalopram>7.5 Citalopram < 68.9
#Escitalopram > 3.1 Escitalopram < 37.3
#Fluoxetine>6.1 Fluoxetine <60
#Fluvoxamine > 28.67 Fluvoxamine < 214.8
#Paroxetine > 7.01 Paroxetine < 58.82
#Sertraline > 20.83 Sertraline < 217.4
drugType.freq.ssri.clean <- drugType.freq.ssri[drugType.freq.ssri$daysPassed.diff!=0 &
  drugType.freq.ssri$daysPassed.diff > 15 & drugType.freq.ssri$daysPassed.diff < 167 &
  (drugType.freq.ssri$generic.Name=="Sertraline" & drugType.freq.ssri$dosePerDay>20.83 & drugType.freq.ssri$dosePerDay < 217.4) | 
  (drugType.freq.ssri$generic.Name=="Paroxetine" & drugType.freq.ssri$dosePerDay>7.01 & drugType.freq.ssri$dosePerDay < 58.82) |
  (drugType.freq.ssri$generic.Name=="Fluvoxamine" & drugType.freq.ssri$dosePerDay>28.67 & drugType.freq.ssri$dosePerDay < 214.8) |
  (drugType.freq.ssri$generic.Name=="Fluoxetine" & drugType.freq.ssri$dosePerDay>6.1 & drugType.freq.ssri$dosePerDay < 60) |
  (drugType.freq.ssri$generic.Name=="Escitalopram" & drugType.freq.ssri$dosePerDay>3.1 & drugType.freq.ssri$dosePerDay < 37.3) |
  (drugType.freq.ssri$generic.Name=="Citalopram" & drugType.freq.ssri$dosePerDay>7.5 & drugType.freq.ssri$dosePerDay < 68.9),]
#to include the smaller days smaller
nrow(drugType.freq.ssri.clean)
#[1] 660352 from 853286
# R6 from 1161182 to 888866
#R6 age filter from 1037840 to 794670

pdf(paste(todaysDate,"_R6_purchaseDistPerGeneric.SSRI.clean.pdf", sep="") )
ggplot(drugType.freq.ssri.clean, 
       aes(x=generic.Name, y=dosePerDay)) + 
  geom_boxplot(varwidth = TRUE, outlier.color = "red", outlier.shape = 1) +
  theme_classic() + 
  labs(title="SSRI Generic Drug per day dose, clean for <10% > 90%", y = "mg/day", x="Generic Name") + 
  annotate("text", x = 3, y = 300, label = "Daily Drug Dosage\nCitalopram=20mg\nEscitalopram=10mg\nFluoxetine=20mg\nFluvoxamine=100mg\nParoxetine=20mg\nSertraline=50mg")

ggplot(drugType.freq.ssri.clean, 
       aes(x=generic.Name, y=dosePerDay)) + 
  geom_boxplot(varwidth = TRUE, outlier.shape = NA) + 
  scale_y_continuous(limits = quantile(drugType.freq.ssri[drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay,
                                       c(0.1, 0.9))) +
  theme_classic() + 
  labs(title="SSRI Generic Drug per day dose, clean for <10% > 90%", y = "mg/day", x="Generic Name") + 
  annotate("text", x = 3, y = 75, label = "Daily Drug Dosage\nCitalopram=20mg\nEscitalopram=10mg\nFluoxetine=20mg\nFluvoxamine=100mg\nParoxetine=20mg\nSertraline=50mg")

dev.off()

ssriDrugDose <- as.data.frame(drugType.freq.ssri.clean %>% 
                                group_by(FINNGENID, generic.Name) %>%
                                summarize(SSRImeanDose=mean(dosePerDay), SSRIsdDose=sd(dosePerDay), 
                                          SSRImedianDose=median(dosePerDay), SSRImadDDD=mad(dosePerDay), SSRInum=n(), 
                                          FirstDate=first(APPROX_EVENT_DAY), LastDate=last(APPROX_EVENT_DAY),
                                          difference= difftime(last(APPROX_EVENT_DAY), first(APPROX_EVENT_DAY), units="days"), numObsv=n()
                                          
                                          ) )
ssriDrugDose <- ssriDrugDose[order(ssriDrugDose$FirstDate, ssriDrugDose$FINNGENID),]
ssriDrugDose <- ssriDrugDose[!duplicated(ssriDrugDose$FINNGENID, fromLast=T),]
nrow(ssriDrugDose)
#45164
#R6  57520 after age filter 48759


#units per day
drugType.freq.ssri.clean$dosePerDayUnits <- drugType.freq.ssri.clean$dosePerDay / drugType.freq.ssri.clean$DDDmg

pdf(paste(todaysDate,"_R6_purchaseDistPerGeneric.SSRI.drugDoseUnit.clean.pdf", sep="") )
ggplot(drugType.freq.ssri.clean, 
       aes(x=generic.Name, y=dosePerDayUnits, fill=generic.Name)) + 
  geom_violin(varwidth = TRUE, outlier.color = "red", outlier.shape = 1) +
  geom_boxplot(width=0.1) + 
  theme_classic() + 
  labs(title="SSRI Generic Drug per day dose, clean for <10% > 90% dose units", y = "mg/day", x="Generic Name") + 
  annotate("text", x = 3, y=5, label = "Daily Drug Dosage\nCitalopram=20mg\nEscitalopram=10mg\nFluoxetine=20mg\nFluvoxamine=100mg\nParoxetine=20mg\nSertraline=50mg")

ggplot(drugType.freq.ssri.clean, 
       aes(x=generic.Name, y=dosePerDayUnits, fill=generic.Name)) + 
  geom_violin(varwidth = TRUE, outlier.color = "red", outlier.shape = 1) +
  geom_boxplot(width=0.1) + 
  scale_y_continuous(limits = quantile(drugType.freq.ssri.clean[drugType.freq.ssri.clean$daysPassed.diff!=0,]$dosePerDayUnits,
                                       c(0.1, 0.9))) +
  theme_classic() + 
  labs(title="SSRI Generic Drug per day dose, clean for <10% > 90%  dose units", y = "mg/day", x="Generic Name") + 
  annotate("text", x = 3, y=5, label = "Daily Drug Dosage\nCitalopram=20mg\nEscitalopram=10mg\nFluoxetine=20mg\nFluvoxamine=100mg\nParoxetine=20mg\nSertraline=50mg")

dev.off()

ssriDrugDose <- as.data.frame(drugType.freq.ssri.clean %>% 
                                group_by(FINNGENID, generic.Name) %>%
                                summarize(SSRImeanDose=mean(dosePerDayUnits), SSRIsdDose=sd(dosePerDayUnits), 
                                          SSRImedianDose=median(dosePerDayUnits), SSRImadDDD=mad(dosePerDayUnits), SSRInum=n(), 
                                          FirstDate=first(APPROX_EVENT_DAY), LastDate=last(APPROX_EVENT_DAY),
                                          difference= difftime(last(APPROX_EVENT_DAY), first(APPROX_EVENT_DAY), units="days"), numObsv=n()
                                          
                                ) )
ssriDrugDose <- ssriDrugDose[order(ssriDrugDose$FirstDate, ssriDrugDose$FINNGENID),]
#infrequently are people given >1 SSRI at a time
ssriDrugDose <- ssriDrugDose[!duplicated(ssriDrugDose$FINNGENID, fromLast=T),]



#I should take the median as opposed to the mean for dosage right? 
ssriDrugDose <- ssriDrugDose[c("FINNGENID", "generic.Name", "SSRImedianDose")]
#dupSSRI <- ssriDrugDose[duplicated(ssriDrugDose$FINNGENID),]

ssriDrugDose <- ssriDrugDose %>%
  gather(key, value, SSRImedianDose) %>%
  spread(generic.Name, value)

ssriDrugDose<- ssriDrugDose[,-2] #remove column key
### For per unit dosage
colnames(ssriDrugDose)<- c("FINNGENID", paste(grep("ine",colnames(ssriDrugDose), value=T ), "_units", sep="" ),
                           paste(grep("ram",colnames(ssriDrugDose), value=T ), "_units", sep="" ) )

nrow(ssriDrugDose)
#48509
# R6 [1] 57571 wth age filter 48759
ssriDrugDose <- ssriDrugDose[order(ssriDrugDose$FINNGENID),]
phenos.cov <- phenos.cov[order(phenos.cov$FINNGENID),]
phenos.cov <- merge(by.x="FINNGENID", by.y="FINNGENID", x=phenos.cov, y=ssriDrugDose, all.x=T)

#replace all NA and 0.0 with 0.5
#leave it for NAs
phenos.cov$Citalopram <- ifelse(is.na(phenos.cov$Citalopram) | phenos.cov$Citalopram < 0.5 , 0.5, phenos.cov$Citalopram)
phenos.cov$Citalopram <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35,  phenos.cov$Citalopram, NA)

phenos.cov$Escitalopram <- ifelse(is.na(phenos.cov$Escitalopram) | phenos.cov$Escitalopram < 0.5 , 0.5, phenos.cov$Escitalopram)
phenos.cov$Escitalopram <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35,  phenos.cov$Escitalopram, NA)

phenos.cov$Fluoxetine <- ifelse(is.na(phenos.cov$Fluoxetine) | phenos.cov$Fluoxetine < 0.5 , 0.5, phenos.cov$Fluoxetine)
phenos.cov$Fluoxetine <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35,  phenos.cov$Fluoxetine, NA)

phenos.cov$Fluvoxamine <- ifelse(is.na(phenos.cov$Fluvoxamine) | phenos.cov$Fluvoxamine < 0.5 , 0.5, phenos.cov$Fluvoxamine)
phenos.cov$Fluvoxamine <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35,  phenos.cov$Fluvoxamine, NA)

phenos.cov$Paroxetine <- ifelse(is.na(phenos.cov$Paroxetine) | phenos.cov$Paroxetine < 0.5 , 0.5, phenos.cov$Paroxetine)
phenos.cov$Paroxetine <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35,  phenos.cov$Paroxetine, NA)

phenos.cov$Sertraline <- ifelse(is.na(phenos.cov$Sertraline) | phenos.cov$Sertraline < 0.5 , 0.5, phenos.cov$Sertraline)
phenos.cov$Sertraline <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35,  phenos.cov$Sertraline, NA)


#IRN 
phenos.cov$Citalopram_IRN <- NA
phenos.cov$Escitalopram_IRN <- NA
phenos.cov$Fluoxetine_IRN <- NA
phenos.cov$Fluvoxamine_IRN <- NA
phenos.cov$Paroxetine_IRN <- NA
phenos.cov$Sertraline_IRN <- NA

phenos.cov[!is.na(phenos.cov$Citalopram),]$Citalopram_IRN <- rankNorm(phenos.cov[!is.na(phenos.cov$Citalopram),]$Citalopram)
phenos.cov[!is.na(phenos.cov$Escitalopram),]$Escitalopram_IRN <- rankNorm(phenos.cov[!is.na(phenos.cov$Escitalopram),]$Escitalopram)
phenos.cov[!is.na(phenos.cov$Fluoxetine),]$Fluoxetine_IRN <- rankNorm(phenos.cov[!is.na(phenos.cov$Fluoxetine),]$Fluoxetine)
phenos.cov[!is.na(phenos.cov$Fluvoxamine),]$Fluvoxamine_IRN <- rankNorm(phenos.cov[!is.na(phenos.cov$Fluvoxamine),]$Fluvoxamine)
phenos.cov[!is.na(phenos.cov$Paroxetine),]$Paroxetine_IRN <- rankNorm(phenos.cov[!is.na(phenos.cov$Paroxetine),]$Paroxetine)
phenos.cov[!is.na(phenos.cov$Sertraline),]$Sertraline_IRN <- rankNorm(phenos.cov[!is.na(phenos.cov$Sertraline),]$Sertraline)

#create matrix of this
#rankNorm vs Quantile norm (distribution of matrix equal)
#average of the first, average of the highest value
#keep the NAs as NAs
phenos.cov$SSRImeanDose_IRN <- 0

#matrix ind by drug quantile normalize that matrix

#phenos.cov$Citalopram_IRN <- ifelse(is.na(phenos.cov$Citalopram_IRN), 0, phenos.cov$Citalopram_IRN)
#phenos.cov$Escitalopram_IRN <- ifelse(is.na(phenos.cov$Escitalopram_IRN), 0, phenos.cov$Escitalopram_IRN)
#phenos.cov$Fluoxetine_IRN <- ifelse(is.na(phenos.cov$Fluoxetine_IRN), 0, phenos.cov$Fluoxetine_IRN)
#phenos.cov$Fluvoxamine_IRN <- ifelse(is.na(phenos.cov$Fluvoxamine_IRN), 0, phenos.cov$Fluvoxamine_IRN)
#phenos.cov$Paroxetine_IRN <- ifelse(is.na(phenos.cov$Paroxetine_IRN), 0, phenos.cov$Paroxetine_IRN)
#phenos.cov$Sertraline_IRN <- ifelse(is.na(phenos.cov$Sertraline_IRN), 0, phenos.cov$Sertraline_IRN)

phenos.cov$SSRImeanDose_IRN <- phenos.cov$Citalopram_IRN + phenos.cov$Escitalopram_IRN + phenos.cov$Fluoxetine_IRN + phenos.cov$Fluvoxamine_IRN + phenos.cov$Paroxetine_IRN + phenos.cov$Sertraline_IRN 

#phenos.cov$SSRImeanDose_IRN <- ifelse(phenos.cov$SSRImeanDose_IRN!=0, phenos.cov$SSRImeanDose_IRN, NA )
#following line previously commented out
#phenos.cov$SSRImeanDose_IRN <- ifelse(is.na(phenos.cov$SSRImeanDose_IRN),phenos.cov$SSRImeanDose_IRN, NA )

pdf(paste(todaysDate,"_R6_doseDist.SSRI.IRN.pdf", sep="") )
ggplot(phenos.cov, aes(x=SSRImeanDose_IRN)) + 
  geom_histogram(color="black", fill="white") + 
  #geom_vline(aes(xintercept=median(SSRImeanDose_IRN)), color="blue", linetype="dashed", size=1) + 
  theme_classic() + 
  labs(title="SSRI Dose careful IRN")

dev.off()

todaysDate <- format.Date(Sys.Date(), "20%y%m%d")
write.table(phenos.cov, paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="\t", quote = F, 
            col.names = T, row.names = F)
command<- paste("gzip -9 -v -f ", paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="")
system(command)
phenos.list <- colnames(phenos.cov[grep("SSRImeanDose_", colnames(phenos.cov))] )
phenos.list <- c(phenos.list,colnames(phenos.cov[grep("SSRIfreq_IRN", colnames(phenos.cov))] ) ) #repeat freq data


write.table(phenos.list, paste(todaysDate, "phenoList.txt", sep="_"), sep="\t", quote=F, col.names = F, row.names = F)
save.image(file=paste(todaysDate,"RData", sep=".") )



#per unit nomalization
#each SSRI is already normalized but  I need to create a new column 
phenos.cov$SSRImedianUnits <- 0.5
phenos.cov$Citalopram_units <- ifelse(is.na(phenos.cov$Citalopram_units) | phenos.cov$Citalopram_units < 0.5, 
                                      0.5, phenos.cov$Citalopram_units)
phenos.cov$Citalopram_units <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35,  phenos.cov$Citalopram_units, NA)

phenos.cov$Escitalopram_units <- ifelse(is.na(phenos.cov$Escitalopram_units) | phenos.cov$Escitalopram_units < 0.5, 
                                      0.5, phenos.cov$Escitalopram_units)
phenos.cov$Escitalopram_units <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35,  phenos.cov$Escitalopram_units, NA)

phenos.cov$Fluoxetine_units <- ifelse(is.na(phenos.cov$Fluoxetine_units) | phenos.cov$Fluoxetine_units < 0.5, 
                                        0.5, phenos.cov$Fluoxetine_units)
phenos.cov$Fluoxetine_units <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35,  phenos.cov$Fluoxetine_units, NA)

phenos.cov$Fluvoxamine_units <- ifelse(is.na(phenos.cov$Fluvoxamine_units) | phenos.cov$Fluvoxamine_units < 0.5, 
                                      0.5, phenos.cov$Fluvoxamine_units)
phenos.cov$Fluvoxamine_units <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35,  phenos.cov$Fluvoxamine_units, NA)

phenos.cov$Paroxetine_units <- ifelse(is.na(phenos.cov$Paroxetine_units) | phenos.cov$Paroxetine_units < 0.5, 
                                       0.5, phenos.cov$Paroxetine_units)
phenos.cov$Paroxetine_units <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35,  phenos.cov$Paroxetine_units, NA)

phenos.cov$Sertraline_units <- ifelse(is.na(phenos.cov$Sertraline_units) | phenos.cov$Sertraline_units < 0.5, 
                                      0.5, phenos.cov$Sertraline_units)
phenos.cov$Sertraline_units <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35,  phenos.cov$Sertraline_units, NA)


drugCols <- grep("units", colnames(phenos.cov), value=T)
pdf(paste(todaysDate, "SSRIunits.perDrug.preNorm.pdf") ) 
par(mfrow=c(3,1))
for(i in 1:length(drugCols)){
  print(drugCols[i])
  rowsCounted <- nrow(ssriDrugDose[!is.na(ssriDrugDose[,drugCols[i]]),])
  a<- ggplot(ssriDrugDose, aes_string(x = drugCols[i])  ) + 
    geom_histogram(fill="white", color="black") + 
    theme_classic() + 
    labs(title=paste("SSRI Dose Unit per Drug ", drugCols[i], "\nNumber of Individuals = ", rowsCounted, sep="" ))
  print(a)
}
dev.off()


phenos.cov$SSRImedianUnits <- 0
phenos.cov$SSRImedianUnits <-  as.numeric(as.character(phenos.cov$Citalopram_units)) +
  as.numeric(as.character(phenos.cov$Escitalopram_units)) + as.numeric(as.character(phenos.cov$Fluoxetine_units)) +
  as.numeric(as.character(phenos.cov$Fluvoxamine_units)) + as.numeric(as.character(phenos.cov$Paroxetine_units)) + 
  as.numeric(as.character(phenos.cov$Sertraline_units))
phenos.cov$SSRImedianUnits <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35, 
                                     phenos.cov$SSRImedianUnits, NA)
phenos.cov$SSRImedianUnits_IRN <- NA
phenos.cov[!is.na(phenos.cov$SSRImedianUnits),]$SSRImedianUnits_IRN<- rankNorm(phenos.cov[!is.na(phenos.cov$SSRImedianUnits),]$SSRImedianUnits)
phenos.cov$SSRImedianUnits_IRN <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35, 
                                     phenos.cov$SSRImedianUnits_IRN, NA)

summary(phenos.cov$SSRImedianUnits)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#0.26    0.88    1.03    1.23    1.50    4.26  174673
#R6
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#0.26    0.87    1.03    1.22    1.49    4.26  204888
#R6 after age
#   0.29    0.88    1.02    1.22    1.48    4.35  207645

summary(phenos.cov$SSRImedianUnits_IRN)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#-4.19   -0.67    0.00    0.00    0.67    3.96  174673
#R6
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#-4.24   -0.67    0.00    0.00    0.67    4.02  204888
#R6 after age filter
#  -4.07   -0.67   -0.01    0.00    0.67    3.79  207645


pdf(paste(todaysDate,"_doseUnitDist.SSRI.pdf", sep="") )
ggplot(phenos.cov, aes(x=SSRImedianUnits)) + 
  geom_histogram(color="black", fill="white") + 
  #geom_vline(aes(xintercept=median(SSRImeanDose_IRN)), color="blue", linetype="dashed", size=1) + 
  theme_classic() + 
  labs(title="SSRI Dose Unit")

ggplot(phenos.cov, aes(x=SSRImedianUnits_IRN)) + 
  geom_histogram(color="black", fill="white") + 
  #geom_vline(aes(xintercept=median(SSRImeanDose_IRN)), color="blue", linetype="dashed", size=1) + 
  theme_classic() + 
  labs(title="SSRI Dose Unit with IRN")


dev.off()


pdf("mitjaPlot.pdf")
#plot(phenos.cov$SSRImedianUnits_IRN, phenos.cov$SSRImeanDose_IRN)

#plot(phenos.cov$Paroxetine_IRN, phenos.cov$Paroxetine_units)
IRN.unit <- rankNorm(phenos.cov[!is.na(phenos.cov$Paroxetine_units),]$Paroxetine_units) 
IRN.median <- phenos.cov[!is.na(phenos.cov$Paroxetine_IRN),]$Paroxetine_IRN
plot(as.numeric(as.character(IRN.median)), as.numeric(as.character(IRN.unit)) )

dev.off()
#plot(phenos.cov$Citalopram_IRN, phenos.cov$Citalopram_units)
IRN.drug <- rankNorm(phenos.cov[!is.na(phenos.cov$Citalopram_units),]$Citalopram_units) 
plot(phenos.cov$Citalopram_IRN, IRN.drug )

#plot(phenos.cov$Escitalopram_IRN, phenos.cov$Escitalopram_units)
IRN.drug <- rankNorm(phenos.cov[!is.na(phenos.cov$Escitalopram_units),]$Escitalopram_units) 
plot(phenos.cov$Escitalopram_IRN,  IRN.drug)

#plot(phenos.cov$Fluoxetine_IRN, phenos.cov$Fluoxetine_units)
IRN.drug <- rankNorm(phenos.cov[!is.na(phenos.cov$Fluoxetine_units),]$Fluoxetine_units) 
plot(phenos.cov$Fluoxetine_IRN, IRN.drug )

#plot(phenos.cov$Fluvoxamine_IRN, phenos.cov$Fluvoxamine_units)
IRN.drug <- rankNorm(phenos.cov[!is.na(phenos.cov$Fluvoxamine_units),]$Fluvoxamine_units) 
plot(phenos.cov$Fluvoxamine_IRN, IRN.drug )

#plot(phenos.cov$Sertraline_IRN, phenos.cov$Sertraline_units)
IRN.drug <- rankNorm(phenos.cov[!is.na(phenos.cov$Sertraline_IRN),]$Sertraline_IRN) 
plot(phenos.cov$Sertraline_IRN, IRN.drug )

dev.off()

v1_irn <- rankNorm(v1)
v2 <- v1 / median(v1)
v2_irn <- rankNorm(v2)
plot(v1_irn,v2_irn)


todaysDate <- format.Date(Sys.Date(), "20%y%m%d")
write.table(phenos.cov, paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="\t", quote = F, 
            col.names = T, row.names = F)
command<- paste("gzip -9 -v -f ", paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="")
system(command)
phenos.list <- colnames(phenos.cov[grep("IRN", colnames(phenos.cov))] )
#meanDose <- colnames(phenos.cov[grep("SSRImeanDose_IRN", colnames(phenos.cov))] )
#freq <- colnames(phenos.cov[grep("SSRIfreq_IRN", colnames(phenos.cov))] )

write.table(phenos.list, paste(todaysDate, "phenoList_quant.txt", sep="_"), sep="\t", quote=F, col.names = F, row.names = F)
save.image(file=paste(todaysDate,"RData", sep=".") )

 ############################## 
 ###### Polygenic Risk S ###### 
 ##############################
 
 #gs://r5_data/PRS/scores/MDD2018_ex23andMe.19fields.sscore
mdd.prs <- read.table("~/workDir/MDD2018_ex23andMe.19fields.sscore", header=T, sep="\t")
mdd.prs <- mdd.prs[,c(1,5)] #finngenid and scores only
colnames(mdd.prs) <- c("FINNGENID", "MDD_PRS")
mdd.prs$MDD_PRS.norm<- ((mdd.prs$MDD_PRS)-mean(mdd.prs$MDD_PRS)/sd(mdd.prs$MDD_PRS))
phenos.cov <- phenos.cov[order(phenos.cov$FINNGENID),]
mdd.prs <- mdd.prs[order(mdd.prs$FINNGENID),]

phenos.cov <- merge(by.x="FINNGENID", by.y="FINNGENID", 
                    x=phenos.cov[order(phenos.cov$FINNGENID),], y=mdd.prs[order(mdd.prs$FINNGENID),], 
                    all.x=T)
#phenos.cov$MDD_PRS <- ifelse(phenos.cov$AGE_AT_DEATH_OR_NOW > 35, phenos.cov$MDD_PRS, NA )

todaysDate <- format.Date(Sys.Date(), "20%y%m%d")
pdf(paste(todaysDate, "_R6_PRS" ,".pdf", sep="") ) 
#plot.list <- c(grep("F5_", colnames(phenos.cov) , value=T), "SSRI")
plot.list <- c(grep("SSRI", colnames(phenos.cov) , value=T))
plot.list <- plot.list[c(1:5, 11,13,14)]
columnClasses.tTest <- c("character", 
                         "numeric",
                         "numeric",
                         "numeric",
                         "numeric"
)
columnNames.tTest <- c("pheno",
                       "stat", 
                       "pVal",
                       "Upper",
                       "Lower"
)

columnNames.table.test.likepTDT <- c("plotName",
                                     "t.test.TDTlike", 
                                     "df", "estimate",
                                     "mean.cases",
                                     "mean.controls",
                                     "lower", 
                                     "upper",
                                     "p.value")
columnClasses.table.test.likepTDT <- c("character",
                                       "numeric", 
                                       "numeric", "numeric",
                                       "numeric", "numeric",
                                       "numeric", 
                                       "numeric",
                                       "numeric")

tTestResults <- read.table(text = "",
                           colClasses = columnClasses.tTest,
                           col.names = columnNames.tTest)


table.test.likepTDT <- read.table(text = "",
                                  col.names = columnNames.table.test.likepTDT,
                                  colClasses = columnClasses.table.test.likepTDT)

for(i in 1:length(plot.list)){
  plot.name <- as.character(plot.list[i])
  print(paste("plotting", plot.name))
  forPlotting <- cbind(phenos.cov[grep(plot.name, colnames(phenos.cov))], 
                       phenos.cov[grep("MDD_PRS", colnames(phenos.cov))] ) 
  
  a <- ggplot(forPlotting[!is.na(forPlotting[,1]),], aes(x = MDD_PRS, fill=as.factor(get(plot.name) ) ) ) +
    geom_density(alpha = .3) +
    theme_classic() + 
    labs(title=paste(todaysDate, "PRS : ", plot.name,  " from MDD summstats", sep="") )
  
  b <- ggplot(forPlotting[!is.na(forPlotting[,1]),], aes(y = MDD_PRS, x=as.factor(get(plot.name))) ) +
    geom_boxplot() +
    labs(title=paste(todaysDate, "PRS : ", plot.name,  " from MDD summstats", sep="") ) +
    theme_classic()
   print(a)
  print(b)
  
  pointStat <- t.test(forPlotting[forPlotting[,1]==0 & !is.na(forPlotting[,1]),]$MDD_PRS, forPlotting[forPlotting[,1]==1 & !is.na(forPlotting[,1]),]$MDD_PRS)$estimate[1] -
    t.test(forPlotting[forPlotting[,1]==0 & !is.na(forPlotting[,1]),]$MDD_PRS, forPlotting[forPlotting[,1]==1 & !is.na(forPlotting[,1]),]$MDD_PRS)$estimate[2]
  sigT <- t.test(forPlotting[forPlotting[,1]==0 & !is.na(forPlotting[,1]),]$MDD_PRS, forPlotting[forPlotting[,1]==1 & !is.na(forPlotting[,1]),]$MDD_PRS)$p.value
  Upper <- t.test(forPlotting[forPlotting[,1]==0 & !is.na(forPlotting[,1]),]$MDD_PRS, forPlotting[forPlotting[,1]==1 & !is.na(forPlotting[,1]),]$MDD_PRS)$conf.int[1]
  Lower <- t.test(forPlotting[forPlotting[,1]==0 & !is.na(forPlotting[,1]),]$MDD_PRS, forPlotting[forPlotting[,1]==1 & !is.na(forPlotting[,1]),]$MDD_PRS)$conf.int[2]
  
  tTestResults <- rbind(tTestResults, cbind(plot.name, pointStat, sigT, Upper, Lower))
  testDeduct <- (na.omit(forPlotting[forPlotting[,1]==1,]$MDD_PRS)  - na.omit(forPlotting[forPlotting[,1]==0,]$MDD_PRS) )
  deviation <- (na.omit(forPlotting[forPlotting[,1]==1,]$MDD_PRS)  - na.omit(forPlotting[forPlotting[,1]==0,]$MDD_PRS) ) /
    sd(na.omit(forPlotting[forPlotting[,1]==1,]$MDD_PRS ) )
  trial <- t.test(deviation, mu=0 )
  
  table.test.likepTDT <- rbind(table.test.likepTDT,   
                               cbind(plot.name, trial$statistic, trial$parameter, trial$estimate,
                                     mean(na.omit(forPlotting[forPlotting[,1]==1,]$MDD_PRS) ),
                                     mean(na.omit(forPlotting[forPlotting[,1]==0,]$MDD_PRS) ),
                                     trial$conf.int[1], trial$conf.int[2], as.numeric(trial$p.value)) )
  
  }
dev.off()

colnames(table.test.likepTDT) <- c("plotName",
                                   "t.test.TDTlike", 
                                   "df", "estimate",
                                   "mean.cases",
                                   "mean.controls",
                                   "lower", 
                                   "upper",
                                   "p.value")


write.table(tTestResults, paste(todaysDate, "MDD.PRS.tTestResults.txt", sep="_"), sep="\t", quote=F, col.names = T, row.names = F)
write.table(table.test.likepTDT, paste(todaysDate, "MDD.PRS.table.test.likepTDT.txt", sep="_"), sep="\t", quote=F, col.names = T, row.names = F)

tTestResults$pointStat <- as.numeric(as.character(tTestResults$pointStat))
tTestResults$Upper <- as.numeric(as.character(tTestResults$Upper))
tTestResults$Lower <- as.numeric(as.character(tTestResults$Lower))

table.test.likepTDT$estimate <- as.numeric(as.character(table.test.likepTDT$estimate))
table.test.likepTDT$upper <- as.numeric(as.character(table.test.likepTDT$upper))
table.test.likepTDT$lower <- as.numeric(as.character(table.test.likepTDT$lower))
table.test.likepTDT$p.value <- as.numeric(as.character(table.test.likepTDT$p.value))
table.test.likepTDT$mean.controls <- as.numeric(as.character(table.test.likepTDT$mean.controls))
table.test.likepTDT$mean.cases <- as.numeric(as.character(table.test.likepTDT$mean.cases))

#head(table.test.likepTDT[order(-table.test.likepTDT$estimate) & table.test.likepTDT$p.value<0.05,])
#summary(table.test.likepTDT$mean.cases - table.test.likepTDT$mean.controls)
table.test.likepTDT$mean.diff <- table.test.likepTDT$mean.cases - table.test.likepTDT$mean.controls


pdf(paste(todaysDate, "MDD.PRS.tTestResults.pdf", sep="_"))
ggplot(tTestResults,aes(x=plot.name, y=pointStat) ) +
  geom_point(shape=21, size=1.0, fill="white") + 
  geom_errorbar(aes(ymin=Lower,ymax=Upper,width=1.2)) +
  theme_classic() +
  theme(axis.title.x = element_text(margin = unit(c(3, 2, 0, 0), "mm"))) + 
  geom_hline(yintercept = 0.0, colour="green", linetype = "longdash") +
  labs(title="2 sample t-test", y="t-statistic", x="outcomes/covariates")+
  coord_flip()

ggplot(table.test.likepTDT,aes(x=plotName, y=mean.diff) ) +
  geom_point(shape=21, size=1.0, fill="white") + 
  geom_errorbar(aes(ymin=lower,ymax=upper,width=1.2)) +
  theme_classic() +
  theme(axis.title.x = element_text(margin = unit(c(3, 2, 0, 0), "mm"))) + 
  geom_hline(yintercept = 0.0, colour="green", linetype = "longdash") + 
  labs(title="Deviation from the case status mean", y="mean.diff deviation", x="outcomes/covariates")+
  coord_flip()

dev.off()

todaysDate <- format.Date(Sys.Date(), "20%y%m%d")
write.table(phenos.cov, paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="\t", quote = F, 
            col.names = T, row.names = F)
command<- paste("gzip -9 -v -f ", paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="")
system(command)
save.image(paste(todaysDate, ".RData", sep=""))

###################################
##### PRS with Pete's script ######
###################################
library(GENESIS)
kinship.mat <- read.table("PRS/kinship/finngen_R5.kin0", header=T, sep="\t")


###################################
##### Phenotype correlations ######
###################################

table(phenos.cov$F5_DEPRESSIO)
#0      1
#225094  27855

#SSRI v MDD
table(phenos.cov$SSRI)
#0      1
#182776  51360


#ANTIDEPRESSANTS
cor.test(phenos.cov$ANTIDEPRESSANTS, phenos.cov$F5_DEPRESSIO)
#0.4438992
length(intersect(phenos.cov[phenos.cov$ANTIDEPRESSANTS==1,]$FINNGENID, phenos.cov[phenos.cov$F5_DEPRESSIO ==1,]$FINNGENID ) ) /
  length(union(phenos.cov[phenos.cov$ANTIDEPRESSANTS==1,]$FINNGENID, phenos.cov[phenos.cov$F5_DEPRESSIO ==1,]$FINNGENID ))
24278/ 77207

length(intersect(phenos.cov[phenos.cov$SSRI==1,]$FINNGENID, phenos.cov[phenos.cov$F5_DEPRESSIO ==1,]$FINNGENID ) )
#[1] 19916

length(union(phenos.cov[phenos.cov$SSRI==1,]$FINNGENID, phenos.cov[phenos.cov$F5_DEPRESSIO ==1,]$FINNGENID ))
# 59301

19916 / 59301
#0.335

#SSRI adherence v MDD 
cor.test(phenos.cov$SSRI, phenos.cov$F5_DEPRESSIO)
#0.5183565
length(intersect(phenos.cov[phenos.cov$SSRIswitchStatus==1,]$FINNGENID, phenos.cov[phenos.cov$F5_DEPRESSIO ==1,]$FINNGENID ) )
length(union(phenos.cov[phenos.cov$SSRIswitchStatus==1,]$FINNGENID, phenos.cov[phenos.cov$F5_DEPRESSIO ==1,]$FINNGENID ))
1813 / 42910
# 0.04225122

#SSRI dosage v MDD
cor.test(phenos.cov$SSRI_increase.surv, phenos.cov$F5_DEPRESSIO)
#0.007031051
length(intersect(phenos.cov[phenos.cov$SSRI_increase.surv==1,]$FINNGENID, phenos.cov[phenos.cov$F5_DEPRESSIO ==1,]$FINNGENID ) )
length(union(phenos.cov[phenos.cov$SSRI_increase.surv==1,]$FINNGENID, phenos.cov[phenos.cov$F5_DEPRESSIO ==1,]$FINNGENID ))
10689/48762
#0.2192076

####################################
######### Survival Analysis ########
####################################

todaysDate <- format.Date(Sys.Date(), "20%y%m%d")
table(data.purchases.vnr$drug.Family)
#Lithium MAO.A.I     MRI  Others    SSRI
#57370   27588  151900  714719  897845

#R6
#Lithium MAO.A.I     MRI  Others    SSRI
#71464   37529  204533  992450 1217716

#vnr subset for SSRI
#data.purchases.vnr.ssri
table(data.purchases.vnr.ssri$drug.Family)
# SSRI
#897845
# R6 1217716

data.purchases.vnr.ssri <- data.purchases.vnr.ssri[order(data.purchases.vnr.ssri$APPROX_EVENT_DAY), ] 
#### Dates
survival.data <- 
  as.data.frame(data.purchases.vnr.ssri %>%
                  group_by(FINNGENID, generic.Name) %>%
                  summarise(FirstDate=first(APPROX_EVENT_DAY), LastDate=last(APPROX_EVENT_DAY),
                 difference= difftime(last(APPROX_EVENT_DAY), first(APPROX_EVENT_DAY), units="days"), numObsv=n()
                                   )
                 )  

nrow(survival.data)
#85318
#R6 110731
  
survival.data.drugType <- 
  as.data.frame(drugType.freq.ssri %>%
                  group_by(FINNGENID) %>%
                  summarise(FirstDate=first(APPROX_EVENT_DAY), LastDate=last(APPROX_EVENT_DAY),
                            difference= difftime(last(APPROX_EVENT_DAY), first(APPROX_EVENT_DAY), units="days")
                  )
  )  
#diff between the two ssri objects above = 6842
#R6 63345

sun.data <- 
  as.data.frame(survival.data %>% 
                  group_by(FINNGENID) %>% 
                  mutate(drugFlow = paste0(generic.Name, collapse = "-"))  )
sun.data.table <- as.data.frame(table(sun.data$drugFlow))
write.table(sun.data.table, "R6_SSRI_sunburstFlow.txt", quote = F, sep="\t", row.names = F, col.names = T)

setwd("~/Dropbox (Partners HealthCare)/Daly_Lab/drugResponseAnalysis/")
sun.data.table <- read.table("R6_SSRI_sunburstFlow.txt", sep="\t", header=T)
sunburst(sun.data.table)
#https://stackoverflow.com/questions/39921971/create-sunburst-plot-in-shiny-using-html-instead-of-sunburstoutput
htmlwidgets::onRender(
  sunburst(sun.data.table),
  '
  function(el,x){
  d3.select(el).select(".sunburst-sidebar").remove()
  }
  '
)

#### Generic Drug switching #####
SSRIswitchGeneric.ID <- survival.data[duplicated(survival.data$FINNGENID, fromLast=T),]$FINNGENID
#28799
#R6 38871

SSRIswitchGeneric.ID.didnt <- survival.data[!duplicated(survival.data$FINNGENID, fromLast=T),]$FINNGENID
SSRIswitchGeneric.ID <- SSRIswitchGeneric.ID[order(SSRIswitchGeneric.ID)]
SSRIswitchGeneric.ID.didnt <- SSRIswitchGeneric.ID.didnt[order(SSRIswitchGeneric.ID.didnt)]
didntSwitch <- setdiff(SSRIswitchGeneric.ID.didnt, unique(SSRIswitchGeneric.ID) )
switched <- intersect(SSRIswitchGeneric.ID.didnt, unique(SSRIswitchGeneric.ID) )
phenos.cov$SSRIswitchStatus_generic <- 3
phenos.cov <- phenos.cov[order(phenos.cov$FINNGENID),]
switched <- switched[order(switched)]
didntSwitch <- didntSwitch[order(didntSwitch)]
#36089
#R6 44924

phenos.cov$SSRIswitchStatus_generic <- ifelse(phenos.cov$FINNGENID %in% as.character(switched), 1,
                                              ifelse(phenos.cov$FINNGENID %in% as.character(didntSwitch), 0, NA))


todaysDate <- format.Date(Sys.Date(), "20%y%m%d")
write.table(phenos.cov, paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="\t", quote = F, 
            col.names = T, row.names = F)
command<- paste("gzip -9 -v -f ", paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="")
system(command)
phenos.list <- "SSRIswitchStatus_generic"
phenos.list
write.table(phenos.list, paste(todaysDate, "phenoList.txt", sep="_"), sep="\t", quote=F, col.names = F, row.names = F)
save.image(file=paste(todaysDate,"RData", sep=".") )


#figuring out duplicates and removing later date ones
survival.data <- survival.data[order(survival.data$FirstDate, survival.data$FINNGENID),]
survival.data <- survival.data[!duplicated(survival.data$FINNGENID, fromLast=T),]

nrow(survival.data[duplicated(survival.data$FINNGENID, fromLast=T),])

survival.data$difference.days <- as.numeric(survival.data$difference)
summary(survival.data$difference.days)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0      18     285    1176    1590    8762
#R6
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0      17     289    1197    1589    9126

#nrow(survival.data.test[survival.data.test$numObsv<2,])
#[1] 8761

#survival.data <- survival.data[survival.data$numObsv>=2,]
#dosage difference
survival.data <- survival.data[order(survival.data$FirstDate,survival.data$FINNGENID, decreasing = T),]
survival.data$FirstDate.ID <- paste(survival.data$FINNGENID, survival.data$FirstDate, sep="_")
survival.data <- survival.data[order(survival.data$LastDate,survival.data$FINNGENID, decreasing = T),]
survival.data$LastDate.ID <- paste(survival.data$FINNGENID, survival.data$LastDate, sep="_")

drugType.freq.ssri.clean <- drugType.freq.ssri.clean[order(drugType.freq.ssri.clean$APPROX_EVENT_DAY, drugType.freq.ssri.clean$FINNGENID, decreasing = T),]
drugType.freq.ssri.clean$Data.ID <- paste(drugType.freq.ssri.clean$FINNGENID, drugType.freq.ssri.clean$APPROX_EVENT_DAY, sep="_")
drugType.freq.ssri.clean$dosePerDayUnits <- drugType.freq.ssri.clean$dosePerDay / drugType.freq.ssri.clean$DDDmg

FirstDate.Dose <- as.data.frame(cbind(drugType.freq.ssri.clean[drugType.freq.ssri.clean$Data.ID %in% survival.data$FirstDate.ID,]$dosePerDayUnits, drugType.freq.ssri.clean[drugType.freq.ssri.clean$Data.ID %in% survival.data$FirstDate.ID,]$FINNGENID) )
colnames(FirstDate.Dose) <- c("dosePerDayUnits.First", "FINNGENID")
FirstDate.Dose$dosePerDayUnits.First <- as.numeric(as.character(FirstDate.Dose$dosePerDayUnits.First))
FirstDate.Dose$dosePerDayUnits.First <- ifelse(FirstDate.Dose$dosePerDayUnits.First=="Inf", -99, FirstDate.Dose$dosePerDayUnits.First)

LastDate.Dose <- as.data.frame(cbind(drugType.freq.ssri.clean[drugType.freq.ssri.clean$Data.ID %in% survival.data$LastDate.ID,]$dosePerDayUnits, drugType.freq.ssri.clean[drugType.freq.ssri.clean$Data.ID %in% survival.data$LastDate.ID,]$FINNGENID) )
colnames(LastDate.Dose) <- c("dosePerDayUnits.Last", "FINNGENID")
LastDate.Dose$dosePerDayUnits.Last <- as.numeric(as.character(LastDate.Dose$dosePerDayUnits.Last))
LastDate.Dose$dosePerDayUnits.Last <- ifelse(LastDate.Dose$dosePerDayUnits.Last=="Inf", -99, LastDate.Dose$dosePerDayUnits.Last)

toClean <- merge(by.x="FINNGENID", by.y="FINNGENID", x=FirstDate.Dose, y=LastDate.Dose, all=T)
#inf looks at the first and only day they purchased?
#When NA, can I 0 it assuming that these patients did not take SSRI after?
### Check if these patients later needed to switch maybe?
### Create list of first and last to see what they switched to? For now, 0 them! 
#
table(toClean$dosePerDayUnits.First=="Inf")

#FALSE  
#48886  
#R6
#FALSE
#11494

table(toClean$dosePerDayUnits.Last=="Inf")

#FALSE
#48499
#R6
#FALSE
#42879
toClean$dosePerDayUnits.Last <- as.numeric(as.character(toClean$dosePerDayUnits.Last))
table(is.na(toClean$dosePerDayUnits.Last))
#FALSE  TRUE
#42879  2401
summary(toClean$dosePerDayUnits.Last)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#-99.0000   0.4444   0.9434  -2.6879   1.6807 300.0000
# FALSE
# 48510
#R6
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#0.2616  0.7317  1.0000  1.2442  1.5935  4.2553    2441

 toClean$dosePerDayUnits.First <- as.numeric(as.character(toClean$dosePerDayUnits.First))
 
 table(is.na(toClean$dosePerDayUnits.First))
 #FALSE  TRUE
 #10966 37544
 #R6
 #FALSE  TRUE
# 11593 33783
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# -99.00    0.02    0.12  -11.75    0.65  100.00   37544
#R6
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#  0.2885  0.7407  1.0000  1.2497  1.6000  4.3478    2401

toClean$dosePerDayUnits.First <- ifelse(is.na(toClean$dosePerDayUnits.First), 0, toClean$dosePerDayUnits.First)

toClean$dosePerDayUnits.First <- ifelse(toClean$dosePerDayUnits.First== -99, NA, toClean$dosePerDayUnits.First)
toClean$dosePerDayUnits.Last <- ifelse(toClean$dosePerDayUnits.Last== -99, NA, toClean$dosePerDayUnits.Last)

 
toClean$doseDiff <-  toClean$dosePerDayUnits.Last - toClean$dosePerDayUnits.First
summary(toClean$doseDiff)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#-99.720   0.074   0.732   1.173   1.327 249.938    3657
#R6
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#-3.5417  0.5102  0.9333  1.0041  1.3889  4.2553    2441
survival.data <- survival.data[order(survival.data$FINNGENID),]
toClean <- toClean[order(toClean$FINNGENID),]

survival.data<-merge(by.x="FINNGENID", by.y="FINNGENID", x=toClean, y=survival.data, all=T)

survival.data <- survival.data[!is.na(survival.data$doseDiff) & !is.na(survival.data$dosePerDayUnits.First) & !is.na(survival.data$dosePerDayUnits.Last),]

summary(survival.data$doseDiff)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#-99.71989   0.07389   0.73171   1.17338   1.32743 249.93773
#R6
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-3.5417  0.5102  0.9333  1.0041  1.3889  4.2553

table(survival.data$doseDiff<0)

#FALSE  TRUE
#43113  4092
#R6
#FALSE  TRUE
#40233  2646
table(survival.data$doseDiff>0)

#FALSE  TRUE
#10436 36769
#R6
#FALSE  TRUE
#5027 37852
table(survival.data$doseDiff==0)

#FALSE  TRUE
#40861  6344
#R6
#FALSE  TRUE
#40498  2381

quantile(survival.data$doseDiff, c(.05, .1, .2, .8, .9, .95))
#5%        10%        20%        80%        90%        95%
# -0.5714286  0.0000000  0.0000000  1.6666667  2.9159864  4.4545455
#R6
#        5%        10%        20%        80%        90%        95%
#-0.1526805  0.0000000  0.4424779  1.6129032  2.1428571  2.7777778

quantile(survival.data[survival.data$doseDiff<0,]$doseDiff, c(.05, .1, .2, .8, .9, .95))
#5%         10%         20%         80%         90%         95%
#-9.59702861 -5.65885592 -2.91531313 -0.18325484 -0.07544806 -0.03338790
#R6
#         5%         10%         20%         80%         90%         95%
#-2.13049342 -1.72341101 -1.23849211 -0.16265173 -0.07843846 -0.04003765

quantile(survival.data[survival.data$doseDiff>0,]$doseDiff, c(.05, .1, .2, .8, .9, .95))
#5%        10%        20%        80%        90%        95%
# 0.1044493 0.2251486 0.4372958 2.0000000 3.4336597 5.3571429
#R6
#       5%       10%       20%       80%       90%       95%
#0.3688525 0.4516129 0.5833333 1.7287675 2.2400000 2.8571429

#Clean up the outliers!!! -5 to 5

nrow(survival.data[survival.data$doseDiff > -5.5 & survival.data$doseDiff < 5.5 ,] )
#[1] 44987 from 45245
survival.data <- survival.data[survival.data$doseDiff > -5.5 & survival.data$doseDiff < 5.5 ,]

#bin creation: censored vs increase
survival.data$SSRIincrease_surv <- ifelse(survival.data$doseDiff>0, 1, 
                                          ifelse(survival.data$doseDiff==0, 0, NA))
#bin creation: censored vs decrease
survival.data$SSRIdecrease_surv <- ifelse(survival.data$doseDiff<0, 1, 
                                          ifelse(survival.data$doseDiff==0, 0, NA))
summary(survival.data$doseDiff)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-5.47619  0.06498  0.69444  0.82419  1.19048  5.49521
#R6
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#-3.5417  0.5102  0.9333  1.0041  1.3889  4.2553

summary(survival.data$difference.days)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#      0      71     421    1329    1901    8761
#R6
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0     119     525    1486    2168    9126

quantile(survival.data$difference.days, c(.05, .1, .2, .8, .9, .95))
#5%    10%    20%    80%    90%    95%
#   0    0   37 2522 4195 5467
#R6
#    5%    10%    20%    80%    90%    95%
#0.0   32.0   88.0 2836.0 4569.2 5889.3

pdf(paste(todaysDate, "_R6_diffDaysVdoseDiff.surv.pdf", sep="") )

plot(survival.data$difference.days, survival.data$doseDiff)
dev.off()


survival.data$obsNum_days <- survival.data$numObsv / survival.data$difference.days
#"Inf" were people who purchased the drug twice in a day
summary(survival.data$obsNum_days)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.0002502 0.0092369 0.0148620       Inf 0.0409005       Inf
#R6
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#0.000270 0.009778 0.014286      Inf 0.028571      Inf
survival.data<-survival.data[survival.data$obsNum_days!="Inf",]

quantile(survival.data$obsNum_days, c(.05, .1, .2, .8, .9, .95))
#5%         10%         20%         80%         90%         95%
#  0.002801120 0.004664904 0.007434944 0.028571429 0.048780488 0.068965517
#R6
#         5%         10%         20%         80%         90%         95%
#0.003894601 0.006025479 0.008660716 0.029411765 0.047619048 0.064516129

survival.data$difference.days_30 <- (survival.data$difference.days / 30) - survival.data$numObsv + 1
quantile(survival.data$difference.days_30, c(.05, .1, .2, .8, .9, .95))
#5%          10%          20%          80%          90%          95%
# -0.200000   0.100000   1.666667  74.353333 116.333333 149.160000
#survival.data <- survival.data[survival.data$difference.days_30 < 180 & survival.data$difference.days_30 > -1,]
#R6
#         5%         10%         20%         80%         90%         95%
#-0.1666667   0.1333333   1.5000000  73.9933333 118.5400000 153.5533333

pdf(paste(todaysDate, "_obsNum_days.surv.pdf", sep="") )
ggplot(survival.data[survival.data$obsNum_days<0.1,], aes(x=obsNum_days)) + 
  geom_histogram(color="black", fill="white") + 
  theme_classic() + 
  labs(title=paste(todaysDate, "_obsNum_days.surv"))
ggplot(survival.data[survival.data$obsNum_days>0.1,], aes(x=obsNum_days)) + 
  geom_histogram(color="black", fill="white") + 
  theme_classic() + 
  labs(title=paste(todaysDate, "_obsNum_days.surv"))

ggplot(survival.data, aes(x=difference.days_30)) + 
  geom_histogram(color="black", fill="white") + 
  theme_classic() + 
  labs(title=paste(todaysDate, "difference.days_30.surv"))

dev.off()

todaysDate <- format.Date(Sys.Date(), "20%y%m%d")

pdf(paste(todaysDate, "_distPurchaseDates.surv.pdf", sep="") )

ggplot(survival.data, aes(x=difference.days)) + 
  geom_histogram(color="black", fill="white") + 
  theme_classic() + 
  labs(title=paste(todaysDate, "distPurchaseDates.surv"))

ggplot(survival.data, aes(x=doseDiff)) + 
  geom_histogram(color="black", fill="white") + 
  theme_classic() + 
  labs(title=paste(todaysDate, "distPurchaseDoses.surv"))

ggplot(survival.data[survival.data$doseDiff<0,], aes(x=doseDiff)) + 
  geom_histogram(color="black", fill="white") + 
  theme_classic() + 
  labs(title=paste(todaysDate, "distPurchaseDoses.surv.decrease"))

ggplot(survival.data[survival.data$doseDiff==0,], aes(x=doseDiff)) + 
  geom_histogram(color="black", fill="white") + 
  theme_classic() + 
  labs(title=paste(todaysDate, "distPurchaseDoses.surv.increase"))

ggplot(survival.data[survival.data$doseDiff>0,], aes(x=doseDiff)) + 
  geom_histogram(color="black", fill="white") + 
  theme_classic() + 
  labs(title=paste(todaysDate, "distPurchaseDoses.surv.noDiff"))

dev.off()


#choosing and renaming columns for merging
survival.data.merging <- as.data.frame(cbind(as.character(survival.data$FINNGENID), as.numeric(survival.data$difference.days),
                                             as.numeric(survival.data$SSRIincrease_surv), as.numeric(survival.data$SSRIdecrease_surv),
                                             as.numeric(as.character(survival.data$doseDiff))) )
colnames(survival.data.merging) <- c("FINNGENID", "SSRI_time.days", "SSRI_increase.surv","SSRI_decrease.surv", "SSRI_doseDiff"  )
survival.data.merging <- survival.data.merging[order(survival.data.merging$FINNGENID),]
phenos.cov <- phenos.cov[order(phenos.cov$FINNGENID),]
#phenos.cov <-  phenos.cov[,c(1:3248,3256:3263)]
phenos.cov <- merge(by.x="FINNGENID", by.y="FINNGENID", x=phenos.cov, y=survival.data.merging, all.x=T)

save.image(file=paste(todaysDate,"RData", sep=".") )

phenos.cov$SSRI_time.days <- as.numeric(phenos.cov$SSRI_time.days)
phenos.cov$SSRI_time.months <- as.numeric(as.character(phenos.cov$SSRI_time.days)) / 30
phenos.cov$SSRI_time.months <- round(phenos.cov$SSRI_time.months, 1)
summary(phenos.cov$SSRI_time.months)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#         0.00   40.10   89.25   98.17  152.50  222.80  181910
#R6
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#0.00   35.70   85.60   97.93  156.10  232.10  221262

phenos.cov$SSRI_time.weeks <- as.numeric(as.character(phenos.cov$SSRI_time.days)) / 7
phenos.cov$SSRI_time.weeks <- round(phenos.cov$SSRI_time.weeks, 1)
summary(phenos.cov$SSRI_time.weeks)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#     0.1   171.9   382.5   420.7   653.5   955.0  181910
#R6
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#0.1   153.0   367.0   419.7   669.0   994.9  221262

phenos.cov$SSRI_time.years <- as.numeric(as.character(phenos.cov$SSRI_time.days)) / 365
phenos.cov$SSRI_time.years <- round(phenos.cov$SSRI_time.years, 1)
summary(phenos.cov$SSRI_time.years)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
# 0.00    3.30    7.30    8.07   12.50   18.30  181910
#R6
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#0.00    2.90    7.00    8.05   12.80   19.10  221262

phenos.cov$SSRI_doseDiff <- as.numeric(as.character(phenos.cov$SSRI_doseDiff))
#phenos.cov$SSRI_increase.surv <- ifelse(phenos.cov$SSRI_increase.surv==0, 1, 10)
phenos.cov$SSRI_increase.surv <- as.numeric(as.character(phenos.cov$SSRI_increase.surv))
#assuming a <0.5 would be little change, and close to 0
phenos.cov$SSRI_increase.surv <-
  ifelse(phenos.cov$SSRI_doseDiff < 0.5 & !is.na(phenos.cov$SSRI_increase.surv), 0, phenos.cov$SSRI_increase.surv)
table(phenos.cov$SSRI_increase.surv)
#0     1
# 4921 31594

#phenos.cov$SSRI_decrease.surv <- ifelse(phenos.cov$SSRI_decrease.surv==1, 0, 1)
phenos.cov$SSRI_decrease.surv <- as.numeric(as.character(phenos.cov$SSRI_decrease.surv))
#assuming a <0.5 would be little change, and close to 0
phenos.cov$SSRI_decrease.surv <- 
  ifelse(phenos.cov$SSRI_doseDiff > -0.5 & !is.na(phenos.cov$SSRI_decrease.surv) , 0, phenos.cov$SSRI_decrease.surv)
table(phenos.cov$SSRI_decrease.surv)
#0    1
#1255 1314

#loop through time series with ~800 increase, break the data down by ~10 parts -- max was 8183 to 40 parts, 204.575 / part
#changed the days to month, check this line below and the plotting line

  
for(i in seq(from=21.78, to= 217.80, by=21.78)){
  tryCatch({ 
    print(paste("filter",i) )
    temp <- phenos.cov[phenos.cov$SSRI_time.months < i,]
    
    cox.increase <- coxph(Surv(SSRI_time.months, SSRI_increase.surv) ~ SSRIswitchStatus, 
                          na.action=na.exclude, data = temp)
    print(paste("cox ph increase " ) )
    print(summary(cox.increase))
    fit.increase <- survfit(cox.increase)
    fit.increase.fit <- survfit(Surv(SSRI_time.months, SSRI_increase.surv) ~ SSRIswitchStatus,
                                na.action=na.exclude, data = temp)
    
    cox.decrease <- coxph(Surv(SSRI_time.months, SSRI_decrease.surv)~SSRIswitchStatus,
                          na.action=na.exclude, data = temp)
    print(paste("cox ph decrease " ) )
    print(summary(cox.decrease))
    #concordance  0.533
    fit.decrease <- survfit(cox.decrease)
    fit.decrease.fit <- survfit(Surv(SSRI_time.months, SSRI_decrease.surv) ~ SSRIswitchStatus,
                                na.action=na.exclude, data = temp)
    
    pdf(paste(todaysDate, "_lessThan", i, ".months_coxPH.dosage.pdf", sep=""))
    p1 <- autoplot(fit.increase) 
    p2 <- autoplot(fit.increase.fit) 
    
    p3 <- autoplot(fit.decrease)  
    p4 <- autoplot(fit.decrease.fit)
    grid.arrange(p1,p2,p3,p4, ncol=2, top=textGrob(paste(todaysDate, "_lessThan", i, ".months_coxPH.dosage", sep=""), 
                                                gp=gpar(fontsize=12, font = 2)))
    
    dev.off()
     }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


##concordance as “the probability of agreement for any two randomly chosen observations, 
#where in this case agreement means that the observation with the shorter survival time of the two 
#also has the larger risk score. The predictor (or risk score) will often be the result of a Cox model 
#or other regression” and notes that: “For continuous covariates concordance is equivalent to Kendall’s tau, 
#and for logistic regression is is equivalent to the area under the ROC curve.”
# fitted models have a concordance index between 0.55 and 0.7 which is due to the noise present in the dat
#+BATCH_AxiomGT1_b01_V2P2.calls+BATCH_AxiomGT1_b02_V2P2.calls+BATCH_AxiomGT1_b03_V2P2.calls+BATCH_AxiomGT1_b04_V2P2.calls+BATCH_AxiomGT1_b05_V2P2.calls+BATCH_AxiomGT1_b06_V2P2.calls+BATCH_AxiomGT1_b07_V2P2.calls+BATCH_AxiomGT1_b08_V2P2.calls+BATCH_AxiomGT1_b09_V2P2.calls+BATCH_AxiomGT1_b10_V2P2.calls+BATCH_AxiomGT1_b11_V2P2.calls+BATCH_AxiomGT1_b12_V2P2.calls+BATCH_AxiomGT1_b13_V2P2.calls+BATCH_AxiomGT1_b14_V2P2.calls+BATCH_AxiomGT1_b15_V2P2.calls+BATCH_AxiomGT1_b17_V2P2.calls+BATCH_AxiomGT1_b18_V2P2.calls+BATCH_AxiomGT1_b19_V2P2.calls+BATCH_AxiomGT1_b20_V2P2.calls+BATCH_AxiomGT1_b21_V2P2.calls+BATCH_AxiomGT1_b22_V2P2.calls+BATCH_AxiomGT1_b23_V2P2.calls+BATCH_AxiomGT1_b24_V2P2.calls+BATCH_AxiomGT1_b25_V2P2.calls+BATCH_AxiomGT1_b26_V2P2.calls+BATCH_AxiomGT1_b27_V2P2.calls+BATCH_AxiomGT1_b28_V2P2.calls+BATCH_AxiomGT1_b29_V2P2.calls+BATCH_AxiomGT1_b30_V2P2.calls+BATCH_AxiomGT1_b31_V2.calls+BATCH_AxiomGT1_b3234_V2.calls+BATCH_AxiomGT1_b33_V2.calls+BATCH_AxiomGT1_b35_V2.calls+BATCH_AxiomGT1_b36_V2.calls+BATCH_AxiomGT1_b37_V2.calls+BATCH_AxiomGT1_b38_V2.calls+BATCH_AxiomGT1_b39_V2.calls+BATCH_AxiomGT1_b40_V2.calls+BATCH_DS10_FINRISK_Palotie_norm+BATCH_DS11_FINRISK_PredictCVD_COROGENE_Tarto_norm+BATCH_DS12_FINRISK_Summit_norm+BATCH_DS13_FINRISK_Bf_norm+BATCH_DS14_GENERISK_norm+BATCH_DS15_H2000_Broad_norm+BATCH_DS16_H2000_Fimm_norm+BATCH_DS17_H2000_Genmets_norm+BATCH_DS18_MIGRAINE_1_norm+BATCH_DS19_MIGRAINE_2_norm+BATCH_DS1_BOTNIA_Dgi_norm+BATCH_DS20_SUPER_1_norm+BATCH_DS21_SUPER_2_norm+BATCH_DS22_TWINS_1_norm+BATCH_DS23_TWINS_2_norm+BATCH_DS24_SUPER_3_norm+BATCH_DS2_BOTNIA_T2dgo_norm+BATCH_DS3_COROGENE_Sanger_norm+BATCH_DS4_FINRISK_Corogene_norm+BATCH_DS5_FINRISK_Engage_norm+BATCH_DS6_FINRISK_FR02_Broad_norm+BATCH_DS7_FINRISK_FR12_norm+BATCH_DS8_FINRISK_Finpcga_norm+BATCH_DS9_FINRISK_Mrpred_norm, 
# ~ SEX_IMPUTED+AGE_AT_DEATH_OR_NOW+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10


todaysDate <- format.Date(Sys.Date(), "20%y%m%d")
write.table(phenos.cov, paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="\t", quote = F, 
            col.names = T, row.names = F)
command<- paste("gzip -9 -v -f ", paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="")
system(command)
phenos.list <- c(colnames(phenos.cov[grep("SSRI_time", colnames(phenos.cov))]), "SSRI_increase.surv", "SSRI_decrease.surv"  )
phenos.list
write.table(phenos.list, paste(todaysDate, "phenoList_time.txt", sep="_"), sep="\t", quote=F, col.names = F, row.names = F)
save.image(file=paste(todaysDate,"RData", sep=".") )

#survival analysis didn't like NAs

phenos.cov.SSRI_surv <- phenos.cov[!is.na(phenos.cov$SSRI_time.months) & !is.na(phenos.cov$SSRI_increase.surv),]
#ID, AGE, BATCH, SEX, PCs, SSRI_increase.surv, SSRI_time.months
#phenos.cov.SSRI_increase.surv <- phenos.cov.SSRI_surv[,c(1:2, 8:70, 90:100, 3258, 3261)]
#R6
phenos.cov.SSRI_increase.surv <- phenos.cov.SSRI_surv[,c(1:2, 8:82, 102:113, 3016, 3013)]
### need to put in SEX!!! 
#phenos.cov.SSRI_increase.surv <- phenos.cov.SSRI_surv[,-c(79)]
write.table(phenos.cov.SSRI_increase.surv, paste(todaysDate, ".SSRI_increase.surv_R6_COV_PHENO_KV.txt", sep=""), sep="\t", quote = F, 
            col.names = T, row.names = F)

phenos.cov.SSRI_surv <- phenos.cov[!is.na(phenos.cov$SSRI_time.months) & !is.na(phenos.cov$SSRI_decrease.surv),]
#phenos.cov.SSRI_decrease.surv <- phenos.cov.SSRI_surv[,c(1:2, 8:70, 90:100, 3259, 3261)]
phenos.cov.SSRI_decrease.surv <- phenos.cov.SSRI_surv[,c(1:2, 8:82, 102:113, 3016, 3014)]

#phenos.cov.SSRI_decrease.surv <- phenos.cov.SSRI_surv[,-c(78)]
write.table(phenos.cov.SSRI_decrease.surv, paste(todaysDate, ".SSRI_decrease.surv_R6_COV_PHENO_KV.txt", sep=""), sep="\t", quote = F, 
            col.names = T, row.names = F)

## Are there deviations from the median
drugType.freq.ssri.clean <- drugType.freq.ssri.clean[drugType.freq.ssri.clean$FINNGENID %in% survival.data.merging$FINNGENID, ]

SSRI_increase.surv.IDs <- phenos.cov.SSRI_increase.surv[,c(1,78,77)]
SSRI_decrease.surv.IDs <- phenos.cov.SSRI_decrease.surv[,c(1,78,77)]

colnames(SSRI_decrease.surv.IDs) <- c("FINNGENID", "SSRI_time.days.decrease", "SSRI_decrease.surv")
colnames(SSRI_increase.surv.IDs) <- c("FINNGENID", "SSRI_time.days.increase", "SSRI_increase.surv")

drugType.freq.ssri.clean <- drugType.freq.ssri.clean[order(drugType.freq.ssri.clean$FINNGENID),]
SSRI_increase.surv.IDs <- SSRI_increase.surv.IDs[order(SSRI_increase.surv.IDs$FINNGENID),]
drugType.freq.ssri.clean <- merge(by.x="FINNGENID", by.y="FINNGENID", x=drugType.freq.ssri.clean, y=SSRI_increase.surv.IDs, all.x=T)

drugType.freq.ssri.clean <- drugType.freq.ssri.clean[order(drugType.freq.ssri.clean$FINNGENID),]
SSRI_decrease.surv.IDs <- SSRI_decrease.surv.IDs[order(SSRI_decrease.surv.IDs$FINNGENID),]
drugType.freq.ssri.clean <- merge(by.x="FINNGENID", by.y="FINNGENID", x=drugType.freq.ssri.clean, y=SSRI_decrease.surv.IDs, all.x=T)

table(drugType.freq.ssri.clean$SSRI_increase.surv)

#0      1
#51190 364389

table(drugType.freq.ssri.clean$SSRI_decrease.surv)

#    0     1
# 6790 11389

ssriDrugUnits.increase <- as.data.frame(drugType.freq.ssri.clean[drugType.freq.ssri.clean$SSRI_increase.surv==1,] %>% 
                                          group_by(FINNGENID) %>%
                                          summarize(SSRImeanDose=mean(dosePerDayUnits), SSRIsdDose=sd(dosePerDayUnits), 
                                                    SSRImedianDose=median(dosePerDayUnits), 
                                                    SSRInum=n(), SSRImadDDD = mad(dosePerDayUnits),
                                                    Q1=quantile(dosePerDayUnits, probs=0.25, na.rm=T), 
                                                    Q3=quantile(dosePerDayUnits, probs=0.75, na.rm=T)
                                          ) ) 


ssriDrugUnits.decrease <- as.data.frame(drugType.freq.ssri.clean[drugType.freq.ssri.clean$SSRI_decrease.surv==1,] %>% 
                                group_by(FINNGENID) %>%
                                summarize(SSRImeanDose=mean(dosePerDayUnits), SSRIsdDose=sd(dosePerDayUnits), 
                                          SSRImedianDose=median(dosePerDayUnits), 
                                          SSRInum=n(), SSRImadDDD = mad(dosePerDayUnits),
                                          Q1=quantile(dosePerDayUnits, probs=0.25, na.rm=T), 
                                          Q3=quantile(dosePerDayUnits, probs=0.75, na.rm=T)
                                ) )

pdf(paste(todaysDate, "_SSRI_doseMAD.pdf", sep=""))
ggplot(ssriDrugUnits.increase, aes(x=SSRImadDDD)) + 
  geom_histogram(color="black", fill="white") + 
  theme_classic() + 
  labs(title=paste(todaysDate, "_SSRI_doseMAD.increase"))

ggplot(ssriDrugUnits.decrease, aes(x=SSRImadDDD)) + 
  geom_histogram(color="black", fill="white") + 
  theme_classic() + 
  labs(title=paste(todaysDate, "_SSRI_doseMAD.decrease"))
dev.off()

decreaseIDs <- 
  nrow(ssriDrugUnits.decrease[ssriDrugUnits.decrease$SSRImadDDD>2.0 & 
                                        !is.na(ssriDrugUnits.decrease$SSRImadDDD) ,]$FINNGENID )
increaseIDs <- ssriDrugUnits.increase[ssriDrugUnits.increase$SSRImadDDD>2.0 &
                                        !is.na(ssriDrugUnits.increase$SSRImadDDD),]$FINNGENID


phenos.cov.SSRI_increase.surv <- phenos.cov.SSRI_increase.surv[!(phenos.cov.SSRI_increase.surv$FINNGENID %in% increaseIDs),]
write.table(phenos.cov.SSRI_increase.surv, paste(todaysDate, ".SSRI_increase.surv_R6_COV_PHENO_KV.txt", sep=""), sep="\t", quote = F, 
            col.names = T, row.names = F)

phenos.cov.SSRI_decrease.surv <- 
  phenos.cov.SSRI_decrease.surv[!(phenos.cov.SSRI_decrease.surv$FINNGENID %in% decreaseIDs),]
write.table(phenos.cov.SSRI_decrease.surv, paste(todaysDate, ".SSRI_decrease.surv_R6_COV_PHENO_KV.txt", sep=""), sep="\t", quote = F, 
            col.names = T, row.names = F)

pdf(paste(todaysDate, "SSRIdosagePRS.pdf", sep="") )

ggplot(phenos.cov[!is.na(phenos.cov$SSRI_increase.surv),], aes(x = MDD_PRS, fill=as.factor(SSRI_increase.surv) ) )  +
  geom_density(alpha = .3) +
  theme_classic() + 
  labs(title=paste(todaysDate, "PRS for SSRI increase : ", plot.name,  " from MDD summstats", sep="") )

ggplot(phenos.cov[!is.na(phenos.cov$SSRI_decrease.surv),], aes(x = MDD_PRS, fill=as.factor(SSRI_decrease.surv) ) )  +
  geom_density(alpha = .3) +
  theme_classic() + 
  labs(title=paste(todaysDate, "PRS for SSRI decrease : ", plot.name,  " from MDD summstats", sep="") )

dev.off()


#ideas for survival analysis?
### increase / decrease dosage
### switch ? 

  
  
  
 #### Ambiguity in R5 vs R6 ####
setwd("~/workDir/r5")
phenos.cov.r5 <- fread("zcat 20200713_R5_COV_PHENO_KV.txt.gz", data.table=F)
phenos.cov.r5.drugs <- fread("zcat R5_COV_PHENO_V1_DRUGS_V1_ALL.txt.gz", data.table=F)

#Difference here is because of the AGE >35
  
  



################################
######### F5 Phenotypes ########
################################

#Because of this
summary(phenos.cov[phenos.cov$SSRIfreq>0,]$F5_DEPRESSIO)
#0      1   NA's
#34248  19036 165508
summary(phenos.cov$F5_DEPRESSIO)
#0      1   NA's
#192220  23424   3148

F5phenos <- grep("F5_*",colnames(phenos.cov), value=T)

for(i in 1:length(F5phenos)){
  print(F5phenos[i])
  temp <- paste(F5phenos[i], "_SSRI",sep="")
  print(temp)
  phenos.cov[,temp] <- NULL
  phenos.cov[,temp] <- 0
  phenos.cov[,temp]<- ifelse(phenos.cov$SSRIfreq>0, as.numeric(as.character(phenos.cov[,F5phenos[i]])), NA)
  print(table(phenos.cov[,temp]))
}


todaysDate <- format.Date(Sys.Date(), "20%y%m%d")
write.table(phenos.cov, paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="\t", quote = F, 
            col.names = T, row.names = F)
command<- paste("gzip -9 -v -f ", paste(todaysDate, "_R6_COV_PHENO_KV.txt", sep=""), sep="")
system(command)
phenos.list <- colnames(phenos.cov[grep("F5_", colnames(phenos.cov))] )

write.table(phenos.list, paste(todaysDate, "phenoList.txt", sep="_"), sep="\t", quote=F, col.names = F, row.names = F)
save.image(file=paste(todaysDate,"RData", sep=".") )


##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### ##### ##### Memory Use Function ##### ##### ##### ##### ###
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

showMemoryUse <- function(sort="size", decreasing=FALSE, limit) {
  
  objectList <- ls(parent.frame())
  
  oneKB <- 1024
  oneMB <- 1048576
  oneGB <- 1073741824
  
  memoryUse <- sapply(objectList, function(x) as.numeric(object.size(eval(parse(text=x)))))
  
  memListing <- sapply(memoryUse, function(size) {
    if (size >= oneGB) return(paste(round(size/oneGB,2), "GB"))
    else if (size >= oneMB) return(paste(round(size/oneMB,2), "MB"))
    else if (size >= oneKB) return(paste(round(size/oneKB,2), "kB"))
    else return(paste(size, "bytes"))
  })
  
  memListing <- data.frame(objectName=names(memListing),memorySize=memListing,row.names=NULL)
  
  if (sort=="alphabetical") memListing <- memListing[order(memListing$objectName,decreasing=decreasing),] 
  else memListing <- memListing[order(memoryUse,decreasing=decreasing),] #will run if sort not specified or "size"
  
  if(!missing(limit)) memListing <- memListing[1:limit,]
  
  print(memListing, row.names=FALSE)
  return(invisible(memListing))
}

# to use function
showMemoryUse(decreasing=TRUE, limit=5)

#### Collapsing Batch Effects
phenos.cov$BATCH_Axiom <- phenos.cov$BATCH_AxiomGT1_b01_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b02_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b03_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b04_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b05_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b06_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b07_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b08_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b09_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b10_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b11_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b12_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b13_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b14_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b15_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b17_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b18_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b19_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b20_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b21_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b22_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b23_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b24_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b25_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b26_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b27_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b28_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b29_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b30_V2P2.calls +
  + phenos.cov$BATCH_AxiomGT1_b31_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b3234_V2.calls + phenos.cov$BATCH_AxiomGT1_b33_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b35_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b36_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b37_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b38_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b39_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b40_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b41_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b42_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b43_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b44_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b45_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b46_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b47_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b48_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b49_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b50_V2.calls +
  + phenos.cov$BATCH_AxiomGT1_b51_V2.calls

################################################################################################
################################################################################################
############################# ARCHIVE ##########################################################
################################################################################################
################################################################################################
#############################
#### without drug dosage ####
#############################
drug.purchase.freq <- as.data.frame(table(data.purchases.vnr.drug$FINNGENID) )
drug.purchase.freq$label <- "drug"
summary(drug.purchase.freq$Freq)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.00    2.00    9.00   22.19   29.00  585.00

unique.drug <- data.purchases.vnr.drug[!duplicated(data.purchases.vnr.drug$FINNGENID),]
data.purchases.vnr.drug <- data.purchases.vnr.drug[order(data.purchases.vnr.drug$APPROX_EVENT_DAY),]
data.purchases.drug.bursts <- as.data.frame(data.purchases.vnr.drug %>% 
                                              group_by(FINNGENID) %>%
                                              mutate(daysPassed.diff=c(NA,diff(APPROX_EVENT_DAY))) )

data.purchases.drug.bursts <- data.purchases.drug.bursts[!is.na(data.purchases.drug.bursts$daysPassed.diff),]
summary(data.purchases.drug.bursts$daysPassed.diff)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.0    40.0    89.0   150.6   111.0  8516.0

quantile(data.purchases.drug.bursts$daysPassed.diff, c(.05, .1, .2, .8, .9, .95))
#   5% 10% 20% 80% 90% 95%
#  14  22  33 120 184 372

median.plotting <- 
  as.data.frame(data.purchases.drug.bursts %>%
                  group_by(FINNGENID) %>%
                  summarise(M=mean(daysPassed.diff), Med=median(daysPassed.diff), Q1=quantile (daysPassed.diff, probs=0.25), 
                            Q2=quantile(daysPassed.diff, probs=0.50), Q3=quantile(daysPassed.diff, probs=0.75), 
                            Q4=quantile(daysPassed.diff, probs=1.00)) ) 


pdf("finngenFreq.pKillerNonSteroid.pdf")

ggplot(drug.purchase.freq, aes(x=reorder(Var1, -Freq), y=Freq)) + geom_col() +
  labs(title="painkillers purchase")  +
  theme_classic()

ggplot(drug.purchase.freq, aes(x=Freq)) + geom_bar() +
  labs(title="painkillers purchase Freq dist") +
  theme_classic()

ggplot(drug.purchase.freq, aes(x=label, y=Freq)) +
  geom_boxplot() +
  stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  geom_jitter(shape=16, position=position_jitter(0.02)) +
  labs(title="painkillers purchase Freq dist") +
  theme_classic()

ggplot(median.plotting, aes(x=reorder(FINNGENID, Med), y=Med)) +
  geom_point(size=.2, shape=21, fill="white") +
  geom_errorbar(aes(ymin=Q1, ymax=Q3), width=.01, position = position_dodge2(width = 0.05, padding = 0.05), linetype="4C88C488") +
  theme_classic()

ggplot(data.purchases.drug.bursts[!is.na(data.purchases.drug.bursts$daysPassed.diff),], 
       aes(x=reorder(FINNGENID, daysPassed.diff, FUN=median), y=daysPassed.diff)) +
  geom_boxplot(outlier.size = 0) +
  # stat_summary(fun.y=median, geom="point", shape=23, size=4) +
  labs(title="painkillers purchase Freq dist") +
  theme_classic()

dev.off()



####NORMALIZING?

drugType.freq.ssri$dosePerDay_IRN <- 0
drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Citalopram" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay_IRN <- 
  rankNorm(drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Citalopram" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay)

drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Escitalopram" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay_IRN <- 
  rankNorm(drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Escitalopram" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay)

drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Fluoxetine" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay_IRN <- 
  rankNorm(drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Fluoxetine" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay)

drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Fluvoxamine" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay_IRN <- 
  rankNorm(drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Fluvoxamine" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay)

drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Paroxetine" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay_IRN <- 
  rankNorm(drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Paroxetine" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay)

drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Sertraline" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay_IRN <- 
  rankNorm(drugType.freq.ssri[drugType.freq.ssri$generic.Name=="Sertraline" & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay)





#1861 CIPRALEX 011314 6028 FGV5DKFMYX  PURCH     25.73       2010-06-06 N06AB10 1861  <NA>     2   <NA>     <NA> 42308499 N06AB10   20 MG          20   100
#CIPRALEX 011314 FGV5DKFMYX

data.purchases.vnr.ssri[data.purchases.vnr.ssri$FINNGENID=="FGV5DKFMYX" & data.purchases.vnr.ssri$valmiste=="CIPRALEX", ]
summary(drugType.freq.ssri[drugType.freq.ssri$dosePerDay>200 & drugType.freq.ssri$daysPassed.diff!=0,]$daysPassed.diff)
drugType.freq.ssri <- drugType.freq.ssri[order(drugType.freq.ssri$APPROX_EVENT_DAY),]

# assume the DDDmg when daysPassed is <3 days,, pkoko_num > 20
drugType.freq.ssri[drugType.freq.ssri$dosePerDay>200 & drugType.freq.ssri$daysPassed.diff< 3 & drugType.freq.ssri$pkoko_num < 20 & drugType.freq.ssri$daysPassed.diff!=0,]
summary(drugType.freq.ssri[drugType.freq.ssri$dosePerDay < 200 & drugType.freq.ssri$daysPassed.diff > 3 & drugType.freq.ssri$pkoko_num > 20 & drugType.freq.ssri$daysPassed.diff!=0,]$dosePerDay)

