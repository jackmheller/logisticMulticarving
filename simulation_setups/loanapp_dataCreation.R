require(wooldridge)
require(dplyr)
data('loanapp')


# Data previously used here: https://www.researchgate.net/profile/William-Hunter-17/publication/5151431_The_Cultural_Gap_Hypothesis_and_Mortgage_Lending_Decisions/links/00b4951fa9f11a10f7000000/The-Cultural-Gap-Hypothesis-and-Mortgage-Lending-Decisions.pdf
# 
# They say that the original data is dirty (missing values, errors).
# The original data has 3,300+722 observations.
# They explain criterions for cleaning the data in section 2.
# 
# Also used here: https://www.jstor.org/stable/pdf/2118254.pdf?casa_token=Bx8AAZ9Rkv0AAAAA:OQquj6at0UQXJiXiOabnTtSM1sW3ziLMNqYzSEIXNCh-IvhVcxk5u1BzL4spW9pkDtbLzXs7GSg31JKjjDKdmEUtoQu19rcFfpSgOAmqour85hxLg1A
# 
# Information about the data is found in section 3.



dim(loanapp)

#creating one variable from reject and approve - label
loanapp$label<- ifelse(loanapp$approve==1, 1, 0)


unique(loanapp$msa)
#there is only one value for msa! We can drop that variable 

#removing action
loanapp<-loanapp[ , -which(names(loanapp) %in% c("action","reject", "approve", "msa"))]
# thick is removed because it is equal to rep
#I am assuming that gdlin==666 means missing. These two observations will be dropped
loanapp$gdlin<-ifelse(loanapp$gdlin==666, NA, loanapp$gdlin)

sum(is.na(loanapp)) #230
#for now just removing missing data 
loanapp<-na.omit(loanapp) 

# sum(loanapp[rowSums(loanapp[,c("black", "white", "hispan")]) == 0,])
#meaning that race is one hot encoded

#Seeing how many unique values each variable has
apply(loanapp, 2, function(x) length(unique(x)))
# Variables we assume to be categorical/ binary: occ, suffolk, typur, married,
# self, gdlin, mortg, cons, pubrec, fixadj, prop, inss, inson, gift, cosign, 
# unver, min30, bd, mi, old, vr, sch, mortno, mortperf, mortlat1, mortlat2, chist,
# multi
# for mortg and cons we are assuming that they are still encoded based on 
# "Mortgage Lending in Boston: Interpreting HMDA Data"

#----- Prep. Cat and binary variables --------#

loanapp$race<-0
loanapp$race<-ifelse(loanapp$black==1, "1",loanapp$race)
loanapp$race<-ifelse(loanapp$hispan==1, "2",loanapp$race)
loanapp$race<-ifelse(loanapp$white==1, "3",loanapp$race)


cols <- c( "occ", "suffolk", "typur", "married","self", "gdlin", "mortg", "cons",
           "pubrec", "fixadj", "prop", "inss", "inson", "gift", "cosign", "unver",
           "min30", "bd", "mi", "old", "vr", "sch", "mortno", "mortperf", "mortlat1", 
           "mortlat2", "chist","multi", "race")
loanapp[cols] <- lapply(loanapp[cols], factor)  ## as.factor() could also be used

#some of them I removed because I made factors out of them, some I removed cause of VIF
drop.cols<-c("black", "hispan", "white") 

loanapp<- loanapp%>% select(-one_of(drop.cols))

#--------- Dummy Encoding ---------#

colsNoNorm<-c("label")

loanapp.norm<-loanapp[ , -which(names(loanapp) %in% colsNoNorm)]

loanapp.norm <- model.matrix( ~ .-1, loanapp.norm)

#--------- Normalizing ---------#

loanapp.norm<-scale(as.data.frame(loanapp.norm))
loanapp.norm<-as.data.frame(loanapp.norm)

#--------- Interactions & Polynomials ---------#

loanapp.inter <- model.matrix( ~.^2, data=as.data.frame(loanapp.norm))[,-1]

# Add target back
loanapp.inter<-cbind(loanapp[,"label"], loanapp.inter)

save(loanapp.inter, file = "loanappInter.Rdata")

