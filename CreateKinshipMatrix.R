#!/usr/bin/env Rscript
#
# Script:   CreateKinshipMatrix.R
# Purpose:  Convert table of family identies twin status into a "kinship" matrix
# Author:   Thomas Nichols
# Version: $Id$
#
Usage='
 Usage:   CreateKinshipMatrix.R In.csv Out.csv

 In.csv must be a plain CSV file, with the first row containing the
 variable names and then one row per subject (no blank rows).
 The following variable names must be in the file: "Subject", 
 "Mother ID", "Father ID", and "Zygosity".  

 Out.csv will be a nSubject by nSubject matrix, with row and column
 names set by the variable "Subject".  The possible values indicate
 the possible type of relationships:

   -1     Self
    0     Unrelated
    1     MZ twins
    2     DZ twins
    3     Siblings

 NOTE:  Assumes there is no more than one twin pair per family.
'

args = commandArgs(TRUE)

if (length(args)!=2)
  stop(Usage)

InFn = args[[1]]
OutFn = args[[2]]

dat=read.csv(InFn)
nSubj=dim(dat)[1]

# Make unique identifier for each family
dat$FamID=factor(paste(dat$Mother.ID,dat$Father.ID))


MZ=DZ=Sib=matrix(0,ncol=nSubj,nrow=nSubj)

for (fam in levels(dat$FamID)) {

  # Any family relationship
  Ifam = (dat$FamID==fam)
  Sib = Sib + Ifam%*%t(Ifam)

  # MZ relationship
  Ifam = (dat$FamID==fam) & (dat$Zygosity == "MZ")
  MZ = MZ + Ifam%*%t(Ifam)

  # DZ relationship
  Ifam = (dat$FamID==fam) & (dat$Zygosity == "NotMZ")
  DZ = DZ + Ifam%*%t(Ifam)

}
  
# Do some nonintuitive arithmatic to get final matrix
Adj=Sib*3-MZ*2-DZ
Adj=Adj-diag(diag(Adj))-1*diag(nSubj)

# Add column names
colnames(Adj)<-as.character(dat$"Subject")
rownames(Adj)<-as.character(dat$"Subject")

write.csv(Adj,OutFn)
