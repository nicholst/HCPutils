#!/usr/bin/env Rscript
#
# Script:   CreateKinshipMatrix.R
# Purpose:  Convert table of family identies twin status into a "kinship" matrix
# Author:   Thomas Nichols
# Version:  http://github.com/nicholst/HCPutils/tree/$Format:%h$
#           $Format:%ci$
#
Usage='
 Usage:   CreateKinshipMatrix.R In.csv Out.csv

 In.csv must be a plain CSV file, with the first row containing the
 variable names and then one row per subject (no blank rows).
 The following variable names must be in the file: "Subject",
 "Mother_ID", "Father_ID", "ZygositySR" and, "ZygosityGT".  
 Zygosity is draw from ZygosityGT preferentially, and then as backup
 from ZygositySR.

 Out.csv will be a nSubject by nSubject matrix, with row and column
 names set by the variable "Subject".  The possible values indicate
 the possible type of relationships:

   -1     Self
    0     Unrelated
    1     MZ twins
    2     DZ twins
    3     Full siblings (FS)
    4     Half siblings (HS)

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
dat$FamID=factor(paste(dat$Mother_ID,dat$Father_ID))

MZ=DZ=FS=HS=matrix(0,ncol=nSubj,nrow=nSubj)

Zyg=as.character(dat$ZygosityGT)

# Fill in missing ZygosityGT values from ZygositySR
BadZ=!(Zyg=="MZ" | Zyg=="DZ")
Zyg[BadZ]=as.character(dat$ZygositySR[BadZ])

old=0
# First loop through mothers
AllFath=c()
Blend=c()
for (moth in levels(factor(dat$Mother_ID))) {
  Imoth = (dat$Mother_ID==moth)
  Faths=c()
  Fams=c()
  for (fath in unique(as.character(dat$Father_ID[Imoth]))) {
    fam=paste(moth,fath)
    Fams=c(Fams,fam)
    Faths=c(Faths,fath)
    AllFath=c(AllFath,fath)

    Ifam=dat$FamID==fam

    # Any family relationship
    Ifam = (dat$FamID==fam)
    FS = FS + Ifam%*%t(Ifam)

    # MZ relationship
    Ifam = (dat$FamID==fam) & (Zyg == "MZ")
    MZ = MZ + Ifam%*%t(Ifam)

    # DZ relationship
    Ifam = (dat$FamID==fam) & ( (Zyg == "DZ")|(Zyg == "NotMZ") )
    DZ = DZ + Ifam%*%t(Ifam)
  }
  if (length(Fams)>1) {
    for (i1 in 1:length(Fams)) {
      fam1=Fams[i1]
      Blend=c(Blend,which(dat$FamID==fam1))
      for (i2 in (i1+1):length(Fams)) {
        fam2 = Fams[i2]
        Ifam1 = dat$FamID==fam1
        Ifam2 = dat$FamID==fam2
        HS = HS + Ifam1%*%t(Ifam2) + Ifam2%*%t(Ifam1)
      }
    }
  }
  if (is.na(HS[90,90]) || HS[90,90]>old) {
    browser()
    old=HS[90,90]
  }

}
# Now catch fathers that appear in 2 or more families
DupFath=AllFath[duplicated(AllFath)]
if (length(DupFath)>0) {
  for (fath in DupFath) {
    Ifath = (dat$Father_ID==fath)
    Fams=unique(as.character(dat$FamID[Ifath]))
    for (i1 in 1:length(Fams)) {
      fam1=Fams[i1]
      Blend=c(Blend,which(dat$FamID==fam1))
      for (i2 in (i1+1):length(Fams)) {
        fam2=Fams[i2]
        Ifam1 = dat$FamID==fam1
        Ifam2 = dat$FamID==fam2
        HS = HS + Ifam1%*%t(Ifam2) + Ifam2%*%t(Ifam1)
      }
    }
  }
}

if length(Blend)>0 {
  cat("WARNING: Blended families:\n")
  print(cbind(dat[Blend,c("Subject","Mother_ID","Father_ID")],Zygosity=Zyg[Blend]))
  cat("\n")
}
browser()
# Do some nonintuitive arithmatic to get final matrix
Adj=FS*3-MZ*2-DZ
Adj=Adj-diag(diag(Adj))-1*diag(nSubj)
Adj=Adj+4*HS
# Add column names
colnames(Adj)<-as.character(dat$"Subject")
rownames(Adj)<-as.character(dat$"Subject")

write.csv(Adj,OutFn)
