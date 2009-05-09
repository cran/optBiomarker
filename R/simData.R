
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License 
## as published by the Free Software Foundation; either version 2 
## of the License, or (at your option) any later version.
##  
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General 
## Public License for more details.
##   
## You may also obtain a copy of the GNU General Public License from 
## the Free Software Foundation by visiting their Web site or by writing to
##   
##
##  Free Software Foundation, Inc.
##  59 Temple Place - Suite 330
##  Boston, MA 02111-1307
##  USA
##
################################################################################

## Function for generating nPt by nBiom data matrix for classification analysis 
## ----------------------------------------------------------------------------

## nTrain = Number of subjects in the training set (nGr1+nGr2)
## nGr1 = Number of subjects in Group 1
## nGr2 = Number of subjects in Group 2 
## nBiom = Number of biomarkers (genes, proteins)
## nRep= Number of technical replications

## sdW= sqrt(experimental or technical variation)
## sdB= sqrt(biological variation)
## rho= common Pearson correlation coefficient between biomarkers 
## foldMin= Minimum value of fold changes
## sigma= Standard deviation of the normal distribution (before truncation)
##        where fold changes are generated from  
## baseExpr = A vector of length nBiom to be used as \mu_g (log2 scale, <16) 

## Main function for simulating data
## ---------------------------------


simData<-function(nTrain=100, nGr1=floor(nTrain/2), nBiom=50,nRep=3,
                  sdW=1.0, sdB=1.0,rho=0,sigma=1,diffExpr=TRUE, foldMin=2,orderBiom=TRUE,baseExpr=NULL)
{
    checkInt<-c(nTrain,nGr1,nBiom, nRep)
    if(!identical(checkInt, floor(checkInt))) stop("non-integer(s) given where argument(s) should be integer valued") 

    if (!is.null(baseExpr) && length(baseExpr)!=nBiom)
    {stop("length of 'baseExpr' should be equal to nBiom") }

    if(nTrain<2) stop("training set size can not be less than 2")

    if( any(c(nGr1, nBiom, nRep)<1)) stop("'nGr1', 'nBiom' and 'nRep' can not be less than 1") 
    nGr2<-nTrain-nGr1

    if (nGr2<1) stop("'nTrain' must be greater than 'nGr1'")
    if (rho<0 || rho>0.95)stop("allowed values of 'rho' are between 0 and 0.95 inclusive")


    ## Simulate  nTrain by nBiom by nRep array of experimental errors
    ## --------------------------------------------------------------
    if (nRep>1) e<-array(rnorm(nTrain*nBiom*nRep,mean=0,sd=sdW),dim=c(nTrain,nBiom,nRep)) else
                e<-array(rnorm(nTrain*nBiom,mean=0,sd=sdW),dim=c(nTrain,nBiom))
    
    ## Average e over replicates
    if (nRep>1) eAvg<-apply(e,c(1,2),mean) else eAvg<-e
    
    ## Simulate  nTrain by nBiom matrix  of biological errors
    ## ------------------------------------------------------
    b<-matrix(rnorm(nTrain*nBiom, mean=0,sd=sdB),nTrain)
    
    ## Add biological and experimental errors
    ## ---------------------------------------
    randErr<-b+eAvg

    ## Introduce correlation structure
    ## -------------------------------

    if (nBiom>1 && rho!=0) {

      ## Define the correlation matrix

      Rho<-diag(rep(1,nBiom))
      Rho[upper.tri(Rho)|lower.tri(Rho)]<-rho
      
      ## Existing standard deviations in the data
      
      sdMat<-diag(sqrt(diag(var(randErr))))

      ## Scale the data to have unit variances
      randErrScaled<-scale(randErr,center=FALSE,scale=TRUE)
   
      ## Covariance matrix with desired Rho
      
      covMat<-sdMat%*%Rho%*%sdMat
      cholRoot<-chol(covMat) ## property: t(cholRoot)%*%cholRoot==covMat

      ## Transformed data to have desired covariance structure

      randErrTrans<-randErrScaled%*%cholRoot

      ## Replace original data by the transformed data

      randErr<- randErrTrans
    }

    ## Give column names
    ## -----------------
    colnames(randErr)<-paste("Biomarker",1:nBiom, sep="")
    
    ## Add \mu_g (baseline expression) to each row of randErr
    ## ------------------------------------------------------
    
    if (is.null(baseExpr))G<-rep(5,nBiom) else G<-baseExpr
    avgData<-randErr+matrix(G, nrow=nTrain, ncol=nBiom, byrow=T)

    
    ## Introducing differential expressions
    ##-------------------------------------
     
    if(diffExpr){
    ## Generate indicators for up/down regulation
    upDown<-sample(c(-1,1),nBiom,replace = TRUE)
        
    ## Generate log2(fold-change) from half normal distribution with mean=0, sd=sigma
    
    foldChange<-rtnorm(nBiom,mean=0,sd=sigma,lower=log2(foldMin))
    if(orderBiom) foldChange<-sort(foldChange)
    
    ## Make the data in the diseased group (D) to be up/down regulated
    ## by the amount foldChange
      
    diffD<-matrix(foldChange*upDown, ncol=nBiom, nrow=nGr2, byrow=T)
    avgData[(nGr1+1):nTrain,]<-avgData[(nGr1+1):nTrain,]+diffD
      }

    ## Prepare data for classification
    ##--------------------------------

    ## Define groups (H=Healthy, D=Diseased)

    class<-factor(c(rep("H",nGr1),rep("D",nGr2)))
    data<-data.frame(class=class,avgData)
    return(data)
}

################################################################################
## End of:                         simData.R
################################################################################

