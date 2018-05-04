ADAM.panessprofile<-function(depMat,display=TRUE,NLIMS=1,
                             main_suffix='genes depleted in at least 1 cell line',
                             xlab='n. cell lines'){
    depMat<-depMat[which(rowSums(depMat)>0),]
    panessprof<-rep(0,ncol(depMat))
    names(panessprof)<-as.character(1:ncol(depMat))
    paness<-summary(as.factor(rowSums(depMat)),maxsum = length(unique(as.factor(rowSums(depMat)))))
    panessprof[as.character(names(paness))]<-paness

    CUMsums<-rev(cumsum(rev(panessprof)))

    names(CUMsums)<-paste('>=',names(CUMsums),sep='')

    if(display){
        par(mfrow=c(2,1))
        par(mar=c(6,4,4,1))

        if(length(main)==0){
            main=c(nrow(depMat),main_suffix)
        }
        barplot(panessprof,ylab='n.genes',xlab=xlab,cex.axis = 0.8,cex.names = 0.8,
                las=2,main=main)

        barplot(CUMsums,ylab='n.genes',xlab=xlab,cex.axis = 0.8,cex.names = 0.6,
                las=2,main='Cumulative sums')

    }
    return(list(panessprof=panessprof,CUMsums=CUMsums))
}
