ADAM.panessprofile<-function(depMat,display=TRUE,
                             main_suffix='fitness genes in at least 1 cell line',
                             xlab='n. dependent cell lines'){
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


        main=paste(nrow(depMat),main_suffix)
        barplot(panessprof,ylab='n.genes',xlab=xlab,cex.axis = 0.8,cex.names = 0.8,
                las=2,main=main)

        barplot(CUMsums,ylab='n.genes',xlab=xlab,cex.axis = 0.8,cex.names = 0.6,
                las=2,main='Cumulative sums')

    }
    return(list(panessprof=panessprof,CUMsums=CUMsums))
}

ADAM.generateNullModel<-function(depMat,ntrials=100,display=TRUE){

    depMat<-depMat[which(rowSums(depMat)>0),]
    nullProf<-matrix(NA,ntrials,ncol(depMat),dimnames = list(1:ntrials,1:ncol(depMat)))
    nullCumSUM<-matrix(NA,ntrials,ncol(depMat),dimnames = list(1:ntrials,paste('≥',1:ncol(depMat),sep='')))
    print('Generating null model...')
    pb <- txtProgressBar(min=1,max=ntrials,style=3)

    for (i in 1:ntrials){
        setTxtProgressBar(pb, i)
        rMat<-
            ADAM.randomisedepMat(depMat)
        Ret<-
            ADAM.panessprofile(rMat,display = FALSE)
        nullProf[i,]<-Ret$panessprof
        nullCumSUM[i,]<-Ret$CUMsums
    }
    Sys.sleep(1)
    close(pb)
    print('')
    print('Done')

    if (display){
        par(mfrow=c(2,1))
        main=c(paste(ntrials,' randomised essentiality profiles of\n',nrow(depMat),' genes across ',ncol(depMat),' cell lines',
                     sep=''))
        boxplot(nullProf,las=2,xlab='n. cell lines',ylab='genes depleted in n cell lines',main=main)
        colnames(nullCumSUM)<-paste(">=",1:ncol(nullCumSUM))
        boxplot(log10(nullCumSUM+1),las=2,main='Cumulative sums',xlab='n. cell lines',
                ylab='log10 [number of genes + 1]',
                cex.axis=0.8)
    }

    return(list(nullProf=nullProf,nullCumSUM=nullCumSUM))
}

ADAM.randomisedepMat<-function(depMat){
    rmat<-apply(depMat,2,sample)
}

ADAM.empiricalOdds<-function(observedCumSum,simulatedCumSum){

    nsamples<-length(observedCumSum)
    ntrials<-nrow(simulatedCumSum)

    odds<-rep(NA,1,nsamples)
    names(odds)<-paste('≥',1:nsamples,sep='')
    for (i in 1:nsamples){

        PDF<-density(simulatedCumSum[,i])


        odds[i]<- log10(observedCumSum[i]/mean(simulatedCumSum[,i]))

    }
    return(odds)
}

ADAM.truePositiveRate<-function(depMat,essentialGeneSet){
    nsamples<-ncol(depMat)

    essentialGeneSet<-intersect(essentialGeneSet,rownames(depMat))

    TPR<-rep(NA,1,nsamples)
    names(TPR)<-paste('≥',1:nsamples,sep='')

    ncells<-rowSums(depMat)

    TP<-rep(NA,1,nsamples)
    names(TP)<-paste('≥',1:nsamples,sep='')

    P<-rep(NA,1,nsamples)
    names(P)<-paste('≥',1:nsamples,sep='')

    for (i in nsamples:1){
        positiveSet<-names(which(ncells>=i))
        P[i]<-length(positiveSet)
        truepositives<-intersect(positiveSet,essentialGeneSet)
        TP[i]<-length(truepositives)
        TPR[i]<-TP[i]/length(essentialGeneSet)
    }

    return(list(P=P,TP=TP,TPR=TPR))
}

ADAM.tradeoffEO_TPR<-function(EO,TPR,test_set_name){
    CCOL<-'red'

    x<-EO
    x[x==Inf]<-max(x[x<Inf])
    x<-(x-min(x))/(max(x)-min(x))

    y<-TPR
    y<-(y-min(y))/(max(y)-min(y))

    orEO<-EO
    orEO[orEO==Inf]<-max(orEO[orEO<Inf])
    orTPR<-TPR

    EO<-x
    TPR<-y
    par(mar=c(4,4,4,4))
    MAIN<-c('log10 (obs/Expct) n.genes [red, left]',
            paste('% covered ',test_set_name,' [blue, right]',sep=''))
    plot(EO,type='l',xlab='genes depleted in >= # cell lines',ylab='',axes=FALSE,lwd=4,main=MAIN,col=CCOL,cex.main=0.8,
         xlim=c(0,length(EO)))
    axis(2,at = seq(0,1,0.2),format(seq(min(orEO),max(orEO),(max(orEO)-min(orEO))/5),digits=2))
    axis(1)
    par(new=TRUE)
    plot(TPR,type='l',xlab='',ylab='',axes=FALSE,lwd=4,col='blue',ylim=c(0,1),xlim=c(0,length(EO)))
    axis(4,at = seq(0,1,0.2),format(seq(min(orTPR),max(orTPR),(max(orTPR)-min(orTPR))/5),digits=2))


    point<-min(which(!y>x))

    abline(v=point)
    abline(h=y[point],lty=2)

    points(point,y[point],pch=16,cex=2)

    legend('top',paste(format(100*orTPR[point],digits=2),'% covered',sep=''),bg = NULL,bty = 'n')

    return(point)
}
