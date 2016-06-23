
# entropy cal
# O 
rm(list=ls())
library(seqinr)
library(ade4)
library(entropy)
b1=character()
xx2b=numeric()
xx2=numeric()
 #x1=read.fasta(file="C://Victor//MATH1//myseq0.fasta")
 x1=read.fasta(file="C://Victor//MATH1//myseq6.fasta")
    gg=length(x1)
  e1=x1[[1]]
 e11=attributes(e1)
 e12=unlist(strsplit(e11[[1]],split="_"))
  cc=as.character(e12[2])  
 seqq2<-lapply(1:gg,
 FUN=function(k){
  b1=x1[[k]]
 b11=attributes(b1)
  b12=unlist(strsplit(b11[[1]],split="_"))
       startmm = as.numeric(b12[3])
	 	     b12end = as.numeric(b12[4])
		  aa2=b12end-startmm
		 		 windowsize=trunc(aa2/10.0)
		 		 if( aa2  > 1000) windowsize=trunc(0.1* aa2)
        if(aa2   >  10000) windowsize=trunc(0.01*aa2)
		if(aa2 >    100000) windowsize=trunc(0.001*aa2)
		if(windowsize < 1) windowsize=trunc(aa2/3.0)
			   b13=length(b1) - windowsize
 starts=seq(1,b13,by=windowsize)
   n=length(starts)
       seqq1<-lapply(1:n,
	 FUN=function(i){
	 		 aa1=starts[i]+windowsize -1
	 chunk=b1[starts[i]:(starts[i]+windowsize-1)]
	 	 	 	 	 aa21=count(chunk,word=3)
					 www3=entropy.empirical(aa21,unit="log2")
		 www1=starts[i]+startmm
	 www2=www1+windowsize-1 
	 	  		 		 yy1=list(www1,www2,www3) 
          return(yy1)		 
		 	 })
	  	 		xxgc1 <- data.frame(matrix(unlist(seqq1), nrow=n, byrow=T))
				if( k == 1) {
				xx2<-xxgc1
				}
				if( k > 1){
				xx2=rbind(xx2,xxgc1)
				}
				 xx2a=(xx2)
					return(xx2a)			})
					
				xx2b=do.call("rbind", seqq2)
								
		png('c://Victor//MATH1//entplot1.png')
		 plot(xx2b$X1,xx2b$X3,type="p",xlab="Nucleotide start position",ylab="entropy content",main=cc)
dev.off()		 
	  plot(xx2b$X1,xx2b$X3,type="p",xlab="Nucleotide start position",ylab="entropy content",main=cc)
  png('c://Victor//MATH1//entboxplot1.png')
		 boxplot(xx2b$X3,ylab="GC content",main=cc)
dev.off()		 
     outliersent= boxplot(xx2b$X3,ylab="entropy content",main=cc)
	head(outliersent$out)
	xxgcf=xx2b[order(xx2b$X3,decreasing=TRUE),]
		head(xxgcf,5)
	tail(xxgcf,5)

	
