#! /usr/bin/Rscript

# INPUT:
#		struct_fn - RNAfold format file with the known structure
#		u - digestion params [A,C,T,G]
#		v - digestion params [A,C,T,G]
#		lower - min fragment length
#		upper - max fragment length
#		num_reads - how many simulated reads to generate
# OUTPUT:
#		*.sim.reads file


# FUNCTIONS
# generate log digestion probability matrix
logdigest_mat<-function(pairvec,u,v,lower,upper,seq) {
    # for each position, compute probability of digestion on 3' side
    k=length(pairvec)
    pm=rbind(v,u)
    # prob digestion for bond<i,i+1>
    dv = {}
    for (i in 1:k) {
    	dv=c(dv, pm[pairvec[i]+1,seq[i]])
    }
    ndv=1-dv
    ldv=log(dv)
    ldv1=c(0,ldv)
    lndv=log(ndv)
    
		# (i,j): read start from nuc i, end at nuc j
    # P[i,j]= P(digest at bond<i-1,i>) P(no digestion in between)
    ldmat=matrix(NA,k,k)
    for (i in 1:(k-1)) {
        di=c(0,cumsum(lndv[ i:(k-1) ]))
        di=di+ldv1[i]
        ldmat[i,i:k]=di
    }
    ldmat[,k]=NA
    
    # size restriction in [lower, upper]
    for (i in 1:k) {
    	for (j in 1:k) {
    		if ((j-i+1 < lower) || (j-i+1 > upper)) {
    			ldmat[i,j] = NA
    		}
    	}
    }
    
    # normalize likelihood to sum to one
    ldmat <- exp(ldmat)
    ldmat <- ldmat/sum(ldmat,na.rm=T)
    ldmat <- log(ldmat)
    
    return(ldmat)
}


# convert from digestion matrix into triple format (start, end, prob)
ldmat2triple<-function(ldmat, remove_NA=TRUE) {
    k=nrow(ldmat)
    ldmat2=as.vector(ldmat)
    imat=matrix(1:k,nrow=k,ncol=k,byrow=F)
    jmat=matrix(1:k,nrow=k,ncol=k,byrow=T)
    d2=data.frame(I=as.vector(imat),J=as.vector(jmat),LOGP=as.vector(ldmat2))
    if (remove_NA) {
    	d2=d2[!is.na(d2$LOGP),]
 		}    
 		d2$P=exp(d2$LOGP)
    return(d2)
}

# convert from triple format (start, end, prob) to digestion matrix
triple2ldmat<-function(triple) {
    k=max(triple$J)
    m=matrix(NA,k,k)
    m[as.matrix(triple[,c("I","J")])]=triple$LOGP
    return(m)
}

ldmat2abundance<-function(ldmat,type=c("5p","3p","read")) {
    dmat=exp(ldmat)
    dmat[is.na(dmat)]=0
    if (type=="5p") {
        return(rowSums(dmat))}
    if (type=="3p") {
        return(colSums(dmat))
    }
    # return read abundance
    k=nrow(dmat)
    rab=rep(0,k)
    triple=ldmat2triple(ldmat)
    for (i in 1:k) {
        rab[i]=sum(triple$P[triple$I<=i & triple$J>=i])
    }
    return(rab)
}

sample.reads<-function(triple,n) {
	options(warn=-1)
	a=table(sample(1:nrow(triple),n,repl=T,prob=triple$P))
	options(warn=0)
	a2=data.frame(triple[as.integer(names(a)),c(1,2)],CNT=as.vector(a))
	a2=a2[order(a2[,2]),]
	a2=a2[order(a2[,1]),]
	a2$LEN=a2$J-a2$I+1
	return(a2)
}

plot.ldmat<-function(ldmat) {
	heatmap.2(exp(ldmat),Rowv=NA,Colv=NA,scale="none",dendrogram="none",
		RowSideColors=pv1col,ColSideColors=pv1col,trace="none",key=F,col=greenred(100),breaks=101,
		   lhei=c(2,9),lwid=c(2,9))
}

plot.abundance<-function(ldmat) {
    layout(matrix(1:3,nrow=3))
    a=ldmat2abundance(ldmat,"5p"); barplot(a,names.arg=pv1,cex.names=0.5,main="5' endpoint abundance",ylab="Abundance")
    a=ldmat2abundance(ldmat,"3p"); barplot(a,names.arg=pv1,cex.names=0.5,main="3' endpoint abundance",ylab="Abundance")
    a=ldmat2abundance(ldmat,"read"); barplot(a,names.arg=pv1,cex.names=0.5,main="Read",ylab="Abundance")
}
# END FUNCTIONS


args <- commandArgs(T)
if (length(args) < 7) {
	cat("USAGE: ./generate_simulated_reads_v2b.R fn u v lower upper num_reads out_fn\n")
	q()
}
fn = args[1]
u = as.numeric(unlist(lapply(unlist(strsplit(args[2],",")), function(x) unlist(strsplit(x,":"))[2])))
names(u) <- unlist(lapply(unlist(strsplit(gsub("T","U",args[2]),",")), function(x) unlist(strsplit(x,":"))[1]))
v = as.numeric(unlist(lapply(unlist(strsplit(args[3],",")), function(x) unlist(strsplit(x,":"))[2])))
names(v) <- unlist(lapply(unlist(strsplit(gsub("T","U",args[3]),",")), function(x) unlist(strsplit(x,":"))[1]))
lower_size = as.numeric(args[4])
upper_size = as.numeric(args[5])
read_thres = as.numeric(args[6])
out_fn = args[7]

# struct as a numeric array: 0=dot, 1=open, 2=close
struct <- {}
pairing <- {}

# scan in the true structure
seq <- unlist(strsplit(gsub("T", "U", scan(fn, what="string", nlines=1, quiet=T)),""))
struct_str <- scan(fn, what="string", skip=1, quiet=T)
struct_str <- struct_str[1]
# convert to vector of 0s and 1s
struct <- gsub("\\.", "0", struct_str, perl=T)
struct <- gsub("\\(", "1", struct, perl=T)
struct <- gsub("\\)", "2", struct, perl=T)
pairing <- gsub("2", "1", struct, perl=T)
struct <- as.numeric(unlist(strsplit(struct,"")))
pairing <- as.numeric(unlist(strsplit(pairing,"")))
true_structure_length <- length(struct)

# generate reads
ldmat=logdigest_mat(pairing,u,v,lower_size,upper_size,seq)
ldtriple=ldmat2triple(ldmat)
ds=sample.reads(ldtriple, read_thres)
#endpoints_5p <- unlist(lapply(sort(unique(ds$I)), function(x) sum(ds[which(ds$I==x),3])))
#endpoints_3p <- unlist(lapply(sort(unique(ds$J)), function(x) sum(ds[which(ds$J==x),3])))

nrid=0
ds[,"I"] <- ds[,"I"]-1
ds[,"J"] <- ds[,"J"]-1
rstrs <- unlist(lapply(1:(dim(ds)[1]), function(x) sprintf("nr%s@@%d@@%d\t%d\t%d\t1\n",nrid+x,ds[x,"CNT"],ds[x,"LEN"],ds[x,"I"],ds[x,"J"])))
# write these reads to file
cat(rstrs, file=out_fn, sep="", append=F)


	
