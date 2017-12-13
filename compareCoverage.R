################################################
# Script: compareCoverage.R
# Author: Bhavna Hurgobin and Michael Renton
# Email:  b.hurgobin@uq.edu.au; michael.renton@uwa.edu.au
# Usage: compareCoverage.R test.cov test
################################################

args<-commandArgs(TRUE)
filename<-args[1]

dat<-read.table(filename,sep='\t',header=T)

nrows<-dim(dat)[1]

######
#start and end of regions that are
#rnc (repeated, non-collapsed) -- r0==0 and r1 >=2 and r2 >= 0.5*overall_median ... how many times repeated
#collapsed -- mean of r0 and r1 is greater than the overall median + 2*sd

###########   DEFINE CONDITIONS HERE
#rnc 
overall_median<-median(rowMeans(dat[,3:4]))
rnc<-(dat$r0==0 & dat$r1 >=2 & dat$r2 >= 0.5*overall_median)

#collapsed 
mrds<-rowMeans(dat[,3:4])
med<-median(mrds)
mn<-mean(mrds)
sdev<-sd(mrds)
coll<-(mrds> med+2*sdev) & dat$r0>3

########## now summarise to give output table with start/end of regions; only regions at least 50 bp long are reported, this can be changed if required.
conds<-rbind(rnc ,coll) 
#any(colSums(conds)>=2)  ## check for double condition - should be FALSE
condnames<-c('rnc','coll')
nconds<-length(condnames) 
results<-NULL

inreg<-0
end<-1
for (i in 1:nrows){
	if (inreg==0 & any(conds[,i])){
		start<-i
		inreg<-which(conds[,i])
	}
	if (inreg!=0 & !any(conds[,i])){
		end<-i-1
		region.length=end-start+1
                r0<-median(dat$r0[start:end])
		r2<-median(dat$r2[start:end])
		r1<-median(dat$r1[start:end])
		r2r1rat<-round(r2/r1,3) #estimates the copy number of the region
		newrow<-c(start,end,region.length,r2,r1,r0,r2r1rat,inreg)
            if (region.length>50) results<-rbind(results,newrow)
		inreg<-0
	}
	if (inreg!=0 & any(conds[,i])){
		if (which(conds[,i])!=inreg){
			end<-i-1
			region.length<-end-start+1
			r0<-median(dat$r0[start:end])
			r2<-median(dat$r2[start:end])
			r1<-median(dat$r1[start:end])
			r2r1rat<-round(r2/r1,3) #estimates the copy number of the region
			newrow<-c(start,end,region.length,r2,r1,r0,r2r1rat,inreg)
			if (region.length>50) results<-rbind(results,newrow)
			inreg<-which(conds[,i])
			start<-i
		}
	}
}

results<-data.frame(results,row.names = NULL)
results<-cbind(1:dim(results )[1] , results )
names(results)<-c('region.number','start','end','region.length','r2','r1','r0','copy number','region.type')
results$region.type<-factor(condnames [results$region.type])
write.table(results,paste(args[2],'_all_regions.out',sep=''),sep='\t',quote=F,row.names = FALSE)
