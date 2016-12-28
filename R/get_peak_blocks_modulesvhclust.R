get_peak_blocks_modulesvhclust <-
function(dataA,simmat,adjacencyfromsimilarity=FALSE,time_step=3,max.rt.diff=10,outloc,column.rm.index=NA,cor.thresh=NA,deepsplit=2,
minclustsize=20,cutheight=0.2,cormethod="spearman",networktype = "unsigned",num_nodes=2,step1log2scale=FALSE,hclustmethod="complete")
{

# browser()
cnames<-colnames(dataA)

cnames[1]<-"mz"
cnames[2]<-"time"

colnames(dataA)<-as.character(cnames)

data_mzrt<-dataA[,c(1:2)]


#1:3,6:11
if(is.na(column.rm.index)==FALSE){
dataA<-dataA[,-c(column.rm.index)]
}


feat_inf<-paste(dataA[,1],dataA[,2],sep="_")
dataA<-dataA[,-c(1:2)]

data_m<-t(dataA)

if(step1log2scale==TRUE){

		data_m<-2^(data_m)
}


sample.col.start<-2

 #which(dataA$mz>145.9 & dataA$mz<146)
 #cor.test(t(as.numeric(dataA[303,-c(1:2)])),t(as.numeric(dataA[774,-c(1:2)])))

powers = c(c(1:10), seq(from = 12, to=20, by=2))



if(adjacencyfromsimilarity==FALSE)
{
	sft = pickSoftThreshold(data=data_m, dataIsExpr=TRUE,powerVector = powers, verbose = 0)
power_val=sft$powerEstimate

if(is.na(power_val)==TRUE){
power_val=6
}

if(cormethod=="pearson"){
ADJdataOne<-adjacency(datExpr=data_m,
                         type = networktype,
                         power = power_val,corOptions = "use = 'p'")
                         }else{
                         	ADJdataOne<-adjacency(datExpr=data_m,
                         type = networktype,
                         power = power_val,corOptions = "use = 'p',method='spearman'")
                         	}
}else{

sft = sft = pickSoftThreshold.fromSimilarity(similarity=simmat, powerVector = powers, verbose = 0)
power_val=sft$powerEstimate

if(is.na(power_val)==TRUE){
power_val=6
}

ADJdataOne<-adjacency.fromSimilarity(similarity=simmat,power=power_val,type=networktype)
}

rnames_simmat<-rownames(ADJdataOne)

dup_names<-which(duplicated(rnames_simmat)==TRUE)

if(length(dup_names)>0){

ADJdataOne<-ADJdataOne[-c(dup_names),-c(dup_names)]
}

#dissTOMCormat=TOMdist(cormat,TOMType="signed")

dissTOMCormat=TOMdist(ADJdataOne)

#if(FALSE)
{
hierTOMCormat = flashClust(as.dist(dissTOMCormat),method=hclustmethod);
par(mfrow=c(1,1))
pdf("plot.pdf")
plot(hierTOMCormat,labels=F,main="Dendrogram")
dev.off()

#save(list=ls(),file="hier_analysis.Rda")

#if(FALSE)
{


colorhdataOne2=cutreeDynamic(hierTOMCormat,distM= dissTOMCormat,deepSplit=deepsplit, minClusterSize=minclustsize, pamRespectsDendro = FALSE, pamStage=TRUE)

#colorhdataOne2=cutreeHybrid(hierTOMCormat,distM= dissTOMCormat,deepSplit=deepsplit, minClusterSize=minclustsize, pamRespectsDendro = FALSE,pamStage=TRUE)

l2colors<-levels(as.factor(colorhdataOne2)) #$labels))


if(length(l2colors)>1){

m2=try(mergeCloseModules(data_m,colors=colorhdataOne2,cutHeight=cutheight),silent=TRUE)

if(is(m2,"try-error")){

	#m2<-doWGCNA(data_m,deepsplit=deepsplit,minclustsize=minclustsize,cutheight=cutheight)
	mod_list<-colorhdataOne2 #as.numeric(m2$colors)
}else{

	mod_list<-as.numeric(m2$colors)
}
#m2=mergeCloseModules(data_m,colors=colorhdataOne2$labels,cutHeight=cutheight)

#print(length(table(m2$colors)))

}else{
	#m2=colorhdataOne2
	#mod_list<-m2 #as.numeric(m2)
	#mod_list<-rep(0,dim(dataA)[1])
	m2<-doWGCNA(data_m,deepsplit=deepsplit,minclustsize=minclustsize,cutheight=cutheight)
	mod_list<-as.numeric(m2$colors)
}
#max.rt.diff<-5


 #mod_list<-as.numeric(m2$colors)
}

}

Alldegrees1=intramodularConnectivity(ADJdataOne,mod_list)

#rownames(Alldegrees1)<-as.character(feat_inf)
#print("mod_list")
#print(mod_list[1:4])

#print("here")

t1<-table(mod_list)
mod_names<-names(t1)
mod_names<-as.numeric(mod_names)

#print(dim(t1))

time_mult_fact<-1

#print("mod_names")
#print(mod_names[1:4])


 diffmatC<-{}
rm(data_m)


dataA<-cbind(data_mzrt,dataA)
dataA<-as.data.frame(dataA)

 #d1<-density(dataA$time,bw="nrd",from=min(dataA$time),to=(max.rt.diff+max(dataA$time,na.rm=TRUE)))

 #time_step<-abs(d1$x[1]-d1$x[time_step+1])

 #t1<-d1$x

 #time_step<-1.5*time_step
 #save(list=ls(),file="current.Rda")


diffmatB<-{}
diffmatB<-lapply(1:length(mod_names),function(i){

	groupA_num<-mod_names[i]

	subdata<-dataA[which(mod_list==groupA_num),]
	subdata<-subdata[order(subdata$time),]

 	groupB<-group_by_rt_hist(subdata,time_step,max_diff_rt=max.rt.diff,groupnum=groupA_num)

	rownames(groupB)<-NULL


	return(groupB)
	  })

	#print("here2")

    diffmatB<-ldply(diffmatB,rbind)
   # save(diffmatB,file="diffmatB.Rda")

	rm(dataA)

#rownames(diffmatB)<-as.character(feat_inf)
	diffmatB<-cbind(Alldegrees1[,c(1:4)],diffmatB)

	 return(diffmatB)
}
