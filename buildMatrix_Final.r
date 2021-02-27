#X-Module
#Copyright Tulika Kakati 2018 
#Distributed under GPL v3 open source license
#This code mines modules of genes using THD-Module Extractor......
buildMatrix_Final <- function(){
		print("THD-Module Extractor is running....")
		state<-readline("type control or disease.....?")
		entity="Example"
		numOfModules<-readline("Give the maximum number of modules......")
		threshold1<-readline("Give the upper bound of threshold range......")
		threshold2<-readline("Give the lower bound of threshold range......")
		threshold1=as.numeric(threshold1)
		threshold2=as.numeric(threshold2)
		thresholdU=as.numeric(threshold1)
		reduceFactor<-readline("Give the reducing factor....")
		reduceFactor=as.numeric(reduceFactor)
		min_g <- readline("Give the minimum number of genes in a cluster......")
		min_g=as.numeric(min_g)
		clusterListIndex=list()
		clusterNo=0
		Files1 <- sprintf("%s-genes.csv",entity)
		s=read.csv(Files1, header=FALSE)
		print("read row sum matrix")
		Files1 <- sprintf("scoreMatrix-%s-%s.csv",entity,state)
		matrixT=read.csv(Files1, header=FALSE)
		matrixT=as.matrix(matrixT)
		Files1 <- sprintf("sumRowMatrix-%s-%s.csv",entity,state)
		sumRowmatrixT=read.csv(Files1, header=FALSE)
		sumRowmatrixT=as.matrix(sumRowmatrixT)
		n=dim(sumRowmatrixT)[1]
		p=n
		sortedmatrixT=sort(sumRowmatrixT[1:n,],decreasing=T,index.return=T)
		print("row sum matrix is sorted..")
		dataFrameSortedmatrixT=data.frame(sortedmatrixT[[1]])
		dataFrameSortedmatrixT=dataFrameSortedmatrixT[1]
		dataFrameSortedmatrixTIndex=data.frame(sortedmatrixT[[2]])
		dataFrameSortedmatrixTIndex=dataFrameSortedmatrixTIndex[1]
		Cluster=as.integer()
		FinalCluster=list()
		FinalSSSimmatrixT=list()
		SSSimmatrixT=as.integer()
		count=1
		i=1
			while(i<n&& threshold1>=threshold2 &&count<=numOfModules){
			threshold_prod=(as.numeric(min_g)*as.numeric(threshold1))
			xx=dataFrameSortedmatrixTIndex[i,]
				while(sumRowmatrixT[xx,1]<threshold_prod && threshold1>= threshold2){
					threshold1=threshold1-reduceFactor
					threshold_prod=(as.numeric(min_g)*as.numeric(threshold1))
				}
			FinalCluster1=as.integer(unlist(FinalCluster))
			xx=dataFrameSortedmatrixTIndex[i,]
			SSSimmatrixT=matrixT[xx,xx]
			Cluster=as.integer(0)
			Cluster=xx
			for(j in 1:p){
				if (matrixT[xx,j]>=threshold1 && j!=xx && (j %in% FinalCluster1)==F){
					Cluster=append(Cluster,j)
					SSSimmatrixT=append(SSSimmatrixT,matrixT[xx,j])
					d=data.frame(Cluster)
					hh=dim(d)[1]
						for (k in 2:hh){
							xy=Cluster[k]
							for (tt in 1:p){
								if(matrixT[xy,tt]>=threshold1 && all.equal(tt,xy)==F && length(Cluster)<min_g){
									Cluster=append(Cluster,tt)
									SSSimmatrixT=append(SSSimmatrixT,matrixT[xy,tt])
								}
							}
						}
				}
			}
			Cluster=unique(unlist(Cluster))
			Cluster=list(Cluster)
			unlistCluster=unlist(Cluster)
			if (dim(as.data.frame(unlistCluster))[1]<min_g){
			a=mean(as.numeric(matrixT[xx,]))
				while(a<threshold1 && threshold1>= threshold2){
					threshold1=as.numeric(threshold1)-as.numeric(reduceFactor)
					threshold_prod=(as.numeric(min_g)*as.numeric(threshold1))
				}
				for(j in 1:p){
					if (matrixT[xx,j]>=threshold1 && j!=xx && (j %in% FinalCluster1)==F){
						Cluster=append(Cluster,j)
						SSSimmatrixT=append(SSSimmatrixT,matrixT[xx,j])
						d=data.frame(Cluster)
						hh=dim(d)[1]
						for (k in 2:hh){
						xy=as.numeric(Cluster[k])
							for (tt in 1:p){
							#print("Second Expansion")
								if(matrixT[xy,tt]>=threshold1 && all.equal(tt,xy)!=TRUE && length(Cluster)<min_g){
									Cluster=append(Cluster,tt)
									SSSimmatrixT=append(SSSimmatrixT,matrixT[xy,tt])
								}
							}
						}
					}
				}

			}
			Cluster=unique(unlist(Cluster))
			Cluster=list(Cluster[1:552])
			unlistCluster=unlist(Cluster)
			if (dim(as.data.frame(unlistCluster))[1]>=min_g){
				print("Cluster size satisfies")
				clusterGenesIndex=as.numeric(unlistCluster)
				clusterGenes=as.numeric()
				d=unlist(as.list(s))
				clusterGenes=paste(d[clusterGenesIndex])
				ht=as.character(clusterGenes)
				clusterGenes=data.frame(clusterGenes)
				count=count+1
				clusterNo=clusterNo+1
				print(dim(clusterGenes)[1])
				clusterListIndex=append(clusterListIndex,i)
				write.table(clusterGenes, paste0(state,"_",entity,"_Module_",clusterNo,".csv"),col.names = F, row.names = F)
				print("cluster found...")
				FinalCluster=append(FinalCluster,Cluster)
				FinalSSSimmatrixT=append(FinalSSSimmatrixT,list(SSSimmatrixT))
			}

			i=i+1
			threshold1=as.numeric(thresholdU)-as.numeric(reduceFactor)
			threshold_prod=(as.numeric(min_g)*as.numeric(threshold1))
			}
		clusterListIndex=unlist(clusterListIndex)
		# Files1 <- sprintf("%s/FCSSSim-%s-%s.RData",getwd(),entity,state)
		# save.image(Files1)
		# #save.image("E:\\DropBox\\Research\\5th Sem\\Module-Correspondence\\Implementation\\23.07.2018\\THD-Module Extractor\\AD\\FCSSSim-AD-control.RData")
		return(list(FinalCluster,FinalSSSimmatrixT))
}
