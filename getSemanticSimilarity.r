#D-CoMEx
#Copyright Tulika Kakati 2018 
#Distributed under GPL v3 open source license
getSemanticSimilarity <- function(clusterGenes,s,dat,d,state,entity){
		library(mygene)
		library("org.Hs.eg.db")
		entrezIDs=as.character()
		genes=s
		#genes=as.character(unlist(genes[,1]))
		print("semantic similarity for entrez IDs....")
		sim<-mgeneSim(clusterGenes,semData=dat,measure="Resnik")
			if(length(sim)>3){
			xx=data.frame(colnames(sim))
			xx=unlist(xx,use.names=FALSE)
			entrezNames=unlist(genes,use.names=FALSE)
			geneIndex=match(xx,entrezNames)
			geneIndex=geneIndex[!is.na(geneIndex)]
			genesNames=as.character(genes[geneIndex,])
			print(length(genesNames))
			sum=sim[1,1]
			for (i in 2: dim(sim)[1]){
				sum=sum+sim[i,1]
			}
			threshold=sum/dim(sim)[1]
			for (i in 2:dim(sim)[1]){
				if (sim[i,1]< threshold){
					genesNames=genesNames[-i]
				}
			}
			clusterWithSemanticSimilarity=genesNames
			#print(length(genesNames))
			write.table(clusterWithSemanticSimilarity, paste0(state,"_",entity,"ClusterWithSemanticSimilarity_",d,".csv"),col.names = F, row.names = F)
			geneID=data.frame(genesNames)
		}
		if (length(sim)<=1){
			clusterWithSemanticSimilarity=NA
		}
		return(clusterWithSemanticSimilarity)
}

