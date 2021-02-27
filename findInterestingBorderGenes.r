		findInterestingBorderGenes <- function(){
		state<-readline("type control or disease.....?")
		entity="Example"
		cat('Finding intersting genes for ', state,' modules of ',entity,' expression profiles \n')
		len<-readline("number of modules....?")
		#Files1 <- sprintf("%s.csv",entity)
		Files1 <- sprintf("%s-genes.csv",entity)
		geneName=read.csv(Files1,header=FALSE)
		#geneName=as.character(matrix(geneName[1:1000,1])))
		geneName=as.character(unlist(geneName[,1]))
		Allgenes=geneName
		library(GOSemSim)
		library(mygene)
		library(hgu95av2.db) 
		library("annotate")
		borderGenesEachCluster=as.character()
		#len=as.numeric(length(FCSSSim[[1]]))
		len=as.numeric(len)
		dat <- godata('org.Hs.eg.db', ont="BP", computeIC=TRUE)	
		for (i in 1:len){
			Files1 <- sprintf("%s_%s_Module_%d.csv",state,entity,i)
			sumIndex=read.csv(Files1,header=FALSE)
			sumIndex=as.matrix(sumIndex)
			geneName=as.character(sumIndex)
			entrez=queryMany(geneName, scopes="symbol", fields="entrezgene", species="human")
			returnall=TRUE 
			geneID=entrez[,c("query","_id")]
			geneID=rename(geneID, c("query"="one", "_id"="two"))
			clusterGenes=geneID$two
			geneName=clusterGenes[!is.na(clusterGenes)]			
			genes=data.frame(geneName)	
			borderGenes=geneName			
			if(length(borderGenes)>=1000){
			borderGenes=borderGenes[1:250]
			}
		    borderGenesEachCluster=borderGenes
			#print(length(borderGenesEachCluster))
			source("getSemanticSimilarity.r")
			borderGenesForEachModule=getSemanticSimilarity(borderGenesEachCluster,genes,dat,i,state,entity)
			write.table(borderGenesForEachModule,paste0(entity,"_interestingEntrezIDs_",state,"_module_",i,".csv"), sep = ",", row.names = FALSE, col.names=FALSE)
			geneID=getGenes(borderGenesForEachModule)
				geneID=rename(geneID, c("_id"="one", "symbol"="two"))
				ht=geneID
				ht=ht[,"two"]
		
			write.table(ht,paste0(entity,"_interestingGenes_",state,"_module_",i,".csv"), sep = ",", row.names = FALSE, col.names=FALSE)
			#write.table(corrScoreMatrix,paste0(state,"-corrScoreMatrix_",i,".csv"), sep = ",", row.names = FALSE, col.names=FALSE)
			}
	}