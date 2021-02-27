moduleCorrespondence <- function(){
#setwd("E:/Research/5th Sem/Module-Correspondence/R package")
library(mygene)
library(reshape2)
library(dplyr)
library(igraph)
library(WGCNA)
library(GOSim)
entity="Example"
nc<- readline(prompt="number of modules for control state: ")
nd<- readline(prompt="Enter number of modules for disease state: ")
# w1=as.numeric(readline(prompt="Enter matching threshold for semantic similarity: "))
# w2=as.numeric(readline(prompt="Enter matching threshold for degree difference: "))
# w3=as.numeric(readline(prompt="Enter matching threshold for correlation of eigenegene module similarity: "))
# w4=as.numeric(readline(prompt="Enter matching threshold for causal gene similarity: "))
hubGeneIndexControlList=list()
hubGeneIndexDiseaseList=list()
expressionMatrixOfControlModuleList=list() 
expressionMatrixOfDiseaseModuleList=list()
geneNameControlModuleList=list()
geneNameDiseaseModuleList=list()
termListOfControlModuleList=list()
termListOfDiseaseModuleList=list()
neighborIndexOfHubGeneOfControlModuleList=list()
neighborIndexOfHubGeneOfDiseaseModuleList=list()
avgSemanticSimilarityOfTermsOfControlModuleList=list()
avgSemanticSimilarityOfTermsOfDiseaseModuleList=list()
eigeneModuleControlList=list()
eigeneModuleDiseaseList=list()
neighborScoreControlModuleList=list()
neighborScoreDiseaseModuleList=list()
hubGeneListMatrixControlDisease=list()
semanticDiffControlDisease=as.numeric()
degreeDiffControlDisease=as.numeric()
corrEigenControlDisease=as.numeric()
commonGeneDiffControlDisease=as.numeric()
infoAboutControlDiseaseModuleList=as.numeric()
controlModuleWithADGenes=as.numeric()
diseaseModuleWithADGenes=as.numeric()
# Data = read.csv("ADCausalGenes.csv",header=FALSE);
# ADCausalGenes=Data[,1]
# ADCausalGenes=as.character(ADCausalGenes)
j=1
for (i in 1:nc){
	 FileC <- sprintf("%s_interestingGenes_control_module_%i.csv",entity, i)
	 read.csv(FileC, header=FALSE)->moduleGenesC
	 moduleGenesC=as.character(moduleGenesC$V1)
	 write.table(moduleGenesC, paste0("filteredModuleControlGenes_",j,".csv"),col.names = F, row.names = F)
	 j=j+1
	 #commonModuleGenes=length(intersect(moduleGenesC,ADCausalGenes))
	 # if (length(intersect(moduleGenesC,ADCausalGenes))>=5)
	 # {
		# controlModuleWithADGenes[length(controlModuleWithADGenes)+1]=i
		# write.table(moduleGenesC, paste0("filteredModuleControlGenes_",j,".csv"),col.names = F, row.names = F)
		# j=j+1
	 # }
 }
 #controlModuleWithADGenes=unique(controlModuleWithADGenes)
 j=1
 for (i in 1:nd){
	 FileC <- sprintf("%s_interestingGenes_disease_module_%i.csv",entity, i)
	 read.csv(FileC, header=FALSE)->moduleGenesC
	 moduleGenesC=as.character(moduleGenesC$V1)
	 write.table(moduleGenesC, paste0("filteredModuleDiseaseGenes_",j,".csv"),col.names = F, row.names = F)
	 j=j+1
	 # commonModuleGenes=length(intersect(moduleGenesC,ADCausalGenes))
	 # if (length(intersect(moduleGenesC,ADCausalGenes))>=5)
	 # {
		# diseaseModuleWithADGenes[length(diseaseModuleWithADGenes)+1]=i
		# write.table(moduleGenesC, paste0("filteredModuleDiseaseGenes_",j,".csv"),col.names = F, row.names = F)
		# j=j+1
	 # }
 }
 #diseaseModuleWithADGenes=unique(diseaseModuleWithADGenes)
 Files1 <- sprintf("Example-control.csv")
read.csv(Files1,header=FALSE)->Data
Files1 <- sprintf("%s-genes.csv",entity)
read.csv(Files1,header=FALSE)->DEG
genesNames=DEG[,1]
genesNames=as.character(genesNames)

 #j=1
 #for (i in 1:length(controlModuleWithADGenes)){
for (i in 1:nc){
	 #j=as.numeric(controlModuleWithADGenes[i])
	 #print(i)
	 FileC <- sprintf("filteredModuleControlGenes_%i.csv", i)
	 read.csv(FileC, header=FALSE)->moduleGenesC
	 moduleGenesC=as.character(moduleGenesC$V1)
	 commonModuleGenes=intersect(moduleGenesC,genesNames)
	 commonModuleGenesIndices=match(commonModuleGenes,genesNames) 
	 expressionMatrixOfModule=Data[commonModuleGenesIndices,2:dim(Data)[2]]
	 geneNameControlModuleList[[i]]=moduleGenesC
	 expressionMatrixOfModule=data.matrix(expressionMatrixOfModule)
	 expressionMatrixOfControlModuleList[[i]]=expressionMatrixOfModule
	 #j=j+1
}
Files1 <- sprintf("Example-disease.csv")
#Files1 <- sprintf("%s_DEGs_Disease_Stage.csv",entity)
read.csv(Files1,header=FALSE)->Data
genesNames=genesNames
 
# j=1
 #for (i in 1:length(diseaseModuleWithADGenes)){
 for (i in 1:nd){
	 #j=as.numeric(diseaseModuleWithADGenes[i])
	 FileD <- sprintf("filteredModuleDiseaseGenes_%i.csv", i)
	 read.csv(FileD, header=FALSE)->moduleGenesD
	 moduleGenesD=as.character(moduleGenesD$V1)
	 commonModuleGenes=intersect(moduleGenesD,genesNames)
	 commonModuleGenesIndices=match(commonModuleGenes,genesNames) 
	 expressionMatrixOfModule=Data[commonModuleGenesIndices,2:dim(Data)[2]]
	 geneNameDiseaseModuleList[[i]]=moduleGenesD
	 expressionMatrixOfModule=data.matrix(expressionMatrixOfModule)
	 expressionMatrixOfDiseaseModuleList[[i]]=expressionMatrixOfModule 
	 #j=j+1
 }
  for (i in 1:length(expressionMatrixOfControlModuleList)){
 	FileC <- sprintf("filteredModuleControlGenes_%i.csv", i)
	read.csv(FileC, header=FALSE)->moduleGenesNamesC
	expressionMatrixOfControlModuleList[[i]]->moduleGenesC
	rownames(moduleGenesC)=NULL
	moduleGenesC=t(moduleGenesC)
	dim(moduleGenesC)
	my_cor_matrix <- cor(moduleGenesC)
	dim(my_cor_matrix)
	#write.table(my_cor_matrix, paste0("control_cor_matrix", i,".txt"),col.names = F, row.names = T)
	my_cor_matrix[upper.tri(my_cor_matrix)] <- 42  
	my_cor_df <- melt(my_cor_matrix)
	head(my_cor_df)
	dim(my_cor_df)   
	my_cor_df <- filter(my_cor_df, value != 42) %>% filter(Var1 != Var2)
	dim(my_cor_df)
	head(my_cor_df)
	my_adj_list <- my_cor_df %>% filter(value > 0.5)
	names(my_adj_list) <- c('from', 'to', 'weight')
	dim(my_adj_list)
	my_adj_list=my_adj_list[1:dim(my_adj_list)[1],1:2]
	my_adj_matrix=get.adjacency(graph.edgelist(as.matrix(my_adj_list), directed=FALSE))
	my_adj_matrix=data.matrix(my_adj_matrix)
	#write.table(my_adj_matrix, paste0("control_adj_matrix", i,".txt"),col.names = F, row.names = T)
	colnames(my_adj_matrix)=paste("V",1:dim(my_adj_matrix)[1])
	mm=rowSums(my_adj_matrix)
	ind=which(mm==max(mm)) 
	hubGeneIndexControlList[[i]]=ind
	neighbor=as.numeric()
	for (d in 1: length(ind)){
		hubGeneRow=my_adj_matrix[ind[d],]
		colnames(hubGeneRow)=NULL
		neighbor= append(neighbor,as.numeric(which(hubGeneRow!=0,arr.ind=T)))
	}
	neighbor=data.matrix(neighbor)
	rownames(neighbor)=NULL
	colnames(neighbor)=paste("V")
	neighborIndexOfHubGeneOfControlModuleList[[i]]=neighbor
	hubGene=moduleGenesNamesC$V1[ind]
	geneList=c(neighbor,ind)
	hubGeneNeighborList=moduleGenesNamesC$V1[geneList]
	res1 <- queryMany(hubGeneNeighborList, scopes='symbol', fields=c('entrezgene', 'go'), species='human')
	cc=res1[1, 'go.CC'][[1]]$id
	mf=res1[1, 'go.MF'][[1]]$id
	bp=res1[1, 'go.BP'][[1]]$id
	terms1=as.character()
	terms1=c(cc,mf,bp)
	termListOfControlModuleList[[i]]=terms1
	neighborScoreControlModuleList[[i]]=length(neighbor)/dim(moduleGenesC)[1]
	if (length(terms1)>0){
		result=getTermSim(terms1, method = "Resnik", verbose = FALSE)
		sumR=as.numeric(rowSums(result, na.rm=TRUE))
		sumR[!is.finite(sumR)]=0
		avgSemanticSimilarityOfTerms=sum(sumR)/dim(result)[1]
	}
	if(length(terms1)==0){
		result=0
		avgSemanticSimilarityOfTerms=0
	}
	print(avgSemanticSimilarityOfTerms)
	avgSemanticSimilarityOfTermsOfControlModuleList[[i]]=as.numeric(avgSemanticSimilarityOfTerms)
	datME=moduleEigengenes(moduleGenesC,1:dim(moduleGenesC)[2])$eigengene
	eigengeneModule=datME[,1]
	eigeneModuleControlList[[i]]=eigengeneModule
 }
  for (i in 1:length(expressionMatrixOfDiseaseModuleList)){
	 FileC <- sprintf("filteredModuleDiseaseGenes_%i.csv", i)
	 read.csv(FileC, header=FALSE)->moduleGenesNamesD
	 expressionMatrixOfDiseaseModuleList[[i]]->moduleGenesD
	 rownames(moduleGenesD)=NULL
	 moduleGenesD=t(moduleGenesD)
	 dim(moduleGenesD)
	 my_cor_matrix <- cor(moduleGenesD)
	 dim(my_cor_matrix)
	 #write.table(my_cor_matrix, paste0("disease_cor_matrix", i,".txt"),col.names = F, row.names = T)
	 my_cor_matrix[upper.tri(my_cor_matrix)] <- 42  
	 my_cor_df <- melt(my_cor_matrix)
	 head(my_cor_df)
	 dim(my_cor_df)   
	 my_cor_df <- filter(my_cor_df, value != 42) %>% filter(Var1 != Var2)
	 dim(my_cor_df)
	 head(my_cor_df)
	 my_adj_list <- my_cor_df %>% filter(value > 0.5)
	 names(my_adj_list) <- c('from', 'to', 'weight')
	 dim(my_adj_list)
	 my_adj_list=my_adj_list[1:dim(my_adj_list)[1],1:2]
	 my_adj_matrix=get.adjacency(graph.edgelist(as.matrix(my_adj_list), directed=FALSE))
	 my_adj_matrix=data.matrix(my_adj_matrix)
	 #write.table(my_adj_matrix, paste0("disease_adj_matrix", i,".txt"),col.names = F, row.names = T)
	 colnames(my_adj_matrix)=paste("V",1:dim(my_adj_matrix)[1])
	 mm=rowSums(my_adj_matrix)
	 ind=which(mm==max(mm)) 
	 hubGeneIndexDiseaseList[[i]]=ind
	 neighbor=as.numeric()
	 for (d in 1: length(ind)){
		hubGeneRow=my_adj_matrix[ind[d],]
		colnames(hubGeneRow)=NULL
		neighbor= append(neighbor,as.numeric(which(hubGeneRow!=0,arr.ind=T)))
	 }
	neighbor=data.matrix(neighbor)
	rownames(neighbor)=NULL
	colnames(neighbor)=paste("V")
	neighborIndexOfHubGeneOfDiseaseModuleList[[i]]=neighbor
	hubGene=moduleGenesNamesD$V1[ind]
	geneList=c(neighbor,ind)
	hubGeneNeighborList=moduleGenesNamesD$V1[geneList]
	res1 <- queryMany(hubGeneNeighborList, scopes='symbol', fields=c('entrezgene', 'go'), species='human')
	#print(res1)
	cc=res1[1, 'go.CC'][[1]]$id
	mf=res1[1, 'go.MF'][[1]]$id
	bp=res1[1, 'go.BP'][[1]]$id
	terms1=as.character()
	terms1=c(cc,mf,bp)
	termListOfDiseaseModuleList[[i]]=terms1
	neighborScoreDiseaseModuleList[[i]]=length(neighbor)/dim(moduleGenesD)[1]
	if (length(terms1)>0){
		result=getTermSim(terms1, method = "Resnik", verbose = FALSE)
		sumR=as.numeric(rowSums(result, na.rm=TRUE))
		sumR[!is.finite(sumR)]=0
		avgSemanticSimilarityOfTerms=sum(sumR)/dim(result)[1]
	}
	if(length(terms1)==0){
		result=0
		avgSemanticSimilarityOfTerms=0
	}
	avgSemanticSimilarityOfTermsOfDiseaseModuleList[[i]]=as.numeric(avgSemanticSimilarityOfTerms)
	datME=moduleEigengenes(moduleGenesD,1:dim(moduleGenesD)[2])$eigengene
	eigengeneModule=datME[,1]
	eigeneModuleDiseaseList[[i]]=eigengeneModule		
 }
 moduleSimilarityMatrix=matrix(nrow=as.numeric(length(expressionMatrixOfControlModuleList)),ncol=as.numeric(length(expressionMatrixOfDiseaseModuleList)))
 for (i in 1:length(expressionMatrixOfControlModuleList)){
	for (j in 1:length(expressionMatrixOfDiseaseModuleList)){
		semanticDiff=as.numeric()
		degreeDiff=as.numeric()
		degreeDiff=as.numeric()
		mmControl=matrix()
		mmDisease=matrix()
		corrEigene=as.numeric()
		commonGenesNumber=as.numeric()
		dataC=expressionMatrixOfControlModuleList[[i]]
		dataD=expressionMatrixOfDiseaseModuleList[[j]]
		sumR=0
	for (d in 1:length(hubGeneIndexControlList[[i]])){
		dd=as.numeric(dataC[hubGeneIndexControlList[[i]][d],])
		mmControl[d]=cor(as.numeric(eigeneModuleControlList[[i]]),dd)
	}
	mmControlValue=max(mmControl)
	dd=mmControl
	hubGeneIndexControlList[[i]]=hubGeneIndexControlList[[i]][which(dd==max(dd))] #the hub list is replaced by the hubgene which hass max correlation with the eigengene of the module
	sumR=0
	for (d in 1:length(hubGeneIndexDiseaseList[[j]])){
		dd=as.numeric(dataD[hubGeneIndexDiseaseList[[j]][d],])
		mmDisease[d]=cor(as.numeric(eigeneModuleDiseaseList[[j]]),dd)
	}
	mmDiseaseValue=max(mmDisease)
	dd=mmDisease
	hubGeneIndexDiseaseList[[j]]=hubGeneIndexDiseaseList[[j]][which(dd==max(dd))]#the hub list is replaced by the hubgene which hass max correlation with the eigengene of the module
			a1=as.numeric()
			a2=as.numeric()
			a=as.numeric()
			b=matrix()
			a=as.numeric(avgSemanticSimilarityOfTermsOfControlModuleList[[i]])
			b=as.matrix(avgSemanticSimilarityOfTermsOfControlModuleList)
			b=as.numeric(b)
			a=(a-min(b))/(max(b)-min(b))
			a1=as.numeric(a)
			a=as.numeric()
			b=matrix()
			a=as.numeric(avgSemanticSimilarityOfTermsOfDiseaseModuleList[[j]])
			b=as.matrix(avgSemanticSimilarityOfTermsOfDiseaseModuleList)
			b=as.numeric(b)
			a=(a-min(b))/(max(b)-min(b))
			a2=as.numeric(a)
			semanticDiff=1/((abs(a1-a2))+1)
			a1=as.numeric()
			a2=as.numeric()
			a=as.numeric()
			b=matrix()
			a=as.numeric(neighborScoreControlModuleList[[i]])
			b=as.matrix(neighborScoreControlModuleList)
			b=as.numeric(b)
			a=(a-min(b))/(max(b)-min(b))
			a1=as.numeric(a)
			a=as.numeric()
			b=matrix()
			a=as.numeric(neighborScoreDiseaseModuleList[[j]])
			b=as.matrix(neighborScoreDiseaseModuleList)
			b=as.numeric(b)
			a=(a-min(b))/(max(b)-min(b))
			a2=as.numeric(a)		
			degreeDiff=1/(abs(a1-a2)+1)
			corrEigene=1/(abs(as.numeric(mmControlValue-mmDiseaseValue))+1)
			a1=as.numeric()
			a2=as.numeric()
			a=as.numeric()
			b=matrix()
			a=as.numeric(length(intersect(as.character(geneNameControlModuleList[[i]]),as.character(geneNameDiseaseModuleList[[j]])))/((as.numeric(dim(expressionMatrixOfControlModuleList[[i]])[1]))+(as.numeric(dim(expressionMatrixOfDiseaseModuleList[[j]])[1]))))
			commonGeneNameControlDiseaseModule=matrix(nrow=as.numeric(length(expressionMatrixOfControlModuleList)),ncol=as.numeric(length(expressionMatrixOfDiseaseModuleList)))
			for (p in 1:length(expressionMatrixOfControlModuleList)){
			for (q in 1:length(expressionMatrixOfDiseaseModuleList)){
			commonGenesNumber=as.numeric(length(intersect(as.character(geneNameControlModuleList[[p]]),as.character(geneNameDiseaseModuleList[[q]])))/((as.numeric(dim(expressionMatrixOfControlModuleList[[p]])[1]))+(as.numeric(dim(expressionMatrixOfDiseaseModuleList[[q]])[1]))))
			commonGeneNameControlDiseaseModule[p,q]<-commonGenesNumber
			}
			}
			b=as.numeric(commonGeneNameControlDiseaseModule)
			a=(a-min(b))/(max(b)-min(b))
			a1=as.numeric(a)
			commonGenesNumber=a1
			semanticDiffControlDisease[length(semanticDiffControlDisease)+1]=semanticDiff
			degreeDiffControlDisease[length(degreeDiffControlDisease)+1]=degreeDiff
			corrEigenControlDisease[length(corrEigenControlDisease)+1]=corrEigene
			commonGeneDiffControlDisease[length(commonGeneDiffControlDisease)+1]=commonGenesNumber
			#res=1
			# source("getDecision.r")
			# res=getDecision(semanticDiff,degreeDiff,corrEigene,commonGenesNumber,w1,w2,w3,w4,i,j)
			# if (res==1){
			# infoAboutControlDiseaseModuleList=as.numeric(semanticDiff+degreeDiff+corrEigene+commonGenesNumber)
			# moduleSimilarityMatrix[i,j]<-infoAboutControlDiseaseModuleList
			
			# }
			# else{
			# infoAboutControlDiseaseModuleList=-1
			# moduleSimilarityMatrix[i,j]<-infoAboutControlDiseaseModuleList			
			
			# }
				
	}
}
  
	semanticDiffControlDisease=matrix(semanticDiffControlDisease,nrow=length(expressionMatrixOfControlModuleList), ncol=length(expressionMatrixOfDiseaseModuleList))
	degreeDiffControlDisease=matrix(degreeDiffControlDisease,nrow=length(expressionMatrixOfControlModuleList), ncol=length(expressionMatrixOfDiseaseModuleList))
	corrEigenControlDisease=matrix(corrEigenControlDisease,nrow=length(expressionMatrixOfControlModuleList), ncol=length(expressionMatrixOfDiseaseModuleList))
	commonGeneDiffControlDisease=matrix(commonGeneDiffControlDisease,nrow=length(expressionMatrixOfControlModuleList), ncol=length(expressionMatrixOfDiseaseModuleList))
	write.table(semanticDiffControlDisease,paste0("semanticDiffControlDisease_",length(expressionMatrixOfControlModuleList),"_",length(expressionMatrixOfDiseaseModuleList),".csv"), sep = ",", row.names = FALSE, col.names=FALSE)
	write.table(degreeDiffControlDisease,paste0("degreeDiffControlDisease_",length(expressionMatrixOfControlModuleList),"_",length(expressionMatrixOfDiseaseModuleList),".csv"), sep = ",", row.names = FALSE, col.names=FALSE)
	write.table(corrEigenControlDisease,paste0("corrEigenControlDisease_",length(expressionMatrixOfControlModuleList),"_",length(expressionMatrixOfDiseaseModuleList),".csv"), sep = ",", row.names = FALSE, col.names=FALSE)
	write.table(commonGeneDiffControlDisease,paste0("commonGeneDiffControlDisease_",length(expressionMatrixOfControlModuleList),"_",length(expressionMatrixOfDiseaseModuleList),".csv"), sep = ",", row.names = FALSE, col.names=FALSE)	
	w1=mean(semanticDiffControlDisease)
	#print(w1)
	w2=mean(degreeDiffControlDisease)
	#print(w2)
	w3=mean(corrEigenControlDisease)
	#print(w3)
	w4=mean(commonGeneDiffControlDisease)
	# print(w4)
	# print(dim(semanticDiffControlDisease))
	# print(dim(degreeDiffControlDisease))
	# print(dim(corrEigenControlDisease))
	# print(dim(commonGeneDiffControlDisease))
	source("getDecision.r")
for (i in 1:length(expressionMatrixOfControlModuleList)){
	for (j in 1:length(expressionMatrixOfDiseaseModuleList)){
		semanticDiff=as.numeric()
		degreeDiff=as.numeric()
		corrEigene=as.numeric()
		commonGenesNumber=as.numeric()
		semanticDiff=as.numeric(semanticDiffControlDisease[i,j])
		#print(semanticDiff)		
		degreeDiff=as.numeric(degreeDiffControlDisease[i,j])
		#print(degreeDiff)		
		corrEigene=as.numeric(corrEigenControlDisease[i,j])
		#print(corrEigene)		
		commonGenesNumber=as.numeric(commonGeneDiffControlDisease[i,j])
		#print(commonGenesNumber)		
		res=0
	
			res=getDecision(semanticDiff,degreeDiff,corrEigene,commonGenesNumber,w1,w2,w3,w4,i,j)
			#print("back")
			#print(res)
			if (res==1){
			infoAboutControlDiseaseModuleList=as.numeric(semanticDiff+degreeDiff+corrEigene+commonGenesNumber)
			moduleSimilarityMatrix[i,j]<-infoAboutControlDiseaseModuleList
			
			}
			else{
			infoAboutControlDiseaseModuleList=-1
			moduleSimilarityMatrix[i,j]<-infoAboutControlDiseaseModuleList			
			
			}
			}
			}
	
	write.table(moduleSimilarityMatrix,paste0("moduleSimilarityMatrix_",length(expressionMatrixOfControlModuleList),"_",length(expressionMatrixOfDiseaseModuleList),".csv"), sep = ",", row.names = FALSE, col.names=FALSE)
	ControlModule=as.numeric()
	DiseaseModule=as.numeric()
	moduleCorrespondence=matrix()
	for (i in 1:dim(moduleSimilarityMatrix)[1]){
		ControlModule[length(ControlModule)+1]=i
		dd=as.numeric(moduleSimilarityMatrix[i,1:dim(moduleSimilarityMatrix)[2]])
		ind=which(dd==max(dd))
		DiseaseModule[length(DiseaseModule)+1]=ind
	}
	moduleCorrespondence=cbind(ControlModule,DiseaseModule)
	write.table(moduleCorrespondence,paste0("moduleCorrespondence_",nc,"_",nd,".csv"), sep = ",", row.names = FALSE, col.names=FALSE)
	cat("Check output file", paste0("moduleCorrespondence_",nc,"_",nd,".csv"))
	
}	