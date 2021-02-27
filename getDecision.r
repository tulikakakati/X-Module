res=getDecision<-function(semanticDiff,degreeDiff,corrEigene,commonGenesNumber,w1,w2,w3,w4,i,j){
decisionList=as.numeric()
decisionMatrix=matrix(nrow=1,ncol=5)
res=as.numeric(0)
if (semanticDiff>=w1 && degreeDiff>=w2 && corrEigene>=w3 && commonGenesNumber>=w4)
{
res=1
decisionList=c(1,1,1,1,res)
}
if (semanticDiff>=w1 && degreeDiff>=w2 && corrEigene>=w3 && commonGenesNumber<w4)
{
res=1
decisionList=c(1,1,1,0,res)
}
if (semanticDiff>=w1 && degreeDiff>=w2 && corrEigene<w3 && commonGenesNumber>=w4)
{
res=1
decisionList=c(1,1,0,1,res)
}
if (semanticDiff>=w1 && degreeDiff<w2 && corrEigene>=w3 && commonGenesNumber>=w4)
{
res=1
decisionList=c(1,0,1,1,res)
}
if (semanticDiff<w1 && degreeDiff>=w2 && corrEigene>=w3 && commonGenesNumber>=w4)
{
res=1
decisionList=c(0,1,1,1,res)
}
if (semanticDiff>=w1 && degreeDiff>=w2 && corrEigene<w3 && commonGenesNumber<w4)
{
res=1
decisionList=c(1,1,0,0,res)
}
if (semanticDiff>=w1 && degreeDiff<w2 && corrEigene>=w3 && commonGenesNumber<w4)
{
res=1
decisionList=c(1,0,1,0,res)
}
if (semanticDiff<w1 && degreeDiff>=w2 && corrEigene>=w3 && commonGenesNumber<w4)
{
res=1
decisionList=c(0,1,1,0,res)
}
if (semanticDiff>=w1 && degreeDiff<w2 && corrEigene<w3 && commonGenesNumber<w4)
{
res=0
decisionList=c(1,0,0,0,res)
}
if (semanticDiff<w1 && degreeDiff>=w2 && corrEigene<w3 && commonGenesNumber>=w4)
{
res=0
decisionList=c(0,1,0,1,res)
}
if (semanticDiff>=w1 && degreeDiff<w2 && corrEigene<w3 && commonGenesNumber<w4)
{
res=0
decisionList=c(1,0,0,0,res)
}
if (semanticDiff<w1 && degreeDiff<w2 && corrEigene<w3 && commonGenesNumber>=w4)
{
res=0
decisionList=c(0,0,0,1,res)
}
if (semanticDiff<w1 && degreeDiff>=w2 && corrEigene<w3 && commonGenesNumber<w4)
{
res=0
decisionList=c(0,1,0,0,res)
}
if (semanticDiff<w1 && degreeDiff<w2 && corrEigene>=w3 && commonGenesNumber<w4)
{
res=0
decisionList=c(0,0,1,0,res)
}
if (semanticDiff<w1 && degreeDiff<w2 && corrEigene<w3 && commonGenesNumber<w4)
{
res=0
decisionList=c(0,0,0,0,res)
}
if (semanticDiff<w1 && degreeDiff<w2 && corrEigene>=w3 && commonGenesNumber>=w4)
{
res=0
decisionList=c(0,0,1,1,res)
}
decisionMatrix=matrix(decisionList,nrow=1,ncol=5)
#write.table(decisionMatrix,paste0("decisionMatrix_",i,"_",j,".csv"), sep = ",", row.names = FALSE, col.names=FALSE)
return (res);
}