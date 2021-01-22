library(data.table)
#library(funr)
library(xgboost)
#script.path=sys.script()
#script.path=strsplit(script.path,"\\/")[[1]]
#script.path=paste(script.path[-length(script.path)],collapse="/")
data.path=paste(script.path,"/data/",sep="")
bedtools=paste(script.path,"/tools/bedtools",sep="")
args=commandArgs(T)
if(length(args)==0)
{
message("\tPlease specify a BED file\n\tUsage:\tXCNV CNV.bed")
quit(save="no")
}else if(args=="--help" | args=="-h")
{
message("\tUsage:\tXCNV CNV.bed")
quit(save="no")
}else if(!file.exists(args[1]))
{
message(paste("\tNo such file:",args[1]))
message("\tUsage:\tXCNV CNV.bed")
quit(save="no")
}


infile=args[1]
indata=fread(infile,header=F)[,1:4]
n.cols=ncol(indata)
if(grepl("chr",indata$V1[1]))
{
indata=fread(infile,header=F)
indata$V1=gsub("chr","",indata$V1)
}
n.cnvs=nrow(indata)
if(n.cnvs==1)
{
indata=rbind(indata,indata)
indata$V3[1]=indata$V2[1]+49
}
indata$V2=as.integer(indata$V2)
indata$V3=as.integer(indata$V3)
check.data=function(mat)
{
test.chr=sum(is.na(match(mat$V1,c(1:22,"X","Y"))))
test.pos=sum(is.na(mat$V2))+sum(is.na(mat$V3))
test.type=sum(is.na(match(mat$V4,c("gain","loss"))))
test.width=sum(mat$V2>mat$V3)
if(test.chr>0 | test.pos>0 | test.type>0 | test.width)
{
message("Your file format is invalid!")
quit(save="no")
}
}

tmp=check.data(indata)

idx=match(indata$V1,c(1:22,"X","Y"))
indata=indata[order(idx,indata$V2,indata$V3),]
infile.sorted=gsub(".bed$",".sort.bed",infile)
fwrite(indata,infile.sorted,sep='\t',row.names=F,col.names=F,quote=F)
labels=paste(indata$V1,indata$V2,indata$V3,indata$V4,sep="_")
cnv.length=indata$V3-indata$V2
genome.file=paste(data.path,"genome.txt",sep="")

#########################
#### ljb26 start
anno.file1=paste(data.path,"hg19_ljb26_all_converted.vcf",sep="")
outfile1=gsub(".bed$",".tmp1.bed",infile)
cmd1=paste(bedtools,"intersect -wao -sorted -g",genome.file ,"-a",infile.sorted,"-b",anno.file1,">",outfile1)
system(cmd1)
anno.data1=fread(outfile1,header=F)
anno.file1.header=fread(anno.file1,nrow=1,header=T)
colnames(anno.data1)[(n.cols+1):(n.cols+ncol(anno.file1.header))]=colnames(anno.file1.header)
labels1=factor(paste(anno.data1$V1,anno.data1$V2,anno.data1$V3,anno.data1$V4,sep="_"),levels=labels)
scores1=sapply(colnames(anno.file1.header)[4:12],function(x){
tapply(anno.data1[[x]],labels1,function(t){
sum(t>0)/length(t)
})
})
missing.values=c(0,0,-12.30,-11.958,-20.000,0.0003)
names(missing.values)=colnames(anno.file1.header)[13:18]
scores2=sapply(colnames(anno.file1.header)[13:18],function(x){
tapply(anno.data1[[x]],labels1,function(t){
t[t=="."]=missing.values[x]
t=as.numeric(t)
mean(t)
})
})
scores=cbind(scores1,scores2)
#### ljb26 end
#########################

#########################
#### CDTS start
anno.file2=paste(data.path,"CDTS_percentile.txt",sep="")
outfile2=gsub(".bed$",".tmp2.bed",infile)
cmd2=paste(bedtools,"intersect -wao -sorted -g", genome.file ,"-a",infile.sorted,"-b",anno.file2,">",outfile2)
system(cmd2)
anno.data2=fread(outfile2,header=F)
CDTS=anno.data2[[n.cols+4]]
CDTS[CDTS=="."]=100
CDTS=as.numeric(CDTS)
labels2=factor(paste(anno.data2$V1,anno.data2$V2,anno.data2$V3,anno.data2$V4,sep="_"),levels=labels)
CDTS_1st=tapply(CDTS,labels2,function(x){sum(x<1)*10})/cnv.length
CDTS_5th=tapply(CDTS,labels2,function(x){sum(x<5)*10})/cnv.length
#### CDTS end
#########################

#########################
#### frequency start
anno.file3=paste(data.path,"merged.cnv.data.bed",sep="")
outfile3=gsub(".bed$",".tmp3.bed",infile)
cmd3=paste(bedtools, "intersect -wao -f 0.5 -F 0.5 -sorted -g", genome.file ,"-a",infile.sorted,"-b",anno.file3,">",outfile3)
system(cmd3)
anno.data3=fread(outfile3,header=F)
labels3=factor(paste(anno.data3$V1,anno.data3$V2,anno.data3$V3,anno.data3$V4,sep="_"),levels=labels)
n.samples=87935
gain.samples=tapply(anno.data3$V9[anno.data3$V8=="gain"],labels3[anno.data3$V8=="gain"],function(x){x=unique(unlist(strsplit(x,",")));x=x[x!="."];y=length(x)})
loss.samples=tapply(anno.data3$V9[anno.data3$V8=="loss"],labels3[anno.data3$V8=="loss"],function(x){x=unique(unlist(strsplit(x,",")));x=x[x!="."];y=length(x)})
gain.samples[is.na(gain.samples)]=0
loss.samples[is.na(loss.samples)]=0
gain.freq=gain.samples/n.samples
loss.freq=loss.samples/n.samples

#### frequency end
#########################

#########################
#### ENCODE start
anno.file4=paste(data.path,"hg19-ccREs.bed",sep="")
outfile4=gsub(".bed$",".tmp4.bed",infile)
cmd4=paste(bedtools,"intersect -wao -sorted -g", genome.file, "-a",infile.sorted,"-b",anno.file4,">",outfile4)
system(cmd4)
anno.data4=fread(outfile4,header=F)
#labels4=factor(paste(anno.data4$V1,anno.data4$V2,anno.data4$V3,anno.data4$V4,sep="_"),levels=labels)
encode=function(overlap,n.cols,labels)
{
cols=c("pELS","CTCF-bound","PLS","dELS", "CTCF-only","DNase-H3K4me3")
tmp=t(sapply(strsplit(as.matrix(overlap[[n.cols+4]]),","),function(x){y=rep(0,6);x=x[x!="."];if(length(x)>0){y[match(x,cols)]=1};y}))
labels=factor(paste(overlap$V1,overlap$V2,overlap$V3,overlap$V4,sep="_"),levels=labels)
len=as.matrix(overlap[[n.cols+5]])
len[len=="."]=0
tmp1=apply(tmp,2,function(x){x*len})
m1=tapply(1:length(labels),factor(labels,levels=unique(labels)),function(idx){y=apply(matrix(tmp1[idx,],ncol=6),2,sum)})
m1=eval(as.call(c(rbind,m1)))
len1=sapply(strsplit(rownames(m1),"_"),function(x){as.numeric(x[3])-as.numeric(x[2])})
m2=apply(m1,2,function(x){x/len1})
colnames(m2)=cols
m2
}
encode.data=encode(anno.data4,n.cols,labels)

###
# haploinsufficiency
anno.file5=paste(data.path,"gencode_v19_features.bed",sep="")
outfile5=gsub(".bed$",".tmp5.bed",infile)
cmd5=paste(bedtools,"intersect -wao -sorted -g", genome.file, "-a",infile.sorted,"-b",anno.file5,">",outfile5)
system(cmd5)
anno.data5=fread(outfile5,header=F)
#labels5=factor(paste(anno.data5$V1,anno.data5$V2,anno.data5$V3,anno.data5$V4,sep="_"),levels=labels)
hap.cal=function(data,n.cols,labels)
{
cols=c("pLI"  ,"Episcore"      ,"GHIS")
labels=factor(paste(data$V1,data$V2,data$V3,data$V4,sep="_"),levels=labels)
m=tapply(1:nrow(data),factor(labels,levels=unique(labels)),function(idx){tmp=as.matrix(data[idx,(n.cols+5):(n.cols+7)]);tmp[tmp=="."]=0;apply(tmp,2,function(t){max(as.numeric(t))})})
m=eval(as.call(c(rbind,m)))
colnames(m)=cols
m
}
hap.data=hap.cal(anno.data5,n.cols,labels)

#####
all.vars=cbind(indata,scores,CDTS_1st,CDTS_5th,gain.freq,loss.freq,encode.data,hap.data)
output.file=gsub(".bed$",".output.csv",infile)
colnames(all.vars)[1]="#V1"
if(n.cnvs==1)
{
all.vars=all.vars[-1,]
}
#fwrite(all.vars,output.file,sep='\t',row.names=F,quote=F)
tmp=file.remove(c(outfile1,infile.sorted,outfile2,outfile3, outfile4, outfile5))
all.vars$Type=ifelse(all.vars$V4=="gain",1,0)
all.vars$Length=all.vars$V3-all.vars$V2+1
all.vars1=data.matrix(all.vars[,-c(1:4)])
load(paste(data.path,"xcnv.model.Rdata",sep=""))
xcnv.score=predict(xcnv.model,newdata=all.vars1)
all.vars=cbind(all.vars,MVP_score=xcnv.score)
colnames(all.vars)[1:4]=c("Chr","Start","End","Type")
fwrite(all.vars,output.file)
