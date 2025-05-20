ReadLen=read.table("RdRp_DB_Read_length.txt",header=F,sep="\t",stringsAsFactors=F)
input="."
samples=list.files(input)
for(s in samples){
	file1=paste(input,"/",s,"/",s,"_rdrp_list.txt",sep="")
	file2=paste(input,"/",s,"/",s,"_RdRp_list.txt",sep="")
	if(!(file.exists(file1)|file.exists(file2))){
		next
	}else{
		print(s)
		if(file.exists(file1)){
			List=read.table(file1,header=F,sep="\t",stringsAsFactors=F)
		}else{
			List=read.table(file2,header=F,sep="\t",stringsAsFactors=F)
		}
		List=List[List[,11]<1e-5,]	
		Len=lapply(List[,2],function(x){y=ReadLen[(ReadLen[,1] %in% c(x,tolower(x),toupper(x))),2];z=ifelse(length(y)==0,"NA",y)})
		Len=unlist(Len)
		Percentage=List[,4]/Len
		Temp=data.frame(List[,1],List[,4],Percentage,List[,11],List[,7:8])
		Temp=Temp[Temp[,3]>0.75,]
		Temp=Temp[rev(order(Temp[,2],Temp[,3])),]
		genes1=unique(Temp[,1])
		genes2=lapply(genes1,function(x){Temp[Temp[,1]==x,][1,c(1,5,6)]})
		Result=do.call(rbind,genes2)
		write.table(Result,paste(input,"/",s,"/",s,"_RdRp_grep_list.bed",sep=""),quote=F,sep="\t",col.names=F,row.names=F)
	}
}			
