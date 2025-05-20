# 读取包含读长信息的文件RdRp_DB_Read_length.txt，以制表符为分隔符，不将字符型数据转换为因子
ReadLen = read.table("virushostdb_cds_list.txt", header = F, sep = "\t", stringsAsFactors = F)

# 设置输入目录为当前目录
input = "."

# 获取目录中的文件列表
samples = list.files(input)

# 遍历每个文件
for (s in samples) {
	# 构建文件路径和文件名
	file1 = paste(input, "/", s, "/", s, "_blastout.out", sep = "")
	file2 = paste(input, "/", s, "/", s, "_Blastout.out", sep = "")
	
	# 如果文件1和文件2都不存在，则跳过当前文件
	if (!(file.exists(file1) | file.exists(file2))) {
		next
	} else {
		print(s)
		
		# 如果存在文件1，则读取文件1；否则读取文件2
		if (file.exists(file1)) {
			List = read.table(file1, header = F, sep = "\t", stringsAsFactors = F)
		} else {
			List = read.table(file2, header = F, sep = "\t", stringsAsFactors = F)
		}
		
		# 从List中筛选出第11列小于1e-5的行
		List = List[List[, 11] < 1e-5,]
		
		# 根据List中的第2列查找对应的读长信息，并处理读长不存在的情况
		Len = lapply(List[, 2], function(x) {
			y = ReadLen[(ReadLen[, 1] %in% c(x, tolower(x), toupper(x))), 2]
			z = ifelse(length(y) == 0, "NA", y)
		})
		Len = unlist(Len)
		
		# 计算List中的第4列除以读长的比例
		Percentage = List[, 4] / Len
		
		# 构建临时数据框Temp，包括List中的第1列、第4列、比例、第11列、和第7、8列
		Temp = data.frame(List[, 1], List[, 4], Percentage, List[, 11], List[, 7:8])
		
		# 从Temp中筛选出比例大于0.75的行
		Temp = Temp[Temp[, 3] > 0.75,]
		
		# 按照第2列和第3列进行降序排列
		Temp = Temp[rev(order(Temp[, 2], Temp[, 3])),]
		
		# 获取Temp中第1列的唯一值
		genes1 = unique(Temp[, 1])
		
		# 遍历每个基因，获取Temp中对应基因的第1列、第5列和第6列的行
		genes2 = lapply(genes1, function(x) {
			Temp[Temp[, 1] == x, ][1, c(1, 5, 6)]
		})
		
		# 将genes2中的结果按行连接起来
		Result = do.call(rbind, genes2)
		
		# 将Result写入到文件s_RdRp_grep_list.bed中，以制表符为分隔符
		write.table(Result, paste(input, "/", s, "/", s, "_RdRp_grep_list.bed", sep = ""), quote = F, sep = "\t", col.names = F, row.names = F)
	}
}
