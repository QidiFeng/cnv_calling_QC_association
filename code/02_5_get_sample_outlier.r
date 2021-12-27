args <- commandArgs(TRUE)
infile <- args[1]
cohortfile <- args[2]
outfile <- args[3]
outfile2 <- args[4]

data_total <- read.table(infile, header=T)
#outlier0 <- data[rowSums(is.na(data))>0,]
#outlier0$OUTLIER <- "value_missing"
##remove those rows with NA.
#data <- data[complete.cases(data), ]
cohortinfo <- read.table(cohortfile)
cohorts <- cohortinfo$V1

lrr_sd <- c()
baf_sd <- c()
gcwf1 <- c()
gcwf2 <- c()
cnv_num <- c()
cnv_prop <- c()
sample_size <- c()
outlier_size <- c()
data_list <- list()
outlier_list <- list()
for (i in 1:length(cohorts)){
	cohort <- cohorts[i]
	print(cohort)
	data <- data_total[which(data_total$Cohort == cohort),]
	data$proportion <- data$CNV_total_length/3.1e9
	sample_size <- c(sample_size,nrow(data))
	criteria1 <- median(data$LRR_SD) + 3*sd(data$LRR_SD)
	outlier_list[[1]] <- data[which(data$LRR_SD > criteria1),]
	if(nrow(outlier_list[[1]])>0){
	outlier_list[[1]]$OUTLIER <- "LRR_SD"
	}
	criteria2 <- median(data$BAF_SD) + 3*sd(data$BAF_SD)
	outlier_list[[2]] <- data[which(data$BAF_SD > criteria2),]
	if(nrow(outlier_list[[2]])>0){
	outlier_list[[2]]$OUTLIER <- "BAF_SD"
	}
	criteria3 <- median(data$GCWF) + 3*sd(data$GCWF)
	criteria4 <- median(data$GCWF) - 3*sd(data$GCWF)
	outlier_list[[3]] <- data[which(data$GCWF > criteria3),]
	outlier_list[[4]] <- data[which(data$GCWF < criteria4),] 
	if(nrow(outlier_list[[3]])>0){
	outlier_list[[3]]$OUTLIER <- "GCWF"
	}
	if(nrow(outlier_list[[4]])>0){
	outlier_list[[4]]$OUTLIER <- "GCWF"
	}

	lrr_sd <- c(lrr_sd,criteria1)
	baf_sd <- c(baf_sd,criteria2)
	gcwf1 <- c(gcwf1,criteria3)
	gcwf2 <- c(gcwf2,criteria4)

	criteria5 <- median(data$CNV_num) + 3*sd(data$CNV_num)
	print(median(data$CNV_num))
	print(sd(data$CNV_num))
	outlier_list[[5]] <- data[which(data$CNV_num > criteria5),]
	if(nrow(outlier_list[[5]])>0){
	outlier_list[[5]]$OUTLIER <- "CNV_num"
	}
	criteria6 <- median(data$proportion) + 3*sd(data$proportion)
	outlier_list[[6]] <- data[which(data$proportion > criteria6),]
	if(nrow(outlier_list[[6]])>0){
	outlier_list[[6]]$OUTLIER <- "CNV_proportion"
	}

	cnv_num <- c(cnv_num,criteria5)
	cnv_prop <- c(cnv_prop,criteria6)
	outlier <- do.call(rbind,outlier_list)
	outlier_size <- c(outlier_size,length(unique(outlier$Sample)))
	data_list[[i]] <- outlier
}
	criteria <- data.frame(BAF_SD=baf_sd, LRR_SD=lrr_sd,GCWF1=gcwf1,GCWF2=gcwf2,CNV_NUM=cnv_num,CNV_PROP=cnv_prop,SAMPLE_SIZE=sample_size,OUTLIER=outlier_size)
	criteria <- round(criteria,4)
	criteria$COHORT <- cohorts
	res <- do.call(rbind,data_list)
	write.table(res, file=outfile, quote=F, row.names=F,sep="\t")
	write.table(criteria, file=outfile2, quote=F, row.names=F,sep="\t")
