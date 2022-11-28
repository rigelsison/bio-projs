library(tidyverse)

#read in raw CFX file
#substitute your file name in "csv"
#ex. rawdata <- read_csv("fam_hex.csv")

rawdata <- read_csv("BGI-LOD-Trial2.1.csv") 

#keep relevant columns
clean <- rawdata[, c("Fluor", "Content", "Cq", "Cq Mean", "Cq Std. Dev")]

#remove standards
clean <- clean[grepl("^Unkn", clean$Content),]

#find samples with no replicates and add dummy duplicates
alone <- clean %>%
  group_by(Fluor, Content) %>%
  filter(n() <2)

alone$Cq <- NA
alone$`Starting Quantity (SQ)` <- NA

#merge together
clean <- rbind(clean, alone)
clean <- clean[order(clean$Fluor, clean$Content), ]

#split between replicates
row_even <- seq_len(nrow(clean)) %% 2
data_row_even <- clean[row_even == 0, ]
data_row_odd <- clean[row_even == 1, ]

data_row_odd$Cq.2 <- data_row_even$Cq

#replicates will now be in same row
data <- data_row_odd

#split dataframe by fluorophore
data_split <- split(data, data$Fluor) 

fam <- data.frame(data_split[1])
hex <- data.frame(data_split[2])

#group by sample
merged <- merge(fam, hex, by.x = 'FAM.Content', by.y = 'HEX.Content')

#remove unnecessary columns
drop <- subset(merged, select = -c(FAM.Fluor, HEX.Fluor))

#create empty final dataframe that matches Benchling draft result table
final <- data.frame(matrix(ncol = 9, nrow = length(unique(data$Content))))
colnames(final) <- c('Sample Description', 'FAM Ct - R1', 'FAM Ct - R2', 'FAM Ct Mean', 'FAM Ct Std. Dev', 'HEX Ct - R1', 'HEX Ct - R2', 'HEX Ct Mean', 'HEX Ct Std. Dev')
final$`Sample Description` <- c(unique(data$Content))

#fill in final dataframe
final$`Sample Description` <- drop$FAM.Content
final$`FAM Ct - R1` <- drop$FAM.Cq
final$`FAM Ct - R2` <- drop$FAM.Cq.2
final$`FAM Ct Mean` <- drop$FAM.Cq.Mean
final$`FAM Ct Std. Dev` <- drop$FAM.Cq.Std..Dev
final$`HEX Ct - R1` <- drop$HEX.Cq
final$`HEX Ct - R2` <- drop$HEX.Cq.2
final$`HEX Ct Mean` <- drop$HEX.Cq.Mean
final$`HEX Ct Std. Dev` <- drop$HEX.Cq.Std..Dev

#output to csv; you can change "promega" to whatever file name you want. This will save to your working directory.
write.csv(final, "fam_hex.csv", row.names = FALSE)

