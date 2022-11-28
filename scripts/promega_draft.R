library(tidyverse)

#read in raw CFX file
#substitute your file name in "csv"
#ex. rawdata <- read_csv("promega.csv")

rawdata <- read_csv(csv) 

#keep relevant columns
clean <- rawdata[, c("Fluor", "Content", "Cq", "Starting Quantity (SQ)", "SQ Mean", "SQ Std. Dev")]

#remove standards
clean <- clean[grepl("^Unkn", clean$Content),]

#remove Quasar 705
clean <- filter(clean, Fluor != 'Quasar 705')

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
data_row_odd$SQ.2 <- data_row_even$`Starting Quantity (SQ)`

#replicates will now be in same row
data <- data_row_odd

#split dataframe by fluorophore
data_split <- split(data, data$Fluor) 

cal_gold <- data.frame(data_split[1])
fam <- data.frame(data_split[2])
quasar_670 <- data.frame(data_split[3])

#group by sample
merged <- merge(cal_gold, fam, by.x = 'Cal.Gold.540.Content', by.y = 'FAM.Content')
merged <- merge(merged, quasar_670, by.x = 'Cal.Gold.540.Content', by.y = 'Quasar.670.Content')

#remove unnecessary columns
drop <- subset(merged, select = -c(Cal.Gold.540.Fluor, FAM.Fluor, Quasar.670.Fluor))

#create empty final dataframe that matches Benchling draft result table
final <- data.frame(matrix(ncol = 19, nrow = length(unique(data$Content))))
colnames(final) <- c('Sample Description', '300 bp Ct - R1', '300 bp Ct - R2', '300 bp SQ - R1', '300 bp SQ - R2', '300 bp SQ - Average', '300 bp SQ - Std. Dev', '150 bp Ct - R1', '150 bp Ct - R2', '150 bp SQ - R1', '150 bp SQ - R2', '150 bp SQ - Average', '150 bp SQ - Std. Dev', '75 bp Ct - R1', '75 bp Ct - R2', '75 bp SQ - R1', '75 bp SQ - R2', '75 bp SQ - Average', '75 bp SQ - Std. Dev')
final$`Sample Description` <- c(unique(data$Content))

#fill in final dataframe
final$`Sample Description` <- drop$Cal.Gold.540.Content
final$`300 bp Ct - R1` <- drop$Cal.Gold.540.Cq
final$`300 bp Ct - R2` <- drop$Cal.Gold.540.Cq.2
final$`300 bp SQ - R1` <- drop$Cal.Gold.540.Starting.Quantity..SQ.
final$`300 bp SQ - R2` <- drop$Cal.Gold.540.SQ.2
final$`300 bp SQ - Average` <- drop$Cal.Gold.540.SQ.Mean
final$`300 bp SQ - Std. Dev` <- drop$Cal.Gold.540.SQ.Std..Dev
final$`150 bp Ct - R1` <- drop$FAM.Cq
final$`150 bp Ct - R2` <- drop$FAM.Cq.2
final$`150 bp SQ - R1` <- drop$FAM.Starting.Quantity..SQ.
final$`150 bp SQ - R2` <- drop$FAM.SQ.2
final$`150 bp SQ - Average` <- drop$FAM.SQ.Mean
final$`150 bp SQ - Std. Dev` <- drop$FAM.SQ.Std..Dev
final$`75 bp Ct - R1` <- drop$Quasar.670.Cq
final$`75 bp Ct - R2` <- drop$Quasar.670.Cq.2
final$`75 bp SQ - R1` <- drop$Quasar.670.Starting.Quantity..SQ.
final$`75 bp SQ - R2` <- drop$Quasar.670.SQ.2
final$`75 bp SQ - Average` <- drop$Quasar.670.SQ.Mean
final$`75 bp SQ - Std. Dev` <- drop$Quasar.670.SQ.Std..Dev

#output to csv; you can change "promega" to whatever file name you want. This will save to your working directory.
write.csv(final, "promega_draft.csv", row.names = FALSE)

