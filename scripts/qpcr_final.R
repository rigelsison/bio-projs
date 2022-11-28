library(tidyverse)

#read in raw CFX file
#substitute your file name in "csv"
#ex. rawdata <- read_csv("qpcr_final.csv")

rawdata <- read_csv(csv)

clean <- rawdata[, c("Content", "Cq", "Cq Mean", "Cq Std. Dev", "Starting Quantity (SQ)", "SQ Mean", "SQ Std. Dev")]

#filter for only unknowns
clean <- clean[grepl("^Unkn", clean$Content),]

#find samples with no replicates and add dummy duplicates
alone <- clean %>%
  group_by(Content) %>%
  filter(n() <2)

alone$Cq <- NA
alone$`Starting Quantity (SQ)` <- NA

#merge together
clean <- rbind(clean, alone)
clean <- clean[order(clean$Content), ]

#split between replicates
row_even <- seq_len(nrow(clean)) %% 2
data_row_even <- clean[row_even == 0, ]
data_row_odd <- clean[row_even == 1, ]

data_row_odd$Cq.2 <- data_row_even$Cq
data_row_odd$SQ.2 <- data_row_even$`Starting Quantity (SQ)`

#replicates will now be in same row
data <- data_row_odd

#create final dataframe
final <- data.frame(matrix(ncol = 7, nrow = length(unique(data$Content))))
colnames(final) <- c("Sample ID", "Target Gene", "Fluorophore", "Average - SQ", "Average - Ct", "SQ Std Deviation", "Ct Std Deviation" )

final$`Sample ID` <- data$Content
final$`Target Gene` <- "16s"
final$`Fluorophore` <- "SYBR"
final$`Average - SQ` <- data$`SQ Mean`
final$`Average - Ct` <- data$`Cq Mean`
final$`SQ Std Deviation` <- data$`SQ Std. Dev`
final$`Ct Std Deviation` <- data$`Cq Std. Dev`

#output to csv; you can change "qpcr_final" to whatever file name you want. This will save to your working directory.
write.csv(final, "qpcr_final.csv", row.names = FALSE)
