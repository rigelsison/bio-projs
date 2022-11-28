library(tidyverse)

#read in raw CFX file
#substitute your file name in "csv"
#ex. rawdata <- read_csv("qpcr_draft.csv")

rawdata <- read_csv(csv)

clean <- rawdata[, c("Content", "Cq", "Cq Mean", "Cq Std. Dev", "Starting Quantity (SQ)", "SQ Mean", "SQ Std. Dev")]

clean <- clean[grepl("^Unkn", clean$Content),]

alone <- clean %>%
  group_by(Content) %>%
  filter(n()<2)

alone$Cq <- NA
alone$`Starting Quantity (SQ)` <- NA

clean <- rbind(clean, alone)
clean <- clean[order(clean$Content), ]

row_even <- seq_len(nrow(clean)) %% 2
data_row_even <- clean[row_even == 0, ]
data_row_odd <- clean[row_even == 1, ]

data_row_odd$Cq.2 <- data_row_even$Cq
data_row_odd$SQ.2 <- data_row_even$`Starting Quantity (SQ)`

data <- data_row_odd

#create final dataframe
final <- data.frame(matrix(ncol = 10, nrow = length(unique(data$Content))))
colnames(final) <- c("Sample", "Target Gene", "Fluorophore", "Rep1 - SQ", "Rep1 - Ct", "Rep2 - SQ", "Rep2 - Ct", "Average - SQ", "Average - Ct", "SQ Std Deviation")

final$`Sample` <- data$Content
final$`Target Gene` <- "16s"
final$`Fluorophore` <- "SYBR"
final$`Rep1 - SQ` <- data$`Starting Quantity (SQ)`
final$`Rep1 - Ct` <- data$Cq
final$`Rep2 - SQ` <- data$SQ.2
final$`Rep2 - Ct` <- data$Cq.2
final$`Average - SQ` <- data$`SQ Mean`
final$`Average - Ct` <- data$`Cq Mean`
final$`SQ Std Deviation` <- data$`SQ Std. Dev`

#output to csv; you can change "qpcr_draft" to whatever file name you want. This will save to your working directory.
write.csv(final, "qpcr_draft.csv", row.names = FALSE)
