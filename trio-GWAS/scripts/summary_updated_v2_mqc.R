#!/usr/bin/env Rscript

errorlist <- list()
warninglist <- list()

library(data.table)
suppressMessages(library(matrixStats))

args <- (commandArgs(TRUE));
phenotype_file <- as.character(args[1]);
cov_file <- as.character(args[2]);
phen_list <- as.character(args[3])
cov_list <- as.character(args[4])
out_file <- as.character(args[5])

ph <- fread(phenotype_file, h=TRUE)
cov <- fread(cov_file, h=TRUE)

## 1) Read lists as character vectors (not data.tables)
plist <- fread(phen_list, h=FALSE)[[1]]
plist <- plist[plist != ""]
clist <- fread(cov_list,  h=FALSE)[[1]]
clist <- clist[clist != ""]

output <- NULL

for (i in seq_along(plist)) {
  ## 2) Take the i-th phenotype name as a string
  phen <- plist[i]
  
  ## 3) Don't shadow the merge() function
  mrg <- merge(ph, cov, by = "IID")
  
  yob_column <- names(cov)[grepl("^Age_", names(cov), ignore.case = TRUE)][1]
  sex_column <- names(cov)[grepl("^sex$", names(cov), ignore.case = TRUE)][1]
  
  ph2 <- subset(mrg, select=c("IID", phen, sex_column, yob_column,
                              "PC1","PC2","PC3","PC4","PC5",
                              "PC6","PC7","PC8","PC9","PC10",
                              "PC11","PC12","PC13","PC14","PC15",
                              "PC16","PC17","PC18","PC19","PC20"))
  names(ph2)[names(ph2) == phen] <- "Outcome"
  names(ph2) <- c("IID","Outcome","Sex","Age",
                  "PC1","PC2","PC3","PC4","PC5",
                  "PC6","PC7","PC8","PC9","PC10",
                  "PC11","PC12","PC13","PC14","PC15",
                  "PC16","PC17","PC18","PC19","PC20")
  
  mean <- mean(ph2$Outcome, na.rm=TRUE)
  sd <- sd(ph2$Outcome, na.rm=TRUE)
  median <- median(ph2$Outcome, na.rm=TRUE)
  
  q1 <- as.numeric(summary(ph2$Outcome)[2])
  q3 <- as.numeric(summary(ph2$Outcome)[5])
  min <- as.numeric(summary(ph2$Outcome)[1])
  max <- as.numeric(summary(ph2$Outcome)[6])
  
  miss <- is.na(ph2$Outcome)
  N <- sum(!is.na(ph2$Outcome))        ## (optional but clearer than summary(miss)[2])
  
  model1 <- lm(Outcome ~ Sex + Age + PC1 + PC2 + PC3 + PC4 + PC5 +
                 PC6 + PC7 + PC8 + PC9 + PC10 +
                 PC11 + PC12 + PC13 + PC14 + PC15 +
                 PC16 + PC17 + PC18 + PC19 + PC20, data = ph2)
  sd_resid <- sd(resid(model1))
  
  stats <- data.frame(phen, N, mean, sd, sd_resid, median, min, max, q1, q3)
  output <- rbind(output, stats)
}

for (i in seq_along(clist)) {
  phen <- clist[i]
  cov2 <- subset(cov, select = c("IID", phen))
  names(cov2) <- c("IID", "Covariate")
  
  mean <- mean(cov2$Covariate, na.rm = TRUE)
  sd <- sd(cov2$Covariate, na.rm = TRUE)
  median <- median(cov2$Covariate, na.rm = TRUE)
  
  q1 <- as.numeric(summary(cov2$Covariate)[2])
  q3 <- as.numeric(summary(cov2$Covariate)[5])
  min <- as.numeric(summary(cov2$Covariate)[1])
  max <- as.numeric(summary(cov2$Covariate)[6])
  
  N <- sum(!is.na(cov2$Covariate))
  sd_resid <- NA
  
  stats <- data.frame(phen, N, mean, sd, sd_resid, median, min, max, q1, q3)
  output <- rbind(output, stats)
}

## 4) Fix write.table args
write.table(output, file = out_file, row.names = FALSE, quote = FALSE, sep = "\t")
