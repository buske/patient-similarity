library(latticeExtra)

data <- read.table('match.out', header=FALSE)

colnames(data) <- c('a', 'b', 'Za', 'Zb', 'MI', 'ICa', 'ICb', 'Gn')

find.matching <- function(data, samples) {
  data$a %in% samples & data$b %in% samples
}

find.containing <- function(data, cohort) {
  term <- paste("^", cohort, "[_-]", sep="")
  grepl(term, data$a) & grepl(term, data$b)
}


data$class <- 0

data$class[find.matching(data,
c('MFDM_2005111',
'MFDM_2005112',
'MFDM_M43377',
'233_56920BD',
'1038_120760Q',
'MFDM_M46095'))] <- 1

for (cohort in c(113, 165, 166, 174, 181, 207, 242, 255, 258, 380, 393)) {
  data$class[find.containing(data, cohort)] <- 1
}

##xyplot(mi ~ max, data, groups=class, col=c("gray", "red"), xlab="Maximum patient IC", ylab="Mutual information", sub="TP: 26, FP: 3, FN: 10, TN: 3201")
##xyplot(score ~ min | a, subset(data, a %in% subset(data, class>0)$a), groups=class, xlab="Min information content", ylab="Score", pch=20, col=c("gray", "red"), par.strip.text=list(cex=0.6))