# input are <infile> <ksize> [<infile> <ksize> ... ]
cat("-- reading arguments\n", sep = "")
Args <- commandArgs();

# R script for ranking #
dat=read.table(Args[3], header=TRUE, sep="\t", comment.char="")
dat$Coverage[dat$Coverage > 100] = 100

# ranking coverage, identity, bitscore, nmismatches, ngaps
ranks = apply(cbind(dat[,c(7:8, 13)], 1/dat[,10:11]),2, rank)
rownames(ranks) = dat$"Hit."

# summing ranks
ranks.s = apply(ranks, 1, sum)

# taking the one with highest ranks, if tie then chose by identity
m = max(ranks.s)
best = which(ranks.s == m)
if (length(best) != 1) {
   hs = names(best)
   best = hs[which.min(dat[dat$"Hit." %in% hs,"Identity"])]
   best_hit = best
} else {
   best_hit = names(best)
}

# check if best are equally ranked #

# write sorted output #
write.table(dat[order(ranks.s, decreasing=TRUE),], file=paste(c(Args[3], ".tab"), sep="", collapse=""), sep="\t", quote=FALSE, row.names=FALSE)

# write out
output = c(paste(c("Best hit is: ", best_hit), sep="", collapse=""),
           paste(c("Rank of hits (best-worst) (ties not properly sorted):", names(ranks.s[order(ranks.s, decreasing=TRUE)])), sep=" ", collapse=" "))
write(output, file=paste(c(Args[3], ".parse.txt"))
