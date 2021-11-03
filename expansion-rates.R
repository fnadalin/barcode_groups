
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
        cat("\nUsage: <reference> <query> <out>\n")
        cat("\n<P0>   group.tsv file for baseline condition (P0)\n")
        cat("<T>    group.tsv file with condition at time T\n")
        cat("<out>  output file containing the expansion rate at time T for each GBC in P0\n\n")
        q()
}

P0_FILE <- args[1]
T_FILE <- args[2]
OUT_FILE <- args[3]

df_p0 <- read.table(P0_FILE, header = TRUE, sep = "\t")
df_t <- read.table(T_FILE, header = TRUE, sep = "\t")

tot_p0 <- sum(df_p0$group_count)
tot_t <- sum(df_t$group_count) 

count_t <- rep(0,nrow(df_p0))
comm <- df_p0$hub[df_p0$hub %in% df_t$hub]
count_t[match(comm,df_p0$hub)] <- df_t$group_count[df_t$hub %in% comm]

cpm_p0 <- 1e6*df_p0$group_count/tot_p0
cpm_t <- 1e6*count_t/tot_t

exp_rate <- log2(cpm_t+1) - log2(cpm_p0+1)
df_exp_rate <- data.frame(GBC = df_p0$hub, exp.rate = exp_rate)
df_exp_rate <- df_exp_rate[order(df_exp_rate$exp.rate, decreasing = TRUE),]

write.table(df_exp_rate, file = OUT_FILE, sep = "\t", row.names = FALSE, quote = FALSE)


q()

