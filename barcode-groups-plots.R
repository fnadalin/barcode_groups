# Plot some stats on the barcodes

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
	cat("\nUsage: Rscript barcode-groups-plots.R <indir> <outdir> <title> <min_group_count>\n\n")
	cat("<indir>     where the results of ./barcode-groups are saved)\n")
	cat("<outdir>    where plots are saved\n")
	cat("<title>     to be put on top of each plot\n\n")
	q()
}

INDIR <- args[1]
OUTDIR <- args[2]
TITLE <- args[3]
MIN_GROUP_COUNT <- as.numeric(args[4])

library("scales")

PARAMS <- file.path(INDIR, "params.txt")

NODE_ATTR <- file.path(INDIR, "node_attr.tsv")
EDGE_ATTR <- file.path(INDIR, "edge_attr.tsv")  
NEIGH_ATTR <- file.path(INDIR, "neighbor_attr.tsv")  
GROUP_SIZE <- file.path(INDIR, "group_size_freq.tsv") 
GROUPS <- file.path(INDIR, "groups.tsv")
GROUPS_ALL <- file.path(INDIR, "groups_all.tsv")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# OUT_NODE_ATTR <- file.path(OUTDIR, "node_attr.pdf")
# OUT_EDGE_ATTR <- file.path(OUTDIR, "edge_attr.pdf")  
# OUT_NEIGH_ATTR <- file.path(OUTDIR, "neighbor_attr.pdf")  
# OUT_EDGE_ATTR_ZOOM <- file.path(OUTDIR, "edge_attr_zoom.pdf")  
# OUT_GROUP_SIZE <- file.path(OUTDIR, "group_size_freq.pdf") 

OUT_NODE_ATTR <- file.path(OUTDIR, "node_attr.pdf")
OUT_EDGE_ATTR <- file.path(OUTDIR, "edge_attr.pdf")
OUT_NEIGH_ATTR <- file.path(OUTDIR, "neighbor_attr.pdf")
OUT_EDGE_ATTR_ZOOM <- file.path(OUTDIR, "edge_attr_zoom.pdf")
OUT_GROUP_SIZE <- file.path(OUTDIR, "group_size_freq.pdf")
OUT_HUB_IN_OUT_DEG <- file.path(OUTDIR, "hub_in_out_deg.pdf")
OUT_GROUP_IN_OUT_DEG <- file.path(OUTDIR, "group_in_out_deg.pdf")
OUT_HUB_VS_GROUP_COUNT <- file.path(OUTDIR, "hub_vs_group_count.pdf")
OUT_GROUP_CUM_FREQ <- file.path(OUTDIR, "cum_freq.pdf")

params <- read.table(PARAMS, sep = " ")
HAMM_DIST <- params[grep("mismatches", params[,1]), 2]

# neighbors / count plot for each barcode 

data <- read.table(NODE_ATTR, header = TRUE)
main <- "Barcode attributes"
sub <- paste("[ hamming distance =", HAMM_DIST, "]")
pdf(OUT_NODE_ATTR, width = 5, height = 5)
# png(OUT_NODE_ATTR, width = 280, height = 300)
plot(data$num_neighbors, data$counts, pch = 16, col = alpha("black", alpha = 0.1), xlab = "# neighbors", ylab = "read count", main = TITLE, sub = sub, log = "y")
dev.off()

# read count of barcode neighbors

data <- read.table(EDGE_ATTR, header = TRUE)
main <- "Barcode neighbors"
sub <- paste("[ hamming distance =", HAMM_DIST, "]")
pdf(OUT_EDGE_ATTR, width = 5, height = 5)
# png(OUT_EDGE_ATTR, width = 280, height = 300)
plot(data$count1, data$count2, pch = 16, col = alpha("black", alpha = 0.1), xlab = "barcode 1 read count", ylab = "barcode 2 read count", main = TITLE, sub = sub)
dev.off()

# read count of barcode neighbors (zoom version)

pdf(OUT_EDGE_ATTR_ZOOM, width = 5, height = 5)
# png(OUT_EDGE_ATTR_ZOOM, width = 280, height = 300)
plot(data$count1, data$count2, pch = 16, col = alpha("black", alpha = 0.1), xlab = "barcode 1 read count", ylab = "barcode 2 read count", main = TITLE, sub = sub, xlim = c(min(data$count2), max(data$count2)))
dev.off()

# read count of top 2 barcode neighbors

data <- read.table(NEIGH_ATTR, header = TRUE)
main <- "Top 2 non-mutual barcode neighbors"
sub <- paste("[ hamming distance =", HAMM_DIST, "]")
pdf(OUT_NEIGH_ATTR, width = 5, height = 5)
# png(OUT_NEIGH_ATTR, width = 280, height = 300)
plot(data$top_neigh_count, data$top2_neigh_count, pch = 16, col = alpha("black", alpha = 0.1), xlab = "top neighbor read count", ylab = "2nd top neighbor read count", main = TITLE, sub = sub)
dev.off()

# histogram of group size frequency

data <- read.table(GROUP_SIZE, header = TRUE)
main <- "Barcode groups"
sub <- paste("[ hamming distance =", HAMM_DIST, "]")
pdf(OUT_GROUP_SIZE, width = 6, height = 5)
# png(OUT_GROUP_SIZE, width = 300, height = 300)
barplot(height = data$num_groups, names.arg = data$size, main = TITLE, sub = sub, xlab = "group size", ylab = "frequency")
box()
dev.off()

# line plot with in-out degree

data <- read.table(GROUPS_ALL, header = TRUE)
main <- "In/out hub degree"
sub <- paste("[ hamming distance =", HAMM_DIST, "]")
pdf(OUT_GROUP_IN_OUT_DEG, width = 5, height = 5)
# png(OUT_GROUP_IN_OUT_DEG, width = 280, height = 300)
plot(data$group_count, data$out_hub_degree/(data$size-1+data$out_hub_degree), ylim = c(0,1), log = "x", type = "h", main = TITLE, sub = sub, xlab = "group read count", ylab = "out degree / hub degree")
abline(v = MIN_GROUP_COUNT, col = "red")
dev.off()

# dot plot with hub read count vs group read count

main <- "In/out hub degree"
sub <- paste("[ hamming distance =", HAMM_DIST, "]")
pdf(OUT_HUB_VS_GROUP_COUNT, width = 5, height = 5)
# png(OUT_HUB_VS_GROUP_COUNT, width = 280, height = 300)
plot(data$hub_count, (data$group_count - data$hub_count) / data$hub_count, ylim = c(0,1), main = TITLE, sub = sub, xlab = "hub read count", ylab = "% read count increase in group")
dev.off()

# dot plot with cumulative GBC groups frequencies

data <- read.table(GROUPS, header = TRUE)
n <- nrow(data)
num_reads <- sum(data$group_count)
cum_gbc_count <- rep(0, n+1)
for (i in 1:n) {
    cum_gbc_count[i+1] <- cum_gbc_count[i] + data$group_count[i]
}
cum_gbc_count <- cum_gbc_count[2:length(cum_gbc_count)]
pdf(OUT_GROUP_CUM_FREQ, width = 5, height = 5)
plot(x = 1:length(cum_gbc_count), y = cum_gbc_count / num_reads, col = "blue", xlab = "Number of groups", ylab = "Fraction of reads", main = TITLE)
dev.off()

q()

