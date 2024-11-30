alpha <- function(x, tree = NULL, base = exp(1)) {
	est <- estimateR(x)
	Richness <- est[1, ]
	Chao1 <- est[2, ]
	ACE <- est[4, ]
	Shannon <- diversity(x, index = 'shannon', base = base)
	Simpson <- diversity(x, index = 'simpson')	#Gini-Simpson 指数
	Pielou <- Shannon / log(Richness, base)
	goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
	result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
	if (!is.null(tree)) {
		PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
		names(PD_whole_tree) <- 'PD_whole_tree'
		result <- cbind(result, PD_whole_tree)
	}
	result
}