nc <- function(adj_matrix) {
	#获取 0-1 矩阵，1 表示节点间存在边，0 表示不存在边
	adj_matrix <- as.matrix(adj_matrix)
	adj_matrix[abs(adj_matrix) != 0] <- 1
	
	#矩阵的特征分解，获取特征值 λ
	lambda <- eigen(adj_matrix, only.values = TRUE)$values
	lambda <- sort(lambda, decreasing = TRUE)
	
	#计算“平均特征根”，获得自然连通度
	lambda_sum <- 0
	N = length(lambda)
	for (i in 1:N) lambda_sum = lambda_sum + exp(lambda[i])
	lambda_average <- log(lambda_sum/N, base = exp(1))
	lambda_average
}