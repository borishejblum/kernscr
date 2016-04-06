do_simulations <- function(num.sim = 100, with_feat = F, data_size = 100, ncol_gene_mat = 5, feat_m = NA, feat_d  = NA, num_pert = 1000, kernel = c("linear", "poly", "gaussian"), ...){

  res <- c()
  res.new <- c()

  for (i in 1:num.sim){
    data <- generate_dataset_advanced(with_feat = with_feat, data_size = data_size, ncol_gene_mat = ncol_gene_mat, feat_m = feat_m, feat_d  = feat_d)
    num <- scoreEval(data, c(5 + ncol_gene_mat,6 + ncol_gene_mat), matrix(0,data_size, data_size), kernel, ...)
    res <- c(res, sum(scoreEval_pert(1000, data, c(5 + ncol_gene_mat,6 + ncol_gene_mat), matrix(0,data_size, data_size), kernel, ...) >= as.numeric(num))/1000)
    num.new <- scoreEval_min_model(data, c(5 + ncol_gene_mat,6 + ncol_gene_mat), matrix(0,data_size, data_size), kernel, ...)
    res.new <- c(res.new, sum(scoreEval_pert_min_model(1000, data, c(5 + ncol_gene_mat,6 + ncol_gene_mat), matrix(0,data_size, data_size), kernel, ...) >= as.numeric(num.new))/1000)
  }

  print(sum(res <= 0.05)/num.sim)
  print(sum(res.new <= 0.05)/num.sim)

  return(list(res,res.new))
}