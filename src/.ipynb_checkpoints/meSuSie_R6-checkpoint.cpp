#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::export]]
double loglik_cpp_R6(arma::vec V, const arma::mat& betahat, const arma::mat& shat2,const arma::vec& prior_weight, const int nancestry, arma::uvec diag_index){		
  mat Xty = betahat/shat2;
  
  
  mat cor_mat(nancestry,nancestry); 
  uvec upper_indices = trimatu_ind(size(cor_mat));
  cor_mat(upper_indices) = V;
  cor_mat = symmatu(cor_mat);
  cor_mat.diag().ones();
  
  mat se_mat = eye(nancestry,nancestry);
  V.elem(diag_index-1) =  sqrt(exp(V.elem(diag_index-1))) ;
  se_mat.diag() = V.elem(diag_index-1);
  
  arma::mat V_mat = se_mat*cor_mat*se_mat;
 // cout<<V_mat<<endl;
 // cout<<V_mat<<endl;
  vec lbf_multi(shat2.n_rows);
  double lbf_part_1,lbf_part_2;
  mat SIGMA,V_SIGMA,post_var,SIGMA_INV;
  rowvec Xty_var;
  //cout<<shat2.n_rows<<endl;
  for(int i=0;i<shat2.n_rows;i++){
    SIGMA_INV = diagmat(shat2.row(i));
    SIGMA = inv(SIGMA_INV);
    V_SIGMA = V_mat*SIGMA;
    // post_var = diagmat(shat2.row(i))-inv_sympd(SIGMA+SIGMA*V_SIGMA);
    post_var = SIGMA_INV-inv(SIGMA+SIGMA*V_SIGMA);
    Xty_var = Xty.row(i)*post_var;	
    lbf_part_1 = 0.5*dot(Xty_var.t(),Xty.row(i));
    lbf_part_2 = 0.5*log(det(V_SIGMA+eye(nancestry,nancestry)));
    lbf_multi(i) = lbf_part_1 - lbf_part_2;
  }
  
  double maxlbf,weighted_sum_w;
  vec lbf_1,lbf_2,w_multi,w_1,w_2,max_vec;
  vec mu = zeros<vec>(shat2.n_rows);

    maxlbf = lbf_multi.max();
    w_multi = exp(lbf_multi - maxlbf)%prior_weight;
    weighted_sum_w = accu(w_multi);
    
    
    Rcpp::Rcout << prior_weight << std::endl;    
    
  double lbf_model = maxlbf + log(weighted_sum_w);
  return(-1.0*lbf_model);
}
// [[Rcpp::export]]
SEXP mvlmm_reg(arma::mat betahat,arma::mat shat2, arma::mat V_mat){
  mat Xty = betahat/shat2;
  double nancestry = shat2.n_cols;
  vec lbf_multi(shat2.n_rows);
  mat post_mean_wmulti(shat2.n_rows,nancestry),post_mean2_wmulti(shat2.n_rows,nancestry);
  
  double lbf_part_1,lbf_part_2;
  mat SIGMA,V_SIGMA,post_var,SIGMA_INV;
  rowvec Xty_var;
  
  for(int i=0;i<shat2.n_rows;i++){
    SIGMA_INV = diagmat(shat2.row(i));
    SIGMA = inv(SIGMA_INV);
    V_SIGMA = V_mat*SIGMA;
    // post_var = diagmat(shat2.row(i))-inv_sympd(SIGMA+SIGMA*V_SIGMA);
    post_var = SIGMA_INV-inv(SIGMA+SIGMA*V_SIGMA);
    Xty_var = Xty.row(i)*post_var;	
    lbf_part_1 = 0.5*dot(Xty_var.t(),Xty.row(i));
    lbf_part_2 = 0.5*log(det(V_SIGMA+eye(nancestry,nancestry)));
    lbf_multi(i) = lbf_part_1 - lbf_part_2;
    post_mean_wmulti.row(i) = Xty_var;
    post_mean2_wmulti.row(i)=square(Xty_var)+post_var.diag().t();
  }
  List res = List::create(Named("lbf") = lbf_multi , Named("post_mean")=post_mean_wmulti,Named("post_mean2")=post_mean2_wmulti);
  
  return(res);

}
/*Note that mu1,mu2 are cube of dimension (p, L, Nancestry) */

// [[Rcpp::export]]
SEXP test_ELBO(arma::mat alpha,arma::cube mu1,arma::cube mu2,arma::cube XtX,arma::mat XtX_diag, arma::mat Xty,arma::vec yty,arma::vec N_vec ,arma::vec sigma2,double KL){
  
  arma::cube B_1 = mu1.each_slice()%alpha; //(p,L,nancestry)
  arma::cube B_2 = mu2.each_slice()%alpha; //(p,L,nancestry)
  arma::mat B_1_bar = sum(B_1,1); //Row sums of the elements returns p*nancestry matrix
  arma::mat B_2_bar = sum(B_2,1); //Row sums of the elements returns p*nancestry matrix
  
  vec BXXB(B_1.n_slices),BbarXXBbar(B_1.n_slices);
  for ( int i=0; i<B_1.n_slices; i++ ){
    BXXB(i) = accu(trans(B_1.slice(i))*XtX.slice(i)%trans(B_1.slice(i)));
    BbarXXBbar(i) =accu(trans(B_1_bar.col(i))*XtX.slice(i)*B_1_bar.col(i));
  }
  
  rowvec XXB2 = sum(XtX_diag%B_2_bar,0);
  rowvec BbarXty = 2*sum(B_1_bar%Xty,0);
  vec intermediate = yty - BbarXty.t()+BbarXXBbar+XXB2.t()-BXXB;
  vec sigma_update = intermediate/N_vec;
  double elbo =accu(-0.5*yty%log(2*datum::pi*sigma2) -0.5/sigma2%intermediate)-KL;
  List res = List::create(Named("ELBO") = elbo , Named("sigma2_update")=sigma_update);
  return(res);
}
// [[Rcpp::export]]
double test_run_loglik_cpp(arma::vec V, const arma::mat& betahat, const arma::mat& shat2,const arma::mat& prior_weight, const int nancestry, arma::uvec diag_index,Rcpp::List config_list){		

  mat Xty = betahat/shat2;
  mat cor_mat(nancestry,nancestry); 
  uvec upper_indices = trimatu_ind(size(cor_mat));
  cor_mat(upper_indices) = V;
  cor_mat = symmatu(cor_mat);
  cor_mat.diag().ones();
  
  mat se_mat = eye(nancestry,nancestry);
  V.elem(diag_index-1) =  sqrt(exp(V.elem(diag_index-1))) ;
  se_mat.diag() = V.elem(diag_index-1);
  
  arma::mat V_mat = se_mat*cor_mat*se_mat;

  arma::mat lbf_mat(shat2.n_rows,config_list.length());

 for(int n_config=0;n_config<config_list.length();n_config++){
	 vec config_vec = config_list[n_config];
	 config_vec = config_vec -1;
	 uvec config_vec_index = conv_to< uvec >::from( config_vec );
	 
	 arma::vec lbf_config(shat2.n_rows);
	 lbf_config.zeros(); 
	 vec mu = zeros<vec>(shat2.n_rows);
	 if(config_vec.n_elem==1){
		 lbf_config = arma::log_normpdf(betahat.col(as_scalar(config_vec)),mu,sqrt(as_scalar(V_mat(config_vec_index,config_vec_index)) + shat2.col(as_scalar(config_vec))))-arma::log_normpdf(betahat.col(as_scalar(config_vec)),mu,sqrt(shat2.col(as_scalar(config_vec))));
	 }else if(config_vec.n_elem!=1){
		 arma::mat V_sub = V_mat(config_vec_index,config_vec_index);
		 arma::mat shat2_sub = shat2.cols(config_vec_index);
		 arma::mat Xty_sub = Xty.cols(config_vec_index);

		 double lbf_part_1,lbf_part_2;
		 mat SIGMA,V_SIGMA,post_var,SIGMA_INV;
		 rowvec Xty_var_sub;
			for(int j=0;j<shat2.n_rows;j++){
							SIGMA_INV = diagmat(shat2_sub.row(j));
							SIGMA = inv(SIGMA_INV);
							V_SIGMA = V_sub*SIGMA;
							post_var = SIGMA_INV-inv(SIGMA+SIGMA*V_SIGMA);
							Xty_var_sub = Xty_sub.row(j)*post_var;	
							lbf_part_1 = 0.5*dot(Xty_var_sub.t(), Xty_sub.row(j));
							lbf_part_2 = 0.5*log(det(V_SIGMA+eye(V_sub.n_cols,V_sub.n_cols)));
							lbf_config(j) = lbf_part_1 - lbf_part_2;
			}	
	 }
	 lbf_mat.col(n_config) = lbf_config;
 }
 
 double maxlbf = lbf_mat.max();
 arma::mat w_multi = exp(lbf_mat - maxlbf)%prior_weight;
 double weighted_sum_w = accu(w_multi);
 double lbf_model = maxlbf + log(weighted_sum_w);
 //Rcpp::Rcout << prior_weight << std::endl;
 return(-1.0*lbf_model);

}

// [[Rcpp::export]]
SEXP test_run_mvlmm_reg(arma::mat betahat,arma::mat shat2, arma::mat V_mat){
  mat Xty = betahat/shat2;
  double nancestry = shat2.n_cols;
  vec lbf_multi(shat2.n_rows);
  mat post_mean_wmulti(shat2.n_rows,nancestry),post_mean2_wmulti(shat2.n_rows,nancestry);
  
  double lbf_part_1,lbf_part_2;
  mat SIGMA,V_SIGMA,post_var,SIGMA_INV;
  rowvec Xty_var;
  
  for(int i=0;i<shat2.n_rows;i++){
    SIGMA_INV = diagmat(shat2.row(i));
    SIGMA = inv(SIGMA_INV);
    V_SIGMA = V_mat*SIGMA;
    // post_var = diagmat(shat2.row(i))-inv_sympd(SIGMA+SIGMA*V_SIGMA);
    post_var = SIGMA_INV-inv(SIGMA+SIGMA*V_SIGMA);
    Xty_var = Xty.row(i)*post_var;	
    lbf_part_1 = 0.5*dot(Xty_var.t(),Xty.row(i));
    lbf_part_2 = 0.5*log(det(V_SIGMA+eye(nancestry,nancestry)));
    lbf_multi(i) = lbf_part_1 - lbf_part_2;
    post_mean_wmulti.row(i) = Xty_var;
    post_mean2_wmulti.row(i)=square(Xty_var)+post_var.diag().t();
  }
  List res = List::create(Named("lbf") = lbf_multi , Named("post_mean")=post_mean_wmulti,Named("post_mean2")=post_mean2_wmulti);
  
  return(res);

}

// [[Rcpp::export]]
arma::vec estimate_prior_annot(const arma::mat& annot_file_subset, const arma::vec& alpha) {
  int n_annot = annot_file_subset.n_cols;             // Number of annotations
  arma::vec w_new = arma::zeros<arma::vec>(n_annot);  // Initialize weights
  arma::vec prior_annot(annot_file_subset.n_rows);    // Prior annotation vector
  
  // Iterative optimization
  for (int iter = 0; iter < 100; ++iter) {
    arma::vec w_old = w_new; // Save the old weights
    
    // Loop over annotations
    for (int a = 0; a < n_annot; ++a) {
      // Compute weights (wt)
      arma::vec w_new_excluded = w_new;
      w_new_excluded.shed_row(a); // Remove the current annotation
      arma::mat annot_excluded = annot_file_subset;
      annot_excluded.shed_col(a); // Remove the current column
      
      arma::vec wt = arma::exp(annot_excluded * w_new_excluded);
      wt /= arma::sum(wt);
      
      // Compute k0, k1, r0, r1
      double k0 = arma::sum(wt % (1 - annot_file_subset.col(a)));
      double k1 = arma::sum(wt % annot_file_subset.col(a));
      double r0 = arma::sum(alpha % (1 - annot_file_subset.col(a)));
      double r1 = arma::sum(alpha % annot_file_subset.col(a));
      
      // Update weight
      w_new[a] = std::log(r1 / r0) - std::log(k1 / k0);
    }
    
    // Check convergence
    double eps = arma::sum(arma::square(w_new - w_old));
    if (eps < 1e-3) {
      break;
    }
  }
  
  // Compute prior annotations
  prior_annot = arma::exp(annot_file_subset * w_new);
  prior_annot /= arma::sum(prior_annot);
  
  return prior_annot;
}



// [[Rcpp::export]]
double llk_prior_causal_cpp(const arma::vec& w, 
                            const arma::mat& annotation, 
                            const arma::mat& alpha) {
  // Compute row sums of alpha
  arma::vec alpha_col = arma::sum(alpha, 1);
  
  // Compute the linear predictor: annotation * w
  arma::vec lin_pred = annotation * w;
  
  // For numerical stability, subtract the maximum value (log-sum-exp trick)
  double max_val = lin_pred.max();
  arma::vec exps = arma::exp(lin_pred - max_val);
  
  // Compute the softmax (i.e., prior) in a numerically stable way
  double sum_exps = arma::sum(exps);
  arma::vec prior = exps / sum_exps;
  
  // Optionally, clamp the prior to avoid log(0) issues, original 1e-20, changed to 1e-10 for pca
  double eps = 1e-10;
  // Using arma::clamp if available (otherwise, you can loop over elements)
  prior = arma::clamp(prior, eps, 1 - eps);
  
  // Compute the log-likelihood component-wise
  arma::vec llk = alpha_col % arma::log(prior); //+ (1 - alpha_col) % arma::log(1 - prior);
  
  //double total_llk = arma::sum(llk);
  //Rcpp::Rcout << "Total log-likelihood: " << total_llk << "\n";
  
  
  // Return the negative total log-likelihood
  return -arma::sum(llk);
}

// Optimizer wrapper using R's optim
// [[Rcpp::export]]
arma::vec optimize_llk(const arma::mat& annotation, const arma::mat& alpha, arma::vec initial_w) {
  // Load R's optim function
  Function optim("optim");

  // Define bounds for the parameters (adjust these as appropriate)
  NumericVector lower = rep(-5.0, initial_w.n_elem);
  NumericVector upper = rep(5.0, initial_w.n_elem);


  // Call optim with the llk_prior_causal_cpp function
  List result = optim(Named("par") = wrap(initial_w),
                      Named("fn") = Rcpp::InternalFunction(&llk_prior_causal_cpp),
                      Named("method") = "L-BFGS-B",  // Specify optimization method
		      Named("lower") = lower,
                      Named("upper") = upper,
		      Named("annotation") = wrap(annotation),
                      Named("alpha") = wrap(alpha));
  
  // Extract the optimized parameters
  NumericVector w_update = result["par"];
  
  // Return the optimized parameters as arma::vec
  return as<arma::vec>(w_update);
}
