
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

#include "utils.hpp"

template <class Type>
Type MMTMB(objective_function<Type>* obj) {
  
  using namespace density;
  //list data
  /*-------------------------------DATA--------------------------------------*/
  DATA_STRUCT(listZs, SparseMatrixList);
  DATA_STRUCT(listBs, SparseMatrixList);
  DATA_STRUCT(listPus, SparseMatrixList);
  DATA_STRUCT(listPvs, SparseMatrixList);

  // solecrop-focal
  DATA_VECTOR(y);                                                     //Vector of yield observations on solecrop
  DATA_MATRIX(X);                                                    // Design matrix for fixed effects on solecrop
  SparseMatrix<Type> Z_DBV = listZs.getName("Z_DBV");                // Design matrix for randoms effects [n x q] DBV
  SparseMatrix<Type> Z_SIGV = listZs.getName("Z_SIGV");              // Design matrix for  randoms effects  [n x q] SIGV
  //spatials data for focal plants | SC
  Eigen::SparseMatrix<Type> Bs_SC_f = listBs.getName("Bs_SC_f");
  Eigen::SparseMatrix<Type> Pu_SC_f = listPus.getName("Pu_SC_f");
  Eigen::SparseMatrix<Type> Pv_SC_f = listPvs.getName("Pv_SC_f");
  //params
  PARAMETER(log_sd_e);                         
  PARAMETER(log_sd_sigv);  
  PARAMETER_VECTOR(beta);                                    // vector of fixed effects
  PARAMETER_VECTOR(SIGV);                                   // vector of random effects [SIGV]
  PARAMETER_VECTOR(cs_SC_f);
  PARAMETER(log_lambda_v_sc_f);
  PARAMETER(log_lambda_u_sc_f);
  
  // solecrop tester
  DATA_VECTOR(y_t);                                               //Vector of yield observations on solecrop
  DATA_MATRIX(X_t);                                              // Design matrix for fixed effects on solecrop
  PARAMETER(log_sd_e_t);                          
  PARAMETER_VECTOR(beta_t);                                    // vector of fixed effects
  PARAMETER_VECTOR(cs_SC_t);
  //spatials data for testers
  Eigen::SparseMatrix<Type> Bs_SC_t = listBs.getName("Bs_SC_t");
  Eigen::SparseMatrix<Type> Pu_SC_t = listPus.getName("Pu_SC_t");
  Eigen::SparseMatrix<Type> Pv_SC_t = listPvs.getName("Pv_SC_t");
  PARAMETER(log_lambda_v_sc_t);
  PARAMETER(log_lambda_u_sc_t);

  // intercropping data
  DATA_ARRAY(Y);                                      // Matrix of observations [n x 2] of yield on intercropping
  DATA_MATRIX(X_mix);                                // Design matrix for fixed effects
  SparseMatrix<Type> Z_BV = listZs.getName("Z_BV");           // Design matrix for random effects [n_1 x q],[DBV, SBV]
  SparseMatrix<Type> Z_SIGV_mix = listZs.getName("Z_SIGV_mix");                         // Design matrix for  randoms effects SIGV
  SparseMatrix<Type> Z_DBV_x_SBV = listZs.getName("Z_DBV_x_SBV");                       //  Design matrix for interactions effects (DBV_w x IGE_p) +( DBV_p x IGE_w)
  //intercropping parameters
  PARAMETER_MATRIX(B);                          // Matrix of fixed effects  [p x 2] on intercropping
  PARAMETER_ARRAY(BV);                          // Matrix of random effects [q x 2][DBV, SBV]
  PARAMETER_ARRAY(DBV_x_SBV);                       //Matrix for interactions effects
  PARAMETER_VECTOR(log_sd_BV);                // vector of variances parameters for U
  PARAMETER_VECTOR(log_sd_DBV_x_SBV);             // vector of variances parameters for U_I
  PARAMETER_VECTOR(unconstr_cor_BV);        //Unconstraint correlation parameter for random effects U
  PARAMETER_VECTOR(unconstr_cor_DBV_x_SBV);      //Unconstraint correlation parameter U_I
  PARAMETER_VECTOR(log_sd_E);                //vector of error variances in intercropping
  PARAMETER_VECTOR(unconstr_cor_E);        // Unconstraint correlation for errors
  DATA_MATRIX(A);                            // Variance-covariance matrix for genetics effects [kinship]
  DATA_MATRIX(Amix);                       //Variance-covariance matrix effects in mixed

  //spatials data for focal plants | IC
  Eigen::SparseMatrix<Type> Bs_IC_f = listBs.getName("Bs_IC_f");
  Eigen::SparseMatrix<Type> Pu_IC_f = listPus.getName("Pu_IC_f");
  Eigen::SparseMatrix<Type> Pv_IC_f = listPvs.getName("Pv_IC_f");
  PARAMETER_VECTOR(cs_IC_f);
  //spatials data for testers plants | IC
  Eigen::SparseMatrix<Type> Bs_IC_t = listBs.getName("Bs_IC_t");
  Eigen::SparseMatrix<Type> Pu_IC_t = listPus.getName("Pu_IC_t");
  Eigen::SparseMatrix<Type> Pv_IC_t = listPvs.getName("Pv_IC_t");
  PARAMETER_VECTOR(cs_IC_t);
  //params for focal+testers
  PARAMETER(log_lambda_u_ic_f);
  PARAMETER(log_lambda_v_ic_f);
  
  PARAMETER(log_lambda_u_ic_t);
  PARAMETER(log_lambda_v_ic_t);

  //optional for SIMULATE with spatials
  DATA_INTEGER(sim_with_spats);
  
  /*----------------------------Local variales -------------------*/
  // size
  int d = 2;
  int q = A.cols(); 
  int n = y.size();
  int n_t = y_t.size();
  //lambdas focal
  Type lambda_u_sc_f = exp(log_lambda_u_sc_f);
  Type lambda_v_sc_f = exp(log_lambda_v_sc_f);
  Type lambda_u_ic_f = exp(log_lambda_u_ic_f);
  Type lambda_v_ic_f = exp(log_lambda_v_ic_f);
  //lambda testers
  Type lambda_u_sc_t = exp(log_lambda_u_sc_t);
  Type lambda_v_sc_t = exp(log_lambda_v_sc_t);
  Type lambda_u_ic_t = exp(log_lambda_u_ic_t);
  Type lambda_v_ic_t = exp(log_lambda_v_ic_t);
  
  /*---------------- DBV & fixed effects (local parameters)-------------------*/
  vector<Type> DBV = BV.matrix().col(0); //extract DBV from U => matrix of [DBV, SBV]
  
  /*----------------------------Objective function-----------------------------*/
  Type nll = Type(0.0);
  
  /*----------------------------contribution to the nll-----------------------*/
  //tester obs in sole crop
  if(!mmtmb::isEmpty(y_t)){
    vector<Type> cs_SC_t_sim(cs_SC_t.size());
    vector<Type> m_t(n_t);
    m_t  = X_t * beta_t;
    if(mmtmb::isSparseValid(Bs_SC_t)){
      m_t += Bs_SC_t*cs_SC_t;
    }
    SparseMatrix<Type> Id_t(n_t, n_t);
    mmtmb::fillMatId(Id_t);
    MVNORM_t<Type> mvn_y_cor_t(Id_t);
    vector<Type> vec_sd_t(n_t);
    vec_sd_t.fill(exp(log_sd_e_t));
    nll += VECSCALE(mvn_y_cor_t, vec_sd_t)(y_t - m_t);
    REPORT(beta_t);
    // contribution of the spatials effects | P-Splines
    if(mmtmb::isSparseValid(Bs_SC_t)){
      Eigen::SparseMatrix<Type> Q_SC_t(cs_SC_t.size(), cs_SC_t.size());
      Q_SC_t = lambda_u_sc_t * Pu_SC_t + lambda_v_sc_t * Pv_SC_t;
      for(int k = 0; k < cs_SC_t.size(); ++k) {
        Q_SC_t.coeffRef(k, k) += Type(1e-6);  
      }
      nll += GMRF(Q_SC_t)(cs_SC_t);
      //simulate spatials effects
      SIMULATE{
        GMRF(Q_SC_t).simulate(cs_SC_t_sim);
      }
    }
    vector<Type> e_sim_t(n_t);
    vector<Type> yobs_t(n_t);
    SIMULATE {
      VECSCALE_t<MVNORM_t<Type>> h_e_t = VECSCALE(mvn_y_cor_t, vec_sd_t);
      h_e_t.simulate(e_sim_t);
      yobs_t  = X_t * beta_t + e_sim_t;
      if(mmtmb::isSparseValid(Bs_SC_t)){
        if(sim_with_spats){
          yobs_t  += Bs_SC_t*cs_SC_t_sim;
        }
      }
      REPORT(yobs_t); 
    }
  }
  //contribution of U to the nll
  UNSTRUCTURED_CORR_t<Type> mvn_u_bv(unconstr_cor_BV);                                               
  vector<Type> sd_BV = exp(log_sd_BV);                                                              
  VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > f_bv = VECSCALE(mvn_u_bv, sd_BV);                           
  MVNORM_t<Type> g_bv(A);                                                                           
  SEPARABLE_t< VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > , MVNORM_t<Type> > h_BV(f_bv, g_bv);
  array<Type> BV_sim(BV.rows(), BV.cols());
  SIMULATE {
    h_BV.simulate(BV_sim);
    REPORT(BV_sim);
  }
  nll += h_BV(BV);
  
  // contribution of the spatials effects | P-Splines | IC
  vector<Type>cs_IC_f_sim (cs_IC_f.size());
  if(mmtmb::isSparseValid(Bs_IC_f)){
    Eigen::SparseMatrix<Type> Q_IC_f(cs_IC_f.size(), cs_IC_f.size());
    Q_IC_f = lambda_u_ic_f * Pu_IC_f + lambda_v_ic_f * Pv_IC_f;
    for(int k = 0; k < cs_IC_f.size(); ++k) {
      Q_IC_f.coeffRef(k, k) += Type(1e-6);  
    }
    nll += GMRF(Q_IC_f)(cs_IC_f);
    //simulate spatials effects
    SIMULATE{
      GMRF(Q_IC_f).simulate(cs_IC_f_sim);
    }
  }
  vector<Type>cs_IC_t_sim (cs_IC_t.size());
  if(mmtmb::isSparseValid(Bs_IC_t)){
    Eigen::SparseMatrix<Type> Q_IC_t(cs_IC_t.size(), cs_IC_t.size());
    Q_IC_t = lambda_u_ic_t * Pu_IC_t + lambda_v_ic_t * Pv_IC_t;
    for(int k = 0; k < cs_IC_t.size(); ++k) {
      Q_IC_t.coeffRef(k, k) += Type(1e-6);  
    }
    nll += GMRF(Q_IC_t)(cs_IC_t);
    //simulate spatials effects
    SIMULATE{
      GMRF(Q_IC_t).simulate(cs_IC_t_sim);
    }
  }
  
  //focal obs in sole crop with sigv
  vector<Type> SIGV_sim(SIGV.size());
  if(!mmtmb::isEmpty(y)){
    vector<Type> cs_SC_f_sim(cs_SC_f.size());
    vector<Type> m(n);
    m  = X * beta + Z_DBV*DBV + Z_SIGV*SIGV;
    //add spatial effects for focal
    if(mmtmb::isSparseValid(Bs_SC_f)){
      m += Bs_SC_f*cs_SC_f ;
    }
    SparseMatrix<Type> Id(n, n);
    mmtmb::fillMatId(Id);
    MVNORM_t<Type> mvn_y_cor(Id);
    vector<Type> vec_sd(n);
    vec_sd.fill(exp(log_sd_e));
    nll += VECSCALE(mvn_y_cor, vec_sd)(y - m);
    vector<Type> e_sim(y.size());
    SIMULATE {
      VECSCALE_t<MVNORM_t<Type>> h_e = VECSCALE(mvn_y_cor, vec_sd);
      h_e.simulate(e_sim);
    }
    // contribution of the SIGV to the nll
    MVNORM_t<Type> mvn_sigv_cor(A);
    vector<Type> sd_sigv(q);
    sd_sigv.fill(exp(log_sd_sigv));
    nll += VECSCALE(mvn_sigv_cor, sd_sigv)(SIGV);
    SIMULATE {
      VECSCALE_t<MVNORM_t<Type>> h_SIGV = VECSCALE(mvn_sigv_cor, sd_sigv);
      h_SIGV.simulate(SIGV_sim);
      REPORT(SIGV_sim);
    }
    // contribution of the spatials effects | P-Splines
    if(mmtmb::isSparseValid(Bs_SC_f)){
      Eigen::SparseMatrix<Type> Q_SC_f(cs_SC_f.size(), cs_SC_f.size());
      Q_SC_f = lambda_u_sc_f * Pu_SC_f + lambda_v_sc_f * Pv_SC_f;
      for(int k = 0; k < cs_SC_f.size(); ++k) {
        Q_SC_f.coeffRef(k, k) += Type(1e-6);  
      }
      nll += GMRF(Q_SC_f)(cs_SC_f);
      //simulate spatials effects
      SIMULATE{
        GMRF(Q_SC_f).simulate(cs_SC_f_sim);
      }
    }
    //simulate yobs
    vector<Type> yobs(n);
    SIMULATE{
      vector<Type> DBV_sim = BV_sim.matrix().col(0);
      yobs  = X * beta + Z_DBV*DBV_sim + Z_SIGV*SIGV_sim + e_sim;
      if(mmtmb::isSparseValid(Bs_SC_f)){
        if(sim_with_spats){
          yobs  += Bs_SC_f*cs_SC_f_sim;
        }
      }
      REPORT(yobs);
    }
    REPORT(beta);
  }
  //contribution of DBV_x_SBV to the nll
  array<Type> DBV_x_SBV_sim(DBV_x_SBV.rows(),  DBV_x_SBV.cols());
  if(mmtmb::isSparseValid(Z_DBV_x_SBV)){
    UNSTRUCTURED_CORR_t<Type> mvn_u_dbv_x_sbv(unconstr_cor_DBV_x_SBV);     
    vector<Type> sd_DBV_x_SBV = exp(log_sd_DBV_x_SBV);                
    VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > f_dbv_x_sbv = VECSCALE(mvn_u_dbv_x_sbv, sd_DBV_x_SBV);
    MVNORM_t<Type> g_dbv_x_sbv(Amix);                                   
    SEPARABLE_t< VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > , MVNORM_t<Type> > h_DBV_x_SBV(f_dbv_x_sbv, g_dbv_x_sbv);
    SIMULATE {
      h_DBV_x_SBV.simulate(DBV_x_SBV_sim);
      REPORT(DBV_x_SBV_sim);
    }
    nll += h_DBV_x_SBV(DBV_x_SBV);
    matrix<Type> Cor_DBV_x_SBV(d,d);
    Cor_DBV_x_SBV = mvn_u_dbv_x_sbv.cov();
    REPORT(Cor_DBV_x_SBV);
  }
  
  //contribution of the intercropping observation to the nll
  UNSTRUCTURED_CORR_t<Type> mvn_y(unconstr_cor_E);                                          
  vector<Type> sd_E = exp(log_sd_E);                                                       
  VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > f_y = VECSCALE(mvn_y, sd_E);
  matrix<Type> Id_mix(Y.rows(), Y.rows());
  mmtmb::fillMatId(Id_mix);
  MVNORM_t<Type> g_y(Id_mix); 
  SEPARABLE_t< VECSCALE_t<UNSTRUCTURED_CORR_t<Type> > , MVNORM_t<Type> > h_Y(f_y, g_y);
  array<Type> E_sim(Y.rows(),  Y.cols());
  SIMULATE {
    h_Y.simulate(E_sim);
  }
  matrix<Type> N(Y.rows(), d);
  N = X_mix*B + Z_BV*BV.matrix();
  if(mmtmb::isSparseValid(Z_DBV_x_SBV)){
    N += Z_DBV_x_SBV*DBV_x_SBV.matrix();
  }
  if(mmtmb::isSparseValid(Bs_IC_f)){
    vector<Type>spatials_effects_IC_f = Bs_IC_f*cs_IC_f;
    for (int i = 0; i < N.rows(); i++) {
      N(i, 0) += spatials_effects_IC_f(i); //first cols
    }
  }
  if(mmtmb::isSparseValid(Bs_IC_t)){
    vector<Type>spatials_effects_IC_t = Bs_IC_t*cs_IC_t;
    for (int i = 0; i < N.rows(); i++) {
      N(i, 1) += spatials_effects_IC_t(i); //2nd cols
    }
  }
  if(!mmtmb::isEmpty(y)){
    vector<Type> m_sigv = Z_SIGV_mix*SIGV;
    for (int i = 0; i < N.rows(); i++) {
      N(i, 0) +=  m_sigv(i);
    }
  }
  nll += h_Y(Y - N.vec());

  //simulate Yobs
  array<Type> Yobs(Y.rows(), d);
  SIMULATE {
    matrix<Type> etaY(Y.rows(), d);
    etaY = X_mix * B + Z_BV * BV_sim.matrix();
    if(mmtmb::isSparseValid(Z_DBV_x_SBV)){
      etaY +=  Z_DBV_x_SBV * DBV_x_SBV_sim.matrix();
    }
    if(mmtmb::isSparseValid(Bs_IC_f)){
      vector<Type>spatials_effects_IC_f_sim = Bs_IC_f*cs_IC_f_sim;
      for (int i = 0; i < etaY.rows(); i++) {
        if(sim_with_spats){
          etaY(i, 0) += spatials_effects_IC_f_sim(i);
        }
      }
    }
    if(mmtmb::isSparseValid(Bs_IC_t)){
      vector<Type>spatials_effects_IC_t_sim = Bs_IC_t*cs_IC_t_sim;
      for (int i = 0; i < etaY.rows(); i++) {
        if(sim_with_spats){
          etaY(i, 1) += spatials_effects_IC_t_sim(i);
        }
      }
    }
    if(!mmtmb::isEmpty(y)){
      matrix<Type> M_sim(Y.rows(), d);
      vector<Type> m_sigv_sim = Z_SIGV_mix * SIGV_sim;
      for (int i = 0; i < M_sim.rows(); i++) { //simplifier avec etaY +=
        M_sim(i, 0) = etaY(i, 0) + m_sigv_sim(i);
        M_sim(i, 1) = etaY(i, 1);
      }
      for (int i = 0; i < Yobs.rows(); i++) {
        for (int j = 0; j < Yobs.cols(); j++) {
          Yobs(i, j) = M_sim(i, j) + E_sim.matrix()(i, j);
        }
      }
    } else {
      for (int i = 0; i < Yobs.rows(); i++) {
        for (int j = 0; j < Yobs.cols(); j++) {
          Yobs(i, j) = etaY(i, j) + E_sim.matrix()(i, j);
        }
      }
      
    }
    REPORT(Yobs);
  }
  
  /*----------------------------report correlation matrix and effects-----------------*/
  matrix<Type> Cor_BV(d,d);
  Cor_BV = mvn_u_bv.cov();
  REPORT(Cor_BV);
  
  matrix<Type> Cor_E(d,d);
  Cor_E = mvn_y.cov();
  REPORT(Cor_E);
  
  return nll;
}

/*-------------------------univariate analysis (solecrop)----------------------------*/
template <class Type>
Type lmm(objective_function<Type>* obj) {
  
  using namespace density;
  
  // Receive list of sparse matrices from R
  DATA_STRUCT(listZs, SparseMatrixList);
  
  /*-------------------------------Data--------------------------------------*/
  DATA_VECTOR(y);                       // Vector of yield observations on solecrop
  DATA_MATRIX(X);                       // Design matrix for fixed effects on solecrop
  SparseMatrix<Type> Z_DBV = listZs.getName("Z_DBV"); // Design matrix for random effects cBV
  DATA_MATRIX(A);                       // dummy matrix => not used
  DATA_SPARSE_MATRIX(Bs);
  DATA_SPARSE_MATRIX(Pu);
  DATA_SPARSE_MATRIX(Pv);

  //optional for SIMULATE with spatials
  DATA_INTEGER(sim_with_spats);
  
  
  /*-----------------------------------Parameters-----------------------------*/
  PARAMETER_VECTOR(beta);               // Vector of fixed effects
  PARAMETER_VECTOR(u);                  // Vector of random effects
  PARAMETER(log_sd_u);                  // Log-standard deviation for random effects
  PARAMETER(log_sd_e);                  // Log-standard deviation for residuals
  PARAMETER_VECTOR(theta);
  PARAMETER(log_lambda_u);
  PARAMETER(log_lambda_v);
  
  /*----------------------------Local variables -----------------------------*/
  int n = y.size();                           // Length of y
  int q = u.size();                          // Number of random parameters
  Type lambda_u = exp(log_lambda_u);
  Type lambda_v = exp(log_lambda_v);
  
  
  /*----------------------------Objective function-----------------------------*/
  Type nll = Type(0.0);                 // Initialize negative log-likelihood
  
  /*----------------------------Contribution to the nll-----------------------*/
  vector<Type> m(n);
  m = X * beta + Z_DBV * u;
  if(mmtmb::isSparseValid(Bs)){
    m += Bs*theta;
  }
  SparseMatrix<Type> Id(n, n);
  mmtmb::fillMatId(Id);
  MVNORM_t<Type> mvn_y_cor(Id);
  vector<Type> vec_sd(n);
  vec_sd.fill(exp(log_sd_e));
  nll += VECSCALE(mvn_y_cor, vec_sd)(y - m);
  vector<Type> e_sim(y.size());
  SIMULATE {
    VECSCALE_t<MVNORM_t<Type>> h_e = VECSCALE(mvn_y_cor, vec_sd);
    h_e.simulate(e_sim);
  }
  
  vector<Type>theta_sim(theta.size());
  if(mmtmb::isSparseValid(Bs)){
    //Spline penalty via GMRF
    Eigen::SparseMatrix<Type> Q(theta.size(), theta.size());
    Q = lambda_u * Pu + lambda_v * Pv;
    for(int k = 0; k < theta.size(); ++k) {
      Q.coeffRef(k, k) += Type(1e-6);  
    }
    nll += GMRF(Q)(theta);
    SIMULATE {
      GMRF(Q).simulate(theta_sim);
    }
  }
  // Contribution of random effects
  vector<Type> unconstr_cor_A(q*(q-1)/2);
  unconstr_cor_A.fill(0.01);
  vector<Type> sd_u(q);                 
  sd_u.fill(exp(log_sd_u));
  nll += VECSCALE(UNSTRUCTURED_CORR(unconstr_cor_A), sd_u)(u);
  vector<Type> u_sim(u.size());
  vector<Type> yobs(n);
  SIMULATE {
    VECSCALE(UNSTRUCTURED_CORR(unconstr_cor_A), sd_u).simulate(u_sim);
    yobs  = X * beta + Z_DBV*u_sim + e_sim;
    if(mmtmb::isSparseValid(Bs)){
      if(sim_with_spats){
        yobs += Bs*theta_sim;
      }
    }
  }
  //report simulate
  SIMULATE{
    REPORT(u_sim);
    REPORT(yobs);
    REPORT(theta_sim);
  }
  
  /*---------------------------Return nll-----------------------------*/
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
