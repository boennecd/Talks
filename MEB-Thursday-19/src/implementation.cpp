#include <TMB.hpp>
#include <cmath>
#include <algorithm> 
#include "gausHermite.h"
#include <future>
#include <limits>

namespace atomic {
/* Computes log CDF of standard normal distribution */
TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  pnorm_log1
  ,
  // OUTPUT_DIM
  1
  ,
  // ATOMIC_DOUBLE
  ty[0] = Rmath::Rf_pnorm5(tx[0], 0, 1, 1, 1);
  ,
  // ATOMIC_REVERSE
  Type const cdf = exp(ty[0]);
  px[0] = dnorm1(tx[0]) / cdf * py[0];
  )
} // namespace atomic

/* Computes the log CDF of normal distribution. 
 * 
 * Args: 
 *   Similar to `stats::pnorm`.
 */
template<class Type>
Type pnorm_log(Type q, Type  mean = 0., Type sd = 1.){
  CppAD::vector<Type> tx(1);
  tx[0] = (q - mean) / sd;
  return atomic::pnorm_log1(tx)[0];
}
VECTORIZE3_ttt(pnorm_log)
VECTORIZE1_t  (pnorm_log)
  
/* quantity needed when mapping from CP to DP w/ SNVA. R functions are 
 func <- function(g){
 g_abs <- abs(g)
 cv <- 2 * g_abs / (4 - pi)
 out <-  cv^(1/3) / sqrt(1 + cv^(2/3))
 ifelse(g < 0, -out, out)
 }
 dfunc <- function(g){
 cv <- 2 * g / (4 - pi)
 cv_2_3 <- (cv * cv)^(1/3)
 1 / (3 * (cv_2_3 + 1)^(3/2) * cv_2_3) * 2 / (4 - pi)
 }
 plot(func, xlim = c(-.9, .9))
 gs <- seq(-.9, .9, length.out = 100)
 all.equal(sapply(gs, numDeriv::grad, func = func),
 dfunc(gs))
 */
inline double gamma_to_nu_func(double const g){
  constexpr double mult = 2. / (4. - M_PI);
  double const g_sign = g < 0 ? -1. : 1., 
               g_abs  =  g * g_sign, 
               cv     = mult * g_abs, 
               out    = pow(cv, 1. / 3.) / sqrt(1. + pow(cv, 2. / 3.));
  return g_sign * out;
}

template<class Type>
Type dgamma_to_nu_func(Type const g){
  Type const mult(2. / (4. - M_PI)),
               cv = mult * g, 
           cv_2_3 = pow(cv * cv, 1. / 3.),
            denom = 3. * pow(cv_2_3 + 1., 3. / 2.) * cv_2_3, 
              out = mult / denom, 
              eps = Type(1. / std::numeric_limits<double>::epsilon());
  
  /* TODO: Bad solution to Inf grad at zero? */
  return CppAD::CondExpLe(out, eps, out, eps);
}

namespace atomic {
TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  gamma_to_nu1
  ,
  // OUTPUT_DIM
  1
  ,
  // ATOMIC_DOUBLE
  ty[0] =  gamma_to_nu_func(tx[0]);
  ,
  // ATOMIC_REVERSE
  px[0] = dgamma_to_nu_func(tx[0]) * py[0];
  )
} // namespace atomic

template<class Type>
Type gamma_to_nu(Type g){
  CppAD::vector<Type> tx(1);
  tx[0] = g;
  return atomic::gamma_to_nu1(tx)[0];
}
VECTORIZE3_ttt(gamma_to_nu)
VECTORIZE1_t  (gamma_to_nu)

/* TODO: check/test */
inline unsigned get_rng_dim(unsigned const n_vcov_params){
  double const n(n_vcov_params);
  return std::lround(std::sqrt(.25 + 2 * n) - .5);
}
template<class Type>
unsigned get_rng_dim(vector<Type> x){
  return get_rng_dim(x.size());
}

/* This function returns the matrix S = R^\top R given an upper triangular 
 * matrix R in column major order */
template<class Type>
matrix<Type> get_vcov_from_trian(Type const *t, unsigned const dim){
  matrix<Type> R(dim, dim);
  R.setZero();
  for(unsigned j = 0; j < dim; ++j)
    for(unsigned i = 0; i <= j; ++i)
      R(i, j) = *t++;
  
  return R.transpose() * R;
}

template<class Type>
matrix<Type> get_vcov_from_trian(vector<Type> const &theta){
  unsigned const dim = get_rng_dim(theta);
  return get_vcov_from_trian(&theta[0], dim);
}

/* TODO: check/test */
template<class Type> 
Type mult_var_dens(vector<Type> const &theta, matrix<Type> const &rngs){
  if(theta.size() < 2L){ /* univariate case */
    Type const sd = exp(theta[0]), 
             half(.5);
    Type out(0);
    Type const * const end = rngs.data() + rngs.size();
    for(Type const *x = rngs.data(); x != end; ++x){
      Type const scaled = *x / sd;
      out += - half * scaled * scaled;
    }
    
    out += - Type(rngs.size()) * log(Type(sqrt(2 * M_PI)) * sd);
    
    return out;
  }
  
  matrix<Type> const vcov = get_vcov_from_trian(theta);
  density::MVNORM_t<Type> norm_dist(vcov);
  Type out(0);
  for(unsigned g = 0; g < rngs.cols(); ++g)
    out -= norm_dist(rngs.col(g));
  
  return out;
}

/* Common arguments used in approximations. 
 * 
 * Args: 
 *   tobs: observed times. 
 *   event: zero/one variables for whether the event is observed.
 *   X: fixed effect design matrix. 
 *   XD: derivative of X wrt time. 
 *   Z: random effect design matrix. 
 *   grp: integer vector with group identifier.
 *   eps: small tolerance for positivity constraint on hazards.
 *   kappa: one-sided L2 penalty for positivity constraint.
 *   b: fixed effect coefficients. 
 *   theta: upper triangular matrix R in column major order such that the 
 *          covariance matrix is R^\top R.
 */
#define COMMON_ARGS                                            \
  vector<Type> const &tobs, vector<Type> const &event,         \
  matrix<Type> const &X, matrix<Type> const &XD,               \
  matrix<Type> const &Z, vector<int> const &grp,               \
  Type const &eps, Type const &kappa, vector<Type> const &b,   \
  vector<Type> const &theta, std::string const &link
#define COMMON_CALL                                            \
  tobs, event, X, XD, Z, grp, eps, kappa, b, theta, link

namespace {
/* Computes the log-likelihood for given random effects. 
 * 
 * Args: 
 *   u: [random effect dim] x [# groups] matrix with random effects. 
 */
template<class Type> 
Type laplace(COMMON_ARGS, matrix<Type> const &u){
  /* checks */
  unsigned const rng_dim = get_rng_dim(theta);
  {
    int const *max_grp = std::max_element(grp.data(), grp.data() + grp.size());
    if(!max_grp or *max_grp != u.cols() - 1L)
      error("Invalid 'grp' or 'u'");  
    if(rng_dim != u.rows())
      error("Invalid 'u'");
  }
  
  /* log-likelihood terms from conditional distribution of the observed 
   * outcomes */
  vector<Type> const eta = ([&](){
    vector<Type> out = X * b;
    for(unsigned i = 0; i < grp.size(); ++i){
      auto const g = grp[i];
      out[i] += Z.row(i) * u.col(g);
    }
    
    return out;
  })();
  vector<Type> const etaD    = XD * b;
         Type  const eps_log = log(eps);
  
  /* compute terms from conditional density */
  Type out(0);
  if(link == "PH")
    for(unsigned i = 0; i < grp.size(); ++i){
      Type const H = exp(eta[i]), 
                 h = etaD[i] * H, 
            if_low = event[i] * eps_log - H - h * h * kappa,
            if_ok  = event[i] * log(h)  - H;
      out += CppAD::CondExpGe(h, eps, if_ok, if_low);
    }
  else if (link == "PO"){
    Type const too_large(30.),
                     one(1.);
    for(unsigned i = 0; i < grp.size(); ++i){
      Type const H = CppAD::CondExpGe(
        eta[i], too_large, eta[i], log(one + exp(eta[i]))), 
                 h = etaD[i] * exp(eta[i] - H),
            if_low = event[i] * eps_log - H - h * h * kappa,
            if_ok  = event[i] * log(h)  - H;
      out += CppAD::CondExpGe(h, eps, if_ok, if_low);
    }
  } else if(link == "probit"){
    Type const tiny(std::numeric_limits<double>::epsilon()), 
               zero(0.), 
                one(1.);
    for(unsigned i = 0; i < grp.size(); ++i){
      Type const H = -pnorm_log(-eta[i]), 
                 h = etaD[i] * dnorm(-eta[i], zero, one) / 
                   (pnorm(-eta[i]) + tiny),
            if_low = event[i] * eps_log - H - h * h * kappa,
            if_ok  = event[i] * log(h)  - H;
      out += CppAD::CondExpGe(h, eps, if_ok, if_low);
      
    }
  } else 
    error("'%s' not implemented", link.c_str());
  
  /* log-likelihood terms from random effect density */
  out += mult_var_dens(theta, u);
  
  return -out;
}

/* 
 Makes an approximation of 
 \begin{align*}
 l(\mu,\sigma) &= 
 \frac 1{\sigma\sqrt{2\pi}}\int
 \exp \left(-\frac{(x-\mu)^2}{2\sigma^2} \right)
 \log\left(1 + \exp(x) \right)dx \\
 &\overset{x = \mu + \sqrt 2\sigma z}{=}
 \frac 1{\sqrt\pi}\int\exp(-z^2)
 \log\left(1 + \exp(\mu + \sqrt 2 \sigma z) \right)dz
 \end{align*}
 
 with (non-adaptive) Gauss–Hermite quadrature
 */
template<class Type>
Type GVA_mlogit_integral
  (Type const mu, Type const sigma, 
   vector<Type> const &GH_x, vector<Type> const &GH_w){
  unsigned const n_nodes = GH_x.size();
  
  Type const mult_sum(sqrt(M_1_PI)), 
                 mult(Type(M_SQRT2) * sigma),
            too_large(30.), 
                 one(1.);
  Type out(0.);
  for(unsigned i = 0; i < n_nodes; ++i){
    Type const eta = mu + mult * GH_x[i];
    out += 
      GH_w[i] * CppAD::CondExpGe(eta, too_large, eta, log(one + exp(eta)));
  }
  
  return mult_sum * out;
}

/* 
 Makes an approximation of 
 \begin{align*}
 l(\mu,\sigma, k, s) & = 
 \frac 1{\sigma\sqrt{2\pi}}\int
 \exp \left(-\frac{(x-\mu)^2}{2\sigma^2} \right)
 \log\left(1 + k \exp(s x) \right)dx \\ 
 &\overset{x = (-\log k+z)/s}{=}
 \frac 1{\lvert s\rvert\sigma\sqrt{2\pi}}\int
 \exp \left(-\frac{(x-\log k-s\mu)^2}{2s^2\sigma^2} \right)
 \log\left(1 + \exp(z) \right)dx
 \end{align*}
 */
template<class Type> 
Type GVA_mlogit_integral
  (Type const mu, Type const sigma, Type const log_k, Type const s,
   vector<Type> const &GH_x, vector<Type> const &GH_w){
  Type const mu_use = s * mu + log_k,
            sig_use = sqrt(s * s) * sigma;
  return GVA_mlogit_integral(mu_use, sig_use, GH_x, GH_w);
}

/* Makes an approximation of 
 l(\mu,\sigma) = 
 \int\phi(x;\mu,\sigma^2)
 (-\log \Phi(x))dx
 */
template<class Type> 
Type GVA_probit_integral
  (Type const mu, Type const sigma,
   vector<Type> const &GH_x, vector<Type> const &GH_w){
  
  Type out(0.);
  Type const mult_sum(sqrt(M_1_PI)), 
                 mult(Type(M_SQRT2) * sigma);
  for(unsigned i = 0; i < GH_x.size(); ++i)
    out += GH_w[i] * (-pnorm_log(mu + mult * GH_x[i]));
  
  return mult_sum * out;
}

/* Makes an approximation of  
 \begin{align*}
 l(\mu,\sigma, k, s) & = 
 \int\phi(x;\mu,\sigma^2)
 (-\log \Phi(k + sx))dx \\
 &= \int\phi(z;\underbrace{s\mu +k}_{\tilde\mu},
 \underbrace{s^2\sigma^2}_{\tilde\sigma^2})
 (-\log \Phi(z))dz
 \end{align*}
 */
template<class Type> 
Type GVA_probit_integral
  (Type const mu, Type const sigma, Type const k, Type const s,
   vector<Type> const &GH_x, vector<Type> const &GH_w){
  return GVA_probit_integral(mu * s + k, sqrt(s * s) * sigma, 
                             GH_x, GH_w);
}

#define GVA_COND_DENS_ARGS                                     \
  Type const eta_fix, Type const etaD_fix, Type const z,       \
  Type const event, Type const va_mean, Type const va_sd,      \
  Type const va_var

/* util class to compute the conditional density of the observed outcomes */
template<class Type>
class GVA_cond_dens {
  /* input dependend objects */
  Type const eps, 
         eps_log = Type(log(eps)), 
           kappa;
  vector<Type> const &GH_x,
                     &GH_w;
  
  /* potentially needed constants */
  Type const mlog_2_pi_half = Type(-log(2 * M_PI) / 2.), 
                        two = Type(2.);
public:
  GVA_cond_dens
  (Type const eps, Type const kappa, vector<Type> const &GH_x, 
   vector<Type> const &GH_w): 
  eps(eps), kappa(kappa), GH_x(GH_x), GH_w(GH_w) { }
  
  /* computes the conditional density term for the PH (log-log) link 
   * function */
  Type ph(GVA_COND_DENS_ARGS) const {
    Type const eta = eta_fix  + va_mean * z, 
                 h = etaD_fix * exp(eta),
          z_scaled = z * va_sd,
                 H = exp(eta + z_scaled * z_scaled / two), 
            if_low = event * eps_log - H - h * h * kappa,
            if_ok  = event * log(h)  - H;
    
    return CppAD::CondExpGe(h, eps, if_ok, if_low);
  }
  
  /* computes the conditional density term for the PO (-logit) link 
   * function */
  Type po(GVA_COND_DENS_ARGS) const {
    Type const eta = eta_fix + va_mean * z,
                 H = GVA_mlogit_integral(
                   va_mean, va_sd, eta_fix, z, GH_x, GH_w),
                 h = etaD_fix * exp(eta - H), 
            if_low = event * eps_log - H - h * h * kappa,
            if_ok  = event * log(h)  - H;
    
    return CppAD::CondExpGe(h, eps, if_ok, if_low);
  }
  
  /* computes the conditional density term for the probit (-probit) link 
   * function */
  Type probit(GVA_COND_DENS_ARGS) const {
    Type const H = GVA_probit_integral(
      va_mean, va_sd, -eta_fix, -z, GH_x, GH_w), 
            diff = eta_fix + z * va_mean,
               h = etaD_fix * exp(
                 mlog_2_pi_half - diff * diff / two - 
                   z * z * va_var / two + H), 
          if_low = event * eps_log - H - h * h * kappa,
          if_ok  = event * log(h)  - H;
    
    return CppAD::CondExpGe(h, eps, if_ok, if_low);
  }
};

/* Computes the lower bound. VA_1D is a specialized case for a 
 * one-dimensional shared effect for each group while VA_MD is for a 
 * multidimensional.
 * 
 * Args: 
 *   theta_VA: vector with VA parameters for each group in terms. For each 
 *             group, the first elements are the mean and the subsequent 
 *             elemenet is an upper diagonal matrix L such that the  
 *             covariance matrix is L^\top L.
 */
template<class Type> 
Type GVA_1D(COMMON_ARGS, vector<Type> const &theta_VA, 
            vector<Type> const &GH_x, vector<Type> const &GH_w){
  /* variational parameters for the means and standard deviations */
  auto get_means_or_sds = [&](bool const get_sd){
    unsigned const shift = get_sd ? 1L : 0L, 
                  n_elem = theta_VA.size() / 2L;
    vector<Type> out(n_elem);
    for(unsigned i = 0; i < n_elem; ++i)
      out[i] = theta_VA[i * 2L + shift];
    
    return out;
  };
  vector<Type> const va_means   = get_means_or_sds(false),
                     va_sds_log = get_means_or_sds(true),
                     va_sds     = exp(va_sds_log), 
                     va_var     = va_sds * va_sds;
  
  /* get random effect variance, fixed effect linear predictors, and 
   * assign constants */
  Type const sigma_log = theta[0], 
             sigma     = exp(sigma_log), 
             half(.5);
  vector<Type> const eta_fix  = X  * b,
                     etaD_fix = XD * b;
  
  /* handle terms from conditional density of observed outcomes */
  Type out(0);
  GVA_cond_dens<Type> const cond_dens(eps, kappa, GH_x, GH_w);
  
#define ADD_COND_DENS(func)                                    \
  for(unsigned i = 0; i < grp.size(); ++i){                    \
    auto const g = grp[i];                                     \
    out += cond_dens.func(                                     \
      eta_fix[i], etaD_fix[i], Z(i, 0), event[i], va_means[g], \
      va_sds[g], va_var[g]);                                   \
  }
  
  if(link == "PH")
    ADD_COND_DENS(ph)
  else if (link == "PO")
    ADD_COND_DENS(po)
  else if (link == "probit")
    ADD_COND_DENS(probit)
  else 
    error("'%s' not implemented", link.c_str());

#undef ADD_COND_DENS
  
  /* handle terms from random effect log density and VA log density */
  {
    Type var_sum(0), mu_sq_sum(0), sds_log_sum(0);
    for(unsigned g = 0; g < va_means.size(); g++){
      var_sum     += va_sds    [g] * va_sds  [g];
      mu_sq_sum   += va_means  [g] * va_means[g];
      sds_log_sum += va_sds_log[g];
      
    }
    Type const var = sigma * sigma, 
          n_groups(va_means.size());
    
    out += 
      - n_groups * sigma_log + sds_log_sum
      + half * (
          - var_sum   / var  
          - mu_sq_sum / var
          + n_groups);
  }
    
  return -out;
}

template<class Type> 
Type GVA_MD(COMMON_ARGS, vector<Type> const &theta_VA, 
            vector<Type> const &GH_x, vector<Type> const &GH_w, 
            unsigned const rng_dim){
  using Eigen::Dynamic;
  using vecT = vector<Type>;
  using std::move;
  
  /* get objects related to model covariance matrix */
  matrix<Type> const vcov = get_vcov_from_trian(theta);
  
  /* maybe not the best idea to matrix multiply by the precision matrix 
   * instead of using solve... */
  Type log_det_vcov;
  matrix<Type> vcov_inv;
  vcov_inv = atomic::matinvpd(vcov, log_det_vcov);
  
  /* get objects from VA distribution */
  unsigned const n_groups = theta_VA.size() / (
    rng_dim + (rng_dim * (rng_dim + 1L)) / 2L);
  std::vector<vecT >         va_means;
  std::vector<matrix<Type> > va_vcovs;
  va_means.reserve(n_groups);
  va_vcovs.reserve(n_groups);
  {
    Type const *t = &theta_VA[0];
    for(unsigned g = 0; g < n_groups; ++g){
      /* insert new mean vector */
      vecT mean_vec(rng_dim);
      for(unsigned i = 0; i < rng_dim; ++i)
        mean_vec[i] = *t++;
      va_means.emplace_back(move(mean_vec));
        
      /* insert new covariance matrix */
      va_vcovs.emplace_back(get_vcov_from_trian(t, rng_dim));
      t += (rng_dim * (rng_dim + 1L)) / 2L;
    }
  }
  
  /* assign constant and fixed effects objects */
  vector<Type> const eta_fix  = X  * b,
                     etaD_fix = XD * b;
  
  /* handle terms from conditional density of observed outcomes */
  Type out(0);
  GVA_cond_dens<Type> const cond_dens(eps, kappa, GH_x, GH_w);
  if(link != "PH")
    error("GVA_MD: link '%s' is not supported", link.c_str());
  
  Type const one(1.);
  for(unsigned i = 0; i < grp.size(); ++i){
    auto const g = grp[i];
    vecT const z = Z.row(i);
    Type const err_mean = (z * va_means[g]).sum(), 
               err_var  = (z * (vecT(va_vcovs[g] * z))).sum(), 
               err_sd   = sqrt(err_var);
    
    out += cond_dens.ph(
      eta_fix[i], etaD_fix[i], one, event[i], err_mean, err_sd, err_var);
  }
  
  /* handle terms from random effect log density and VA log density */
  {
    matrix<Type> va_cov_sum(rng_dim, rng_dim);
    va_cov_sum.setZero();
    Type lb_term(0.);
    for(unsigned g = 0; g < n_groups; ++g){
      lb_term += atomic::logdet(va_vcovs[g]) -
        (va_means[g] * vecT(vcov_inv * va_means[g])).sum();
      va_cov_sum += va_vcovs[g];
      
    }
    lb_term -= (va_cov_sum * vcov_inv).trace();
    
    out += lb_term / Type(2.);
  }
  
  /* add final log determinant term and constants */
  out += Type(n_groups) * (-log_det_vcov + Type(rng_dim)) / Type(2.);
  
  return -out;
}

template<class Type> 
Type GVA(COMMON_ARGS, vector<Type> const &theta_VA, unsigned const n_nodes){
  /* checks */
  unsigned const rng_dim = get_rng_dim(theta);
  {
    int const *max_grp = 
      std::max_element(grp.data(), grp.data() + grp.size());
    /* require mean and full covariance matrix */
    unsigned const expe_size = (*max_grp + 1L) * 
      (rng_dim + (rng_dim * (rng_dim + 1L)) / 2L); 
    if(theta_VA.size() != expe_size)
      error("theta_VA.size(): %i. Expected %i.", 
            theta_VA.size(), expe_size);
  }
  
  auto const GH_xw = GaussHermite::gaussHermiteData(n_nodes);
  vector<Type> GH_x(n_nodes), GH_w(n_nodes);
  for(unsigned i = 0; i < n_nodes; ++i){
    GH_x[i] = Type(GH_xw.x[i]);
    GH_w[i] = Type(GH_xw.w[i]);
  }
  
  if(rng_dim == 1L)
    return GVA_1D(COMMON_CALL, theta_VA, GH_x, GH_w);
  
  return   GVA_MD(COMMON_CALL, theta_VA, GH_x, GH_w, rng_dim);
}


/* 
 Compute psi function w/ adaptive Gauss–Hermite quadrature. That is, 
 
 \begin{align*} 
 f(\sigma^2) &= \int 2 \phi(z; \sigma^2)\Phi(z)\log \Phi(z)dz \\
 &= \frac 2{\sqrt{2\pi\sigma^2}}
 \int \exp\left(-\frac {z^2}{2\sigma^2}\right)\Phi(z)\log \Phi(z)dz \\
 &= \frac 2{\sqrt{2\pi\sigma^2}}
 \int \exp\left(-\frac {z^2}
 {2\gamma^2\sigma^2/(\gamma^2 + \sigma^2)}\right)
 \underbrace{\exp\left(\frac {z^2}{2\gamma^2}\right)
 \Phi(z)\log \Phi(z)}_{f(z;\gamma)}dz \\
 &\overset{z = s \sqrt{2}\gamma\sigma/\sqrt{\gamma^2 + \sigma^2}}{=} 
 \frac 2{\sigma\sqrt{2\pi}}
 \frac {\sqrt 2\sigma\gamma}{\sqrt{\gamma^2 + \sigma^2}}
 \int \exp\left(-s^2\right)
 f(s \sqrt{2}\gamma\sigma/\sqrt{\gamma^2 + \sigma^2};\gamma)ds \\
 &\approx 
 \frac {2\gamma}{\sqrt{\pi(\gamma^2 + \sigma^2)}}
 \sum_{i = 1}^n  
 w_i f(x_i \sqrt{2}\gamma\sigma/\sqrt{\gamma^2 + \sigma^2})
 \end{align*}
 
 with \gamma(\sigma) = 1.
*/
template<class Type>
Type snva_entropy_term
  (Type const sigma_sq, vector<Type> const &GH_x, vector<Type> const &GH_w){
  unsigned const n_nodes = GH_x.size();
  Type const gamma(1.), 
          gamma_sq = gamma * gamma, 
               two(2.);
  
  Type out(0.);
  Type const
    mult_sum(Type(M_2_SQRTPI) * gamma / sqrt(sigma_sq + gamma_sq)), 
        mult(mult_sum * sqrt(sigma_sq) / Type(sqrt(M_2_PI)));
  
  for(unsigned i = 0; i < n_nodes; ++i){
    Type const xi = GH_x[i] * mult;
    out += 
      GH_w[i] * exp(xi * xi / two / gamma_sq) * pnorm(xi) * pnorm_log(xi);
  }
  
  return mult_sum * out;
}

/* Computes the approximate mode of the skew-normal distribution and the 
 * Hessian at the approximation */
template<class Type>
struct SNVA_mode_n_Hess {
  Type mode, Hess;
};

template<class Type> 
SNVA_mode_n_Hess<Type> get_SNVA_mode_n_Hess
  (Type const mu, Type const sigma, Type const rho){
  Type const one(1.), 
             two(2.), 
            zero(0.),
             eps(std::numeric_limits<double>::epsilon()),
           alpha = sigma * rho, 
          a_sign = CppAD::CondExpLe(alpha, zero, -one, one),
              nu = Type(sqrt(M_2_PI)) * alpha / sqrt(1. + alpha * alpha), 
           nu_sq = nu * nu,
           gamma = Type((4. - M_PI) / 2.) * nu_sq * nu / 
             pow(one - nu_sq, Type(3./2.)), 
            mode = mu + sigma * (
              nu - gamma * sqrt(one - nu_sq) / two -
                a_sign / two * exp(- Type(2. * M_PI) / (
                    a_sign * alpha + eps))), 
               z = rho * (mode - mu), 
             phi = dnorm(z, zero, one), 
             Phi = pnorm(z), 
            Hess = - one / sigma / sigma - rho * rho * phi * (
              z * Phi + phi) / (Phi * Phi + eps);
  
  return { mode, Hess };
}

/* Makes an approximation of
 l(\mu,\sigma, \rho) =\frac{2}{\sigma\sqrt{2\pi}}
 \int \exp\left(-\frac{(z - \mu)^2}{2\sigma^2}\right)
 \Phi(\rho(z - \mu))
 \log(1 + \exp(z)) dz
 
 and the mode of the random effect density. The approximation seems to 
 poor if there is a large skew.
 */
template<class Type> 
Type SNVA_mlogit_integral
  (Type const mu, Type const sigma, Type const rho, 
   vector<Type> const &GH_x, vector<Type> const &GH_w){
  /* compute xi and lambda and assign constants */
  auto const dvals = get_SNVA_mode_n_Hess(mu, sigma, rho);
  Type const one(1.), 
             two(2.),
        too_large(30.),
         sigma_sq = sigma * sigma,
               xi = dvals.mode, 
           lambda = one / sqrt(-dvals.Hess);
  
  /* compute approximate integral with Gauss–Hermite quadrature */
  Type out(0.);
  Type const 
    mult_sum(lambda / sigma * Type(M_2_SQRTPI)),
    mult    (Type(M_SQRT2) * lambda);
  
  for(unsigned i = 0; i < GH_x.size(); ++i){
    Type const o = GH_x[i],
              xo = xi + mult * o, 
         xo_diff = (xo - mu),  
             fac = CppAD::CondExpGe(
               xo, too_large, xo, log(one + exp(xo)));
    
    out += 
      GH_w[i] * exp(o * o - xo_diff * xo_diff / two / sigma_sq) * 
      pnorm(rho * (xo - mu)) * fac;
  }
  
  return mult_sum * out;
}
  
/* Makes an approximatin of 
 \begin{align*}
 (\tilde\mu, \tilde \sigma^2, \tilde \rho) &= 
 (s\mu + \log k, s^2\sigma^2,\rho/s)\\ 
 l(\mu,\sigma, \rho, k, s) &=  
 2\int\phi(x;\mu,\sigma^2)
 \Phi(\rho(x - \mu))
 \log(1 + k\exp(sx)) dx \\
 &= 2\int\phi(z;\tilde\mu,\tilde\sigma^2)
 \Phi(\tilde\rho(z - \tilde\mu))
 \log(1 + \exp(z)) dz \\
 &=\frac{2}{\tilde\sigma\sqrt{2\pi}}
 \int \exp\left(-\frac{(z - \tilde\mu)^2}{2\tilde\sigma^2}\right)
 \Phi(\tilde\rho(z - \tilde\mu))
 \log(1 + \exp(z)) dz 
 \end{align*}
 */
template<class Type> 
Type SNVA_mlogit_integral
  (Type const mu, Type const sigma, Type const rho, Type const log_k, 
   Type const s, vector<Type> const &GH_x, vector<Type> const &GH_w){
  return SNVA_mlogit_integral(
    s * mu + log_k, sqrt(s * s) * sigma, rho / s, GH_x, GH_w);
}

/* Makes an approximation of
 l(\mu,\sigma, \rho) =  
 2\int\phi(z;\mu,\sigma^2)
 \Phi(\rho(z - \mu))
 (-\log \Phi(z)) dz
 
 and the mode of the random effect density. The approximation seems to 
 poor if there is a large skew.
 */
template<class Type> 
Type SNVA_probit_integral
  (Type const mu, Type const sigma, Type const rho, 
   vector<Type> const &GH_x, vector<Type> const &GH_w){
  /* compute xi and lambda and assign constants */
  auto const dvals = get_SNVA_mode_n_Hess(mu, sigma, rho);
  Type const one(1.), 
             two(2.),
        sigma_sq = sigma * sigma,
              xi = dvals.mode, 
          lambda = one / sqrt(-dvals.Hess);
  
  /* compute approximate integral with Gauss–Hermite quadrature */
  Type out(0.);
  Type const 
    mult_sum(lambda / sigma * Type(M_2_SQRTPI)),
    mult    (Type(M_SQRT2) * lambda);
  
  for(unsigned i = 0; i < GH_x.size(); ++i){
    Type const o = GH_x[i],
              xo = xi + mult * o, 
         xo_diff = (xo - mu);
    
    out += 
      GH_w[i] * exp(o * o - xo_diff * xo_diff / two / sigma_sq) * 
      pnorm(rho * (xo - mu)) * (-pnorm_log(xo));
  }
  
  return mult_sum * out;
}

/* Makes an approximation of
 \begin{align*}
 (\tilde\mu, \tilde \sigma^2, \tilde \rho) &= 
 (s\mu + k, s^2\sigma^2,\rho/s)\\ 
 l(\mu,\sigma, \rho, k, s) &=  
 2\int\phi(z;\tilde\mu,\tilde\sigma^2)
 \Phi(\tilde\rho(z - \tilde\mu))
 (-\log \Phi(z)) dz \\
 &=\frac{2}{\tilde\sigma\sqrt{2\pi}}
 \int \exp\left(-\frac{(z - \tilde\mu)^2}{2\tilde\sigma^2}\right)
 \Phi(\tilde\rho(z - \tilde\mu))
 (-\log \Phi(z)) dz
 \end{align*}
 */
template<class Type> 
Type SNVA_probit_integral
  (Type const mu, Type const sigma, Type const rho, Type const k, 
   Type const s, vector<Type> const &GH_x, vector<Type> const &GH_w){
  return SNVA_probit_integral(s * mu + k, sqrt(s * s) * sigma, rho / s, 
                              GH_x, GH_w);
}

#define SNVA_COND_DENS_ARGS                                      \
  Type const eta_fix, Type const etaD_fix, Type const z,         \
  Type const event, Type const va_mu, Type const va_sd,          \
  Type const va_rho, Type const va_d, Type const va_var,         \
  Type const dist_mean, Type const dist_var

/* util class to compute the conditional density of the observed outcomes */
template<class Type>
class SNVA_cond_dens {
  /* input dependend objects */
  Type const eps, 
         eps_log = Type(log(eps)), 
           kappa;
  vector<Type> const &GH_x,
                     &GH_w;
  
  /* potentially needed constants */
  Type const mlog_2_pi_half = Type(-log(2 * M_PI) / 2.), 
                        two = Type(2.);
            
public:
  SNVA_cond_dens
  (Type const eps, Type const kappa, vector<Type> const &GH_x, 
   vector<Type> const &GH_w): 
  eps(eps), kappa(kappa), GH_x(GH_x), GH_w(GH_w) { }
  
  /* computes the conditional density term for the PH (log-log) link 
  * function */
  Type ph(SNVA_COND_DENS_ARGS) const {
    Type const h = etaD_fix * exp(eta_fix + z * dist_mean),
        z_scaled = z * va_sd,
               H = two * exp(
                 eta_fix + z * va_mu + 
                   z_scaled * z_scaled / two) * pnorm(z * va_d),
          if_low = event * eps_log - H - h * h * kappa,
          if_ok  = event * log(h)  - H;
    
    return CppAD::CondExpGe(h, eps, if_ok, if_low);
  }
  
  /* computes the conditional density term for the PO (-logit) link 
   * function */
  Type po(SNVA_COND_DENS_ARGS) const {
    Type const H = SNVA_mlogit_integral(
      va_mu, va_sd, va_rho, eta_fix, z, GH_x, GH_w),
               h = etaD_fix * exp(eta_fix + z * dist_mean - H), 
          if_low = event * eps_log - H - h * h * kappa,
          if_ok  = event * log(h)  - H;
    
    return CppAD::CondExpGe(h, eps, if_ok, if_low);
  }
  
  /* computes the conditional density term for the probit (-probit) link 
   * function */
  Type probit(SNVA_COND_DENS_ARGS) const {
    Type const H = SNVA_probit_integral(
        va_mu, va_sd, va_rho, -eta_fix, -z, GH_x, GH_w),
            diff = (eta_fix + z * dist_mean),
               h = etaD_fix * exp(
                 mlog_2_pi_half - diff * diff / two - 
                   z * z * dist_var / two + H), 
            if_low = event * eps_log - H - h * h * kappa,
            if_ok  = event * log(h)  - H;
    
    return CppAD::CondExpGe(h, eps, if_ok, if_low);
  }
};

/* The following functions maps from the input vector to the 
 * parameterization used in 
 * > Ormerod, J. T. (2011). Skew-normal variational approximations for 
 * > Bayesian inference. Unpublished article. */
template<class Type>
struct SNVA_1D_input {
  vector<Type> va_mus, va_sds_log, va_rhos;
};

/* maps from direct parameters vector to direct parameters */
template<class Type>
SNVA_1D_input<Type> SNVA_1D_theta_DP_to_DP(vector<Type> const &theta_VA){
  auto get_va_param_vecs = [&](unsigned const shift){
    unsigned const n_elem = theta_VA.size() / 3L;
    vector<Type> out(n_elem);
    for(unsigned i = 0; i < n_elem; ++i)
      out[i] = theta_VA[i * 3L + shift];
    
    return out;
  };
  
  return { get_va_param_vecs(0L), get_va_param_vecs(1L), 
           get_va_param_vecs(2L) };
}

/* maps from mean, log standard deviation, and Pearson's moment coefficient 
 * of skewness to direct parameters. The first function overrides the input
 * vectors values and moves them at the end */
template<class Type>
SNVA_1D_input<Type> SNVA_1D_CP_to_DP
  (vector<Type> &mus, vector<Type> &sigma_logs, vector<Type> &gammas){
  using std::move;
  
  /* needed constants */
  Type const one(1.),
             two(2.), 
             Tpi(M_PI);
  
  unsigned const n_elem = mus.size();
  for(unsigned i = 0; i < n_elem; ++i){
    Type const mu = mus[i], 
         sigma_sq = exp(2 * sigma_logs[i]), 
            gamma = gammas[i];
    
    /* compute needed intermediaries */
    Type const mu_z = gamma_to_nu(gamma), 
            mu_z_sq = mu_z * mu_z,
           omega_sq = sigma_sq / (one - mu_z_sq);
    
    /* insert values */
    mus[i]        =  mu - sqrt(omega_sq) * mu_z;
    sigma_logs[i] = log(omega_sq) / two;
    gammas[i]     = sqrt(Tpi / omega_sq) * mu_z / sqrt(two - Tpi * mu_z_sq);
    
  }
  return { move(mus), move(sigma_logs), move(gammas) };
}

template<class Type>
SNVA_1D_input<Type> SNVA_1D_theta_CP_to_DP
  (vector<Type> const &theta_VA){
  unsigned const n_elem = theta_VA.size() / 3L;
  vector<Type>     mus(n_elem), 
            sigma_logs(n_elem), 
                gammas(n_elem);
  
  /* extract parameters */
  for(unsigned i = 0; i < n_elem; ++i){
    unsigned const grp = i * 3L;
    mus[i]        = theta_VA[grp     ], 
    sigma_logs[i] = theta_VA[grp + 1L],
    gammas[i]     = theta_VA[grp + 2L];
  }
  
  return SNVA_1D_CP_to_DP(mus, sigma_logs, gammas);
}

template<class Type>
class get_gamma{
  Type const c1 = Type(0.99527), 
             c2 = Type(2.) * c1, 
            one = Type(1.);
public:
  /* inverse of 2 * c1 * logit(gamma) - c1 */
  Type operator()(Type const gtrans) const {
    return c2 / (one + exp(-gtrans)) - c1;
  }
};

/* maps from mean, log standard deviation, and __transformed__  
 * Pearson's moment coefficient of skewness to direct parameters */
template<class Type>
SNVA_1D_input<Type> SNVA_1D_theta_CP_trans_to_DP
  (vector<Type> const &theta_VA){
  using std::move;
  unsigned const n_elem = theta_VA.size() / 3L;
  vector<Type> mus(n_elem), 
        sigma_logs(n_elem), 
            gammas(n_elem);
  
  /* extract parameters */
  get_gamma<Type> trans_g;
  for(unsigned i = 0; i < n_elem; ++i){
    unsigned const grp = i * 3L;
    mus[i]        =         theta_VA[grp     ];
    sigma_logs[i] =         theta_VA[grp + 1L];
    gammas[i]     = trans_g(theta_VA[grp + 2L]);
    
  }
  
  return SNVA_1D_CP_to_DP(mus, sigma_logs, gammas);
}

/* TODO: describe SNVA implementation */
template<class Type> 
Type SNVA_1D
  (COMMON_ARGS, vector<Type> const &theta_VA, vector<Type> const &GH_x, 
   vector<Type> const &GH_w, std::string const &param_type){
  using std::move;
  
  /* needed constants */
  Type const half(.5), 
              one(1.),
              two(2.);
  
  /* assign vectors w/ variational parameters */
#define SET_PARAMS(meth_use)                                   \
  auto params = meth_use(theta_VA);                            \
  va_mus      = move(params.va_mus);                           \
  va_sds_log  = move(params.va_sds_log);                       \
  va_rhos     = move(params.va_rhos)                            

  vector<Type> va_mus, va_sds_log, va_rhos;
  if(param_type == "DP"){
    SET_PARAMS(SNVA_1D_theta_DP_to_DP);
  } else if (param_type == "CP_trans"){
    SET_PARAMS(SNVA_1D_theta_CP_trans_to_DP);
  } else if(param_type == "CP"){
    SET_PARAMS(SNVA_1D_theta_CP_to_DP);
  } else
    error("SNVA_1D: param_type '%s' is not implemented", 
          param_type.c_str());
#undef SET_PARAMS
  
  vector<Type> const va_sds  = exp(va_sds_log),
                     va_vars = va_sds * va_sds;
  
  /* get random effect variance, fixed effect linear predictors, and 
   * assign constants */
  Type const sigma_log = theta[0], 
             sigma     = exp(sigma_log),
             T_sqrt_2_pi(sqrt(M_2_PI));
  vector<Type> const eta_fix  = X  * b,
                     etaD_fix = XD * b;
  
  auto const va_d = ([&](){
    unsigned const n_elem = va_rhos.size(); 
    vector<Type> out(n_elem);
    for(unsigned g = 0; g < n_elem; ++g){
      Type const rho_var = va_rhos[g] * va_vars[g];
      out[g] = rho_var / sqrt(one + rho_var * va_rhos[g]);
    }
    
    return out;
  })();
  
  /* TODO: we may already have these value with CP parameterization */
  vector<Type> const dist_means = va_mus + T_sqrt_2_pi * va_d, 
                      dist_vars = va_vars - Type(M_2_PI) * va_d * va_d;
  
  /* handle terms from conditional density of observed outcomes */
  Type out(0);
  SNVA_cond_dens<Type> const cond_dens(eps, kappa, GH_x, GH_w);
  
#define ADD_COND_DENS(func)                                    \
  for(unsigned i = 0; i < grp.size(); ++i){                    \
    auto const g = grp[i];                                     \
    out += cond_dens.func(                                     \
      eta_fix[i], etaD_fix[i], Z(i, 0), event[i], va_mus[g],   \
      va_sds[g], va_rhos[g], va_d[g], va_vars[g],              \
      dist_means[g], dist_vars[g]);                            \
  }
  
  if(link == "PH")
    ADD_COND_DENS(ph)
  else if (link == "PO")
    ADD_COND_DENS(po)
  else if (link == "probit")
    ADD_COND_DENS(probit)
  else 
    error("'%s' not implemented", link.c_str());
#undef ADD_COND_DENS
  
  /* handle terms from random effect log density and VA log density */
  {
    Type mu_sq_sum(0), mu_d_sum(0), var_sum(0.);
    for(unsigned g = 0; g < va_mus.size(); g++){
      mu_sq_sum += va_mus [g] * va_mus[g];
      mu_d_sum  += va_mus [g] * va_d  [g];
      var_sum   += va_vars[g];
      
      /* add non-constant entropy terms */
      out       += va_sds_log[g] -
        snva_entropy_term(va_rhos[g] * va_rhos[g] * va_vars[g], 
                          GH_x, GH_w);
    }
    
    /* terms from expected log normal density */
    Type const n_groups(va_mus.size());
    out += 
      - n_groups * sigma_log
      - ((mu_sq_sum + var_sum) / two + T_sqrt_2_pi * mu_d_sum) / 
        sigma / sigma;
    
    /* add constant terms from entropy */
    out += n_groups * (half - Type(log(2)));
  }
  
  return -out;
}

/* The following functions maps from the input vector to the 
 * parameterization used in 
 * > Ormerod, J. T. (2011). Skew-normal variational approximations for 
 * > Bayesian inference. Unpublished article. */
template<class Type>
struct SNVA_MD_input {
  std::vector<vector<Type> > va_mus, 
                            va_rhos;
  std::vector<matrix<Type> > va_lambdas;
};

/* maps from direct parameters vector to direct parameters though from C 
 * matrix where Lambda = C^\top C */
template<class Type>
SNVA_MD_input<Type> SNVA_MD_theta_DP_to_DP
  (vector<Type> const &theta_VA, unsigned const rng_dim){
  using vecT = vector<Type>;
  using std::move;
  unsigned const n_groups = theta_VA.size() / (
    rng_dim * 2L + (rng_dim * (rng_dim + 1L)) / 2L);
  
  SNVA_MD_input<Type> out;
  std::vector<vecT > &va_mus = out.va_mus, 
                     &va_rhos = out.va_rhos;
  std::vector<matrix<Type> > &va_lambdas = out.va_lambdas;
  
  va_mus    .reserve(n_groups);
  va_rhos   .reserve(n_groups);
  va_lambdas.reserve(n_groups);
  
  Type const *t = &theta_VA[0];
  for(unsigned g = 0; g < n_groups; ++g){
    /* insert new mu vector */
    vecT mu_vec(rng_dim);
    for(unsigned i = 0; i < rng_dim; ++i)
      mu_vec[i] = *t++;
    va_mus.emplace_back(move(mu_vec));
    
    /* insert new lambda matrix */
    va_lambdas.emplace_back(get_vcov_from_trian(t, rng_dim));
    t += (rng_dim * (rng_dim + 1L)) / 2L;
    
    /* insert new rho vector */
    vecT rho_vec(rng_dim);
    for(unsigned i = 0; i < rng_dim; ++i)
      rho_vec[i] = *t++;
    va_rhos.emplace_back(move(rho_vec));
  }
  
  return out;
}

/* maps from mean, C matrix with Cov = C^\top C, and __transformed__  
 * Pearson's moment coefficient of skewness to direct parameters */
template<class Type>
SNVA_MD_input<Type> SNVA_MD_theta_CP_trans_to_DP
  (vector<Type> const &theta_VA, unsigned const rng_dim){
  using vecT = vector<Type>;
  using std::move;
  unsigned const n_mu = rng_dim, 
             n_lambda = (rng_dim * (rng_dim + 1L)) / 2L,
              n_per_g = n_mu + n_lambda + rng_dim,
             n_groups = theta_VA.size() / n_per_g;
  
  SNVA_MD_input<Type> out;
  std::vector<vecT > &va_mus = out.va_mus, 
                     &va_rhos = out.va_rhos;
  std::vector<matrix<Type> > &va_lambdas = out.va_lambdas;
  
  va_mus    .reserve(n_groups);
  va_rhos   .reserve(n_groups);
  va_lambdas.reserve(n_groups);
  
  get_gamma<Type> trans_g;
  Type const *t = &theta_VA[0], 
            one(1.), 
            two(2.), 
            Tpi(M_PI), 
        sqrt_pi(sqrt(M_PI));
  
  /* intermediaries */
  vecT gamma(rng_dim), 
          nu(rng_dim), 
       omega(rng_dim);
  
  for(unsigned g = 0; g < n_groups; ++g, t += n_per_g){
    /* get gamma parameters */
    Type const *gi = t + n_mu + n_lambda;
    for(unsigned i = 0; i < rng_dim; ++i)
      gamma[i] = trans_g(*gi++);
    
    /* Compute intermediaries and rho */
    vecT rhos(rng_dim);
    matrix<Type> Sigma = get_vcov_from_trian(t + n_mu, rng_dim);
    for(unsigned i = 0; i < rng_dim; ++i){
      Type const &gv = gamma[i];
      nu[i]    = gamma_to_nu(gv);
      omega[i] = sqrt(Sigma(i, i) / (one - nu[i] * nu[i]));
      rhos[i]  = 
        sqrt_pi * nu[i] / omega[i] / sqrt(two - Tpi * nu[i] * nu[i]);
      
      /* replace nu by nu * omega */
      nu[i] *= omega[i];
    }
    va_rhos.emplace_back(move(rhos));
    
    /* assign mu and Lambda */
    vecT mu(rng_dim);
    Type const *mi = t;
    for(unsigned i = 0; i < rng_dim; ++i)
      mu[i] = *mi++;
    mu -= nu;
    va_mus.emplace_back(move(mu));
    
    auto const nu_vec = nu.matrix();
    Sigma += nu_vec * nu_vec.transpose();
    va_lambdas.emplace_back(move(Sigma));
  }
  
  return out;
}

template<class Type> 
Type SNVA_MD
  (COMMON_ARGS, vector<Type> const &theta_VA, vector<Type> const &GH_x, 
   vector<Type> const &GH_w, std::string const &param_type, 
   unsigned const rng_dim){
  using vecT = vector<Type>;
  using std::move;
  
  /* get objects related to model covariance matrix */
  matrix<Type> const vcov = get_vcov_from_trian(theta);
  
  /* maybe not the best idea to matrix multiply by the precision matrix 
  * instead of using solve... */
  Type log_det_vcov;
  matrix<Type> vcov_inv;
  vcov_inv = atomic::matinvpd(vcov, log_det_vcov);
  
  /* get objects from VA distribution */
  std::vector<vecT > va_mus, 
                    va_rhos;
  std::vector<matrix<Type> > va_lambdas;
  
#define SET_PARAMS(meth_use)                                   \
  auto const input = meth_use(theta_VA, rng_dim);              \
  va_mus     = move(input.va_mus);                             \
  va_rhos    = move(input.va_rhos);                            \
  va_lambdas = move(input.va_lambdas)                          
  
  if(param_type == "DP"){
    SET_PARAMS(SNVA_MD_theta_DP_to_DP);
  } else if(param_type == "CP_trans"){
    SET_PARAMS(SNVA_MD_theta_CP_trans_to_DP);
  } else
    error("SNVA_MD: param_type '%s' is not implemented", 
          param_type.c_str());
#undef SET_PARAMS
  
  unsigned const n_groups = va_mus.size();
  
  /* assign constant and fixed effect objects */
  vecT const eta_fix = X  * b,
            etaD_fix = XD * b;
  Type const sqrt_2_pi(sqrt(M_2_PI)), 
                  one(1.),
                small(std::numeric_limits<double>::epsilon());
  
  /* assign object used in VA distribution */
  std::vector<vecT> va_ds;
  va_ds.reserve(n_groups);
  for(unsigned g = 0; g < n_groups; ++g){
    vecT new_d = va_lambdas[g] * va_rhos[g];
    Type const denom = sqrt(one + (va_rhos[g] * new_d).sum());
    new_d /= denom;
    va_ds.emplace_back(move(new_d));
  }
  
  /* handle terms from conditional density of observed outcomes */
  Type out(0);
  SNVA_cond_dens<Type> const cond_dens(eps, kappa, GH_x, GH_w);
  if(link != "PH")
    error("SNVA_MD: link '%s' is not supported", link.c_str());
  
  for(unsigned i = 0; i < grp.size(); ++i){
    auto const g = grp[i];
    vecT const z = Z.row(i);
    
    Type const mu = (z * va_mus     [g]).sum(), 
            sd_sq = (z * (va_lambdas[g] * z)).sum(),
               sd = sqrt(sd_sq),
                d = (z * va_ds      [g]).sum(), 
              rho = d / sd_sq / sqrt(one - d * d / sd_sq),
         d_scaled = sqrt_2_pi * d, 
        dist_mean = mu + d_scaled, 
         dist_var = sd_sq - d_scaled * d_scaled;
    
    out += cond_dens.ph(
      eta_fix[i], etaD_fix[i], one, event[i],
      mu, sd, rho, d, sd_sq, dist_mean, dist_var);
  }
  
  /* handle terms from random effect log density and VA log density */
  {
    matrix<Type> va_lambda_sum(rng_dim, rng_dim);
    va_lambda_sum.setZero();
    Type lb_t_mult_half(0.), lb_t_mult_other(0.);
    for(unsigned g = 0; g < n_groups; ++g){
      lb_t_mult_half += atomic::logdet(va_lambdas[g]) -
        (va_mus[g] * (vcov_inv * va_mus[g])).sum();
      va_lambda_sum += va_lambdas[g];
      lb_t_mult_other -= (va_mus[g] * (vcov_inv * va_ds[g])).sum();
      
      auto const llt_mat = va_lambdas[g].llt();
      vecT va_rho_scaled = 
        (llt_mat.matrixU() * va_rhos[g].matrix()).array() + small;
      Type const r_L_r = (va_rho_scaled * va_rho_scaled).sum();
      out -= snva_entropy_term(r_L_r, GH_x, GH_w);
      
    }
    lb_t_mult_half -= (va_lambda_sum * vcov_inv).trace();
    
    out += lb_t_mult_half / Type(2.) + lb_t_mult_other * sqrt_2_pi;
  }
  
  /* add final log determinant term and constants */
  out += Type(n_groups) * (
    -log_det_vcov /  Type(2.) + Type(rng_dim) /  Type(2.) - Type(M_LN2));
  
  return -out;
}

template<class Type> 
Type SNVA(COMMON_ARGS, vector<Type> const &theta_VA, 
          unsigned const n_nodes, std::string const &param_type){
  /* checks */
  unsigned const rng_dim = get_rng_dim(theta);
  {
    int const *max_grp = 
      std::max_element(grp.data(), grp.data() + grp.size());
    /* require mean and full covariance matrix */
    unsigned const expe_size = (*max_grp + 1L) * 
    (rng_dim * 2L + (rng_dim * (rng_dim + 1L)) / 2L); 
    if(theta_VA.size() != expe_size)
      error("theta_VA.size(): %i. Expected %i.", 
            theta_VA.size(), expe_size);
  }
  
  auto const GH_xw = GaussHermite::gaussHermiteData(n_nodes);
  vector<Type> GH_x(n_nodes), GH_w(n_nodes);
  for(unsigned i = 0; i < n_nodes; ++i){
    GH_x[i] = Type(GH_xw.x[i]);
    GH_w[i] = Type(GH_xw.w[i]);
  }
  
  if(rng_dim == 1L)
    return SNVA_1D(COMMON_CALL, theta_VA, GH_x, GH_w, param_type);
  
  return   SNVA_MD(COMMON_CALL, theta_VA, GH_x, GH_w, param_type, rng_dim);
}

/* function to test for approximate equality */
template<typename T> 
bool is_all_equal(T const test, T const expect, T const eps){
  return std::abs(expect) > eps ? 
    std::abs(test - expect) / std::abs(expect) < eps : 
    std::abs(test - expect)                    < eps;
}

/* small class to make testing easier */
struct test_res {
  int code; 
  double test, expect;
  std::string const testname = "";
  
  bool operator()(double const eps) const {
    return !is_all_equal(test, expect, eps);
  }
  
  void call_error() const {
    error("'%s' failed with C: %d, R: %.10f, E: %.10f", 
          testname.c_str(), code, test, expect);
  }
};

#ifdef TEST_PNORM_LOG1
/* test function for pnorm_log. It has to be executed on another thread due
 * to the taping I think */
test_res test_pnorm_log1(){
  using std::abs;
  using std::vector;
  /*
   x <- 2; mu <- 1; sd <- 2
   dput(pnorm(x, mu, sd, log.p = TRUE))
   numericDeriv(quote(pnorm(x, mu, sd, log.p = TRUE)), theta = "x")
   dput(dnorm(x, mu, sd) / pnorm(x, mu, sd))
   */
  constexpr double eps = 1e-8;
  vector<AD<double> > x(1);
  x[0] = 3;
  Independent(x);
  AD<double> const mu = 1, sd = 2;
  vector<AD<double> > y(1);
  y[0] = pnorm_log(x[0], mu, sd);
  CppAD::ADFun<double> func(x, y);
  
  vector<double> xx(1);
  xx[0] = 2;
  vector<double> yy = func.Forward(0, xx);
  test_res t1 { 1L, yy[0], -0.368946415288657, "pnorm_log1" };
  if(t1(eps))
    return t1;
  
  vector<double> w(1);
  w[0] = 1;
  vector<double> dy = func.Reverse(1, w);
  
  test_res t2 { 2L, dy[0], 0.254580216918517, "pnorm_log1" };
  if(t2(eps))
    return t2;
  
  return { 0L, 0., 0. };
}
#endif


/* tests to check cpp functions */
template<class Type>
Type cpp_tests(vector<Type> const &x){
  using std::vector;
  using std::abs;
  {
    /* 
     xw <- fastGHQuad::gaussHermiteData(5)
     paste("{", 
     paste0(sprintf("%20.16f", xw$x, xw$x), collapse = ", "),
     "}")
     paste("{", 
     paste0(sprintf("%20.16f", xw$w, xw$x), collapse = ", "),
     "}")
     */
    constexpr double eps = 1e-10;
    constexpr unsigned n = 5L;
    auto xw = GaussHermite::gaussHermiteData(n);
    std::vector<double> const x = 
      {  -2.0201828704560856,  -0.9585724646138187,   0.0000000000000000,   0.9585724646138187,   2.0201828704560856 };
    std::vector<double> const w = 
      {   0.0199532420590459,   0.3936193231522411,   0.9453087204829413,   0.3936193231522414,   0.0199532420590459 };
    
    for(unsigned i = 0; i < n; ++i){
      if(!is_all_equal(xw.x[i], x[i], eps) or 
           !is_all_equal(xw.w[i], w[i], eps))
        error("incorrect value from 'gaussHermiteData'");
    }
    {
#ifdef TEST_PNORM_LOG1
      /* works the first four times and then it fails. I do not know why... */
      auto fu = std::async(std::launch::async, test_pnorm_log1);
      auto const res = fu.get();
      if(res.code != 0)
        res.call_error();
#endif
    }
  }
  {
    /*
     psi <- function(sigma)
     ifelse(
     sigma < .Machine$double.eps^(1/2), 
     -log(2), 
     sapply(sigma, function(s)
     integrate(function(x) 
     2 * dnorm(x, sd = s) * pnorm(x) * pnorm(x, log.p = TRUE),
     -Inf, Inf, rel.tol = .Machine$double.eps^(3/4))$value))
     dput(psi(0:4))
     */
    vector<double> sigmas = { 0., 1., 2., 100. };
    vector<double> expect = {
      -0.693147180559945, -0.5, -0.319885055354923, -0.00720608477031823 };
    constexpr unsigned const n_nodes = 20L;
    auto xw = GaussHermite::gaussHermiteData(n_nodes);
    ::vector<double> x(n_nodes), w(n_nodes);
    for(unsigned i = 0; i < n_nodes; ++i){
      x[i] = xw.x[i];
      w[i] = xw.w[i];
    }
    
    for(unsigned i = 0L; i < 4L; ++i){
      double res = snva_entropy_term(sigmas[i] * sigmas[i], x, w);
      test_res test_result { (int)i, res, expect[i], "snva_entropy_term" };
      if(test_result(1e-8))
        test_result.call_error();
    }
    
    double z = 1e-4;
    for(unsigned i = 0; i < 9L; ++i, z /= 10.){
      double res = snva_entropy_term(z * z, x, w);
      test_res test_result { (int)i, res, -log(2), 
                             "snva_entropy_term (small)" };
      if(test_result(1e-8))
        test_result.call_error();
    }
    {
      double res = snva_entropy_term(0., x, w);
      test_res test_result { 0L, res, -log(2), 
                             "snva_entropy_term (zero)" };
      if(test_result(1e-8))
        test_result.call_error();
    }
  }
  {
    /*
     integrand <- function(x, mu, sigma, k, s, use_log = FALSE){
     f1  <- dnorm(x, mean = mu, sd = sigma, log = use_log) 
     eta <- log(k) + s * x
     f2  <- ifelse(eta > 30, eta, log(1 + exp(eta)))
     f1 * f2
     }
     library(fastGHQuad)
     psi <- function(mu, sigma, k, s, rule){
     x      <- mu + sqrt(2) * sigma * rule$x
     f <- integrand(x, mu, sigma, k, s) * exp((x - mu)^2 / 2 / sigma^2)
     sigma * sqrt(2) * drop(rule$w %*% f) 
     }
     wx <- gaussHermiteData(15L)
     dput(mapply(psi, 
     mu    = c(-1,   1,  0,  2), 
     sigma = c(.01,  1, 10,  5), 
     k     = c(  2, .5,  5,  6), 
     s     = c(  1,  2, -1, -4), 
     MoreArgs = list(rule = wx)))
     */
    unsigned const n_nodes = 15L;
    auto xw = GaussHermite::gaussHermiteData(n_nodes);
    ::vector<double> x(n_nodes), w(n_nodes);
    for(unsigned i = 0; i < n_nodes; ++i){
      x[i] = xw.x[i];
      w[i] = xw.w[i];
    }
    vector<double> const mu = { -1,   1,  0,  2 }, 
                      sigma = { .01,  1, 10,  5 }, 
                          k = {   2, .5,  5,  6 }, 
                          s = {   1,  2, -1, -4 }, 
                     ex_res = { 0.551456924101031, 1.84776400976526, 4.88998305310309, 5.42588929528414 };
    for(unsigned i = 0; i < mu.size(); ++i){
      double const intval = GVA_mlogit_integral(
        mu[i], sigma[i], log(k[i]), s[i], x, w);
      test_res res { (int)i, intval, ex_res[i], "GVA_mlogit_integral" };
      if(res(1e-8))
        res.call_error();
    }
  }
  {
    /*
     integrand <- function(x, mu, sigma, rho, k, s, use_log = FALSE){
     f1 <- dnorm(x, mean = mu, sd = sigma, log = use_log) 
     f2 <- -pnorm(k + s * x, log.p = TRUE)
     if(use_log) f1 + log(f2) else f1 * f2
     }
     psi <- function(mu, sigma, rho, k, s)
     integrate(
     integrand, mu = mu, sigma = sigma, rho = rho, k = k, s = s, 
     lower = -Inf, upper = Inf, rel.tol = 1e-10, abs.tol = 0)$value
     # the approximation should work well for these values
     dput(mapply(
     psi, 
     mu    = c(0, -1, 1, 0, -1), 
     sigma = c(1, .5, 1, 2, .7), 
     k     = c(1, -1, 2, 1, .8), 
     s     = c(1, .2, 1, 1, -1)))
     */
    unsigned const n_nodes = 40L;
    auto xw = GaussHermite::gaussHermiteData(n_nodes);
    ::vector<double> x(n_nodes), w(n_nodes);
    for(unsigned i = 0; i < n_nodes; ++i){
      x[i] = xw.x[i];
      w[i] = xw.w[i];
    }
    vector<double> const 
      mu     = { 0, -1, 1, 0, -1 }, 
      sigma  = { 1, .5, 1, 2, .7 },
      k      = { 1, -1, 2, 1, .8 }, 
      s      = { 1, .2, 1, 1, -1 }, 
      ex_res = { 0.35934760083045, 2.16633048203693, 0.0187669055995364, 0.942178782481152, 
                 0.078745703164216 };
    for(unsigned i = 0; i < mu.size(); ++i){
      double const intval = GVA_probit_integral(
        mu[i], sigma[i], k[i], s[i], x, w);
      test_res res { (int)i, intval, ex_res[i], "GVA_probit_integral" };
      if(res(1e-8))
        res.call_error();
    }
  }
  {
    /*
     integrand <- function(x, mu, sigma, rho, k, s, use_log = FALSE){
     f1 <- dnorm(x, mean = mu, sd = sigma, log = use_log) 
     f2 <- 2 * pnorm(rho * (x - mu), log.p = use_log)
     eta <- log(k) + s * x
     f3  <- ifelse(eta > 30, eta, log(1 + exp(eta)))
     if(use_log) f1 + f2 + log(f3) else f1 * f2 * f3
     }
     psi <- function(mu, sigma, rho, k, s)
     integrate(
     integrand, mu = mu, sigma = sigma, rho = rho, k = k, s = s, 
     lower = -Inf, upper = Inf, rel.tol = 1e-10, abs.tol = 0)$value
     # the approximation should work well for these values
     dput(mapply(
     psi, 
     mu    = c(0, -1, 1, 0, -1), 
     sigma = c(1, .5, 1, 2, .7), 
     rho   = c(0, -1, 1, 1, -2), 
     k     = c(1, .3, 2, 1, .8), 
     s     = c(1, .2, 1, 1, .4)))
     */
    unsigned const n_nodes = 40L;
    auto xw = GaussHermite::gaussHermiteData(n_nodes);
    ::vector<double> x(n_nodes), w(n_nodes);
    for(unsigned i = 0; i < n_nodes; ++i){
      x[i] = xw.x[i];
      w[i] = xw.w[i];
    }
    vector<double> const 
        mu     = { 0, -1, 1, 0, -1 }, 
        sigma  = { 1, .5, 1, 2, .7 },
        rho    = { 0, -1, 1, 1, -2 },
        k      = { 1, .3, 2, 1, .8 }, 
        s      = { 1, .2, 1, 1, .4 }, 
        ex_res = { 0.80605918334744, 0.213372360274554, 2.38729257810521, 1.78136403451249, 
                   0.374363761000257 };
    for(unsigned i = 0; i < mu.size(); ++i){
      double const intval = SNVA_mlogit_integral(
        mu[i], sigma[i], rho[i], log(k[i]), s[i], x, w);
      test_res res { (int)i, intval, ex_res[i], "SNVA_mlogit_integral" };
      if(res(1e-8))
        res.call_error();
    }
  }
  {
    /*
     integrand <- function(x, mu, sigma, rho, k, s, use_log = FALSE){
     f1 <- dnorm(x, mean = mu, sd = sigma, log = use_log) 
     f2 <- pnorm(rho * (x - mu), log.p = use_log)
     f3 <- -2 * pnorm(k + s * x, log.p = TRUE)
     f1 * f2 * f3
     }
     psi <- function(mu, sigma, rho, k, s)
     integrate(
     integrand, mu = mu, sigma = sigma, rho = rho, k = k, s = s, 
     lower = -Inf, upper = Inf, rel.tol = 1e-10, abs.tol = 0)$value
     # the approximation should work well for these values
     dput(mapply(
     psi, 
     mu    = c(0, -1,  1,  0,  -1), 
     sigma = c(1, .5,  1,  2,  .7), 
     rho   = c(0, -1,  1,  1,  -2), 
     k     = c(1, .3, -2,  1, -.8), 
     s     = c(1, .2,  1, -1,  .4)))
     */
    unsigned const n_nodes = 40L;
    auto xw = GaussHermite::gaussHermiteData(n_nodes);
    ::vector<double> x(n_nodes), w(n_nodes);
    for(unsigned i = 0; i < n_nodes; ++i){
      x[i] = xw.x[i];
      w[i] = xw.w[i];
    }
    vector<double> const 
        mu     = { 0, -1,  1,  0,  -1 }, 
        sigma  = { 1, .5,  1,  2,  .7 },
        rho    = { 0, -1,  1,  1,  -2 },
        k      = { 1, .3, -2,  1, -.8 }, 
        s      = { 1, .2,  1, -1,  .4 }, 
        ex_res = { 0.35934760083045, 0.645851484652832, 1.34185746494562, 1.79323321188706, 
                   2.50177463219666 };
    for(unsigned i = 0; i < mu.size(); ++i){
      double const intval = SNVA_probit_integral(
        mu[i], sigma[i], rho[i], k[i], s[i], x, w);
      test_res res { (int)i, intval, ex_res[i], "SNVA_probit_integral" };
      if(res(1e-8))
        res.call_error();
    }
  }
  {
    /* 
     dput(mu <- c(4, 3, 8))
     dput(Sig <- matrix(c(1, 1, 2, 1, 4, 3, 2, 3, 6), 3))
     dput(rho <- c(5, 7, 9))
     dput(c(mu, chol(Sig)[upper.tri(Sig, TRUE)], rho))
     */
    ::vector<Type> theta(12);
    theta << 4., 3., 8., 1., 1., 1.73205080756888, 2., 0.577350269189626, 1.29099444873581, 
             5., 7., 9.;
    
    auto input = SNVA_MD_theta_DP_to_DP(theta, 3L);
    vector<double> const mu = { 4, 3, 8 }, 
                        Sig = { 1, 1, 2, 1, 4, 3, 2, 3, 6 }, 
                        rho = { 5, 7, 9 };
    {
      test_res res { 0L, (double)input.va_mus.size(), 1., 
                     "SNVA_MD_theta_DP_to_DP" };
      if(res(1e-8))
        res.call_error();
    }
    
#define TEST_ELE(res, ex)                                    \
{                                                            \
  test_res r1 { 0L, (double)input.res[0].size(),             \
                (double)ex.size(),                           \
                "SNVA_MD_theta_DP_to_DP (" #ex ")" };        \
  if(r1(1e-8))                                               \
    r1.call_error();                                         \
  for(unsigned i = 0; i < ex.size(); ++i){                   \
    test_res r2 { (int)(i + 1L),                             \
                   asDouble(*(input.res[0].data() + i)),     \
                   ex[i],                                    \
                   "SNVA_MD_theta_DP_to_DP (" #ex ")" };     \
    if(r2(1e-8))                                             \
      r2.call_error();                                       \
  }                                                          \
}
  
    TEST_ELE(va_mus, mu)
    TEST_ELE(va_lambdas, Sig)
    TEST_ELE(va_rhos, rho)
#undef TEST_ELE
  }
  {
    /*
     cp_to_dp <- function(mu, Sigu, gamma){
     Sig <- matrix(0., length(mu), length(mu))
     Sig[upper.tri(Sig, TRUE)] <- Sigu
     Sig <- crossprod(Sig)
     
     gamma <- 2 * 0.99527 / (1 + exp(-gamma)) - 0.99527
     cv <- 2 * abs(gamma) / (4 - pi)
     nu <- ifelse(gamma < 0, -1, 1) * cv^(1/3) / sqrt(1 + cv^(2/3))
     omegas <- sqrt(diag(Sig) / (1 - nu^2))
     rhos <- sqrt(pi) * nu / sqrt(2 - pi * nu * nu) / omegas
     onu <- omegas * nu
     Lambda <- Sig + onu %o% onu
     
     list(xi = mu - onu, Lambda = Lambda, rhos = rhos)
     }
     dput(cp_to_dp(
     mu = c(-1, 0, 1), 
     Sigu = c(4, 2, 1, 1, 2, 3), 
     gamma = c(-2, .5, 0)))
     */
    ::vector<Type> theta(12);
    theta << -1., 0., 1.,
             4., 2., 1., 1., 2., 3., 
             -2., .5, 0.;
    
    auto input = SNVA_MD_theta_CP_trans_to_DP(theta, 3L);
    vector<double> const mu = { 3.83496892053337, -1.85176038792831, 1 }, 
      Sig = { 39.3769244625236, -0.953203923908205, 4, -0.953203923908205, 
              8.42901653430041, 4, 4, 4, 14 }, 
      rho = { -0.592481457389101, 0.458273315832141, 0 };
    {
        test_res res { 0L, (double)input.va_mus.size(), 1., 
                       "SNVA_MD_theta_CP_trans_to_DP" };
        if(res(1e-8))
          res.call_error();
    }
    
#define TEST_ELE(res, ex)                                     \
{                                                             \
  test_res r1 { 0L, (double)input.res[0].size(),              \
                (double)ex.size(),                            \
                "SNVA_MD_theta_CP_trans_to_DP (" #ex ")" };   \
  if(r1(1e-8))                                                \
    r1.call_error();                                          \
  for(unsigned i = 0; i < ex.size(); ++i){                    \
    test_res r2 { (int)(i + 1L),                              \
                  asDouble(*(input.res[0].data() + i)),       \
                  ex[i],                                      \
                  "SNVA_MD_theta_CP_trans_to_DP (" #ex ")" }; \
    if(r2(1e-8))                                              \
      r2.call_error();                                        \
  }                                                           \
}

  TEST_ELE(va_mus, mu)
  TEST_ELE(va_lambdas, Sig)
  TEST_ELE(va_rhos, rho)
#undef TEST_ELE
  }
  
  return x[0];
}
} // namespace

template<class Type>
Type objective_function<Type>::operator() ()
{
  /* common data objects and parameters */
  DATA_STRING(app_type);
  DATA_VECTOR(tobs);
  
  /* maybe run tests at this point */
  if(app_type == "run tests")
    return cpp_tests(tobs);
  
  DATA_VECTOR(event);
  DATA_MATRIX(X);
  DATA_MATRIX(XD);
  DATA_MATRIX(Z);
  DATA_IVECTOR(grp);
  DATA_STRING(link);
  
  /* They are marked as parameter such that the user can change them 
   * later */
  PARAMETER(eps);
  PARAMETER(kappa);
 
  PARAMETER_VECTOR(b);
  PARAMETER_VECTOR(theta);
  
  /* checks */
  {
    unsigned const n = tobs.size();
    auto check_rows = [n](matrix<Type> const &x, char const *msg){
      if(x.rows() != n)
        error(msg);
    };
    check_rows(X , "invalid 'X'");
    check_rows(XD, "invalid 'XD'");
    check_rows(Z , "invalid 'Z'");
    if(n != event.size())
      error("invalid 'event'");
    if(n != grp.size())
      error("invalid 'grp'");
    
    if(b.size()  != X.cols())
      error("invalid 'b'");
    if(XD.cols() != X.cols())
      error("invalid 'XD' (# columns)");
    
    if(Z.cols() != get_rng_dim(theta))
      error("invalid 'Z' (# columns: %d %d)", Z.cols(), get_rng_dim(theta));
  }
  
  /* perform approximations using method descibed at 
   *   https://github.com/kaskr/adcomp/issues/233#issuecomment-306032192 
   * 
   * each branch may contain further parameters or data objects */
  if(app_type == "Laplace"){
    PARAMETER_MATRIX(u);
    return laplace(COMMON_CALL, u);
    
  } else if(app_type == "GVA"){
    PARAMETER_VECTOR(theta_VA);
    DATA_INTEGER(n_nodes);
    return GVA(COMMON_CALL, theta_VA, n_nodes);
    
  } else if(app_type == "SNVA"){
    PARAMETER_VECTOR(theta_VA);
    DATA_INTEGER(n_nodes);
    DATA_STRING(param_type)
    return SNVA(COMMON_CALL, theta_VA, n_nodes, param_type);
    
  }
  
  error("approximation method '%s' is not implemented", app_type.c_str());
  return Type(0);
}
