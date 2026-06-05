//
// HBomics: Hierarchcal bayesian model for subclone driver omics events 
// Allele imblance ratio theta from enrichment SNP sites were modeled and re-parameterized by a selection strength a.
//
// This model incorporate multi-level overdispersion and combined with empirical bayes methods that borrow from 
// global enrichment data and genomics data (WGS).
//
// segmental and subclone overdispersions were considered and global nucleotide mapping bias and variance were considered as
// in hierarchical prior as hyperparameters with empirical bayes.
// 
// level 1 (top level to down): global mapping bias (tech bias) and allele ratio variance (tech + bio)
// level 2: segmental/subclone variance shared and flow to different feature/segments
// level 3: individual SNP with segmental hierarchical prior;
//
// this model also consider the nature of het. SNPs captured from features, a Gamma function used to model overdispersion was
// re-parameterized with a sigmoid SNP density function, i.e., segments with more SNP density would more rely on likelihood 
// overdispersion, while segments with extreme low SNP density would more rely on the subclone variance (which served as empirical bayes) 
// 
// theta (allele imbalance ratio upon enrichment) was re-parameterized by WGS allele ratio, which potentially remove the bias brought by
// purity, ploidy and other confounders
//
// We assume purity and ploidy are shared among WGS and other enrichment sequencing data;
// segmental variance was consider consistent among sequencing technologies;

// theta reparameterization
// incoporate with allele frequencies from WGS data to reparameterize
// observed theta with latent true selection strength on SNP site, where a is the selection strength on the SNP site and b is the 
// allele frequency observed from WGS data (no enrichment)
functions {
  real f_a_b(real a, real b) {
    return (a * b) / (a * b + (1 - a) * (1 - b));
  }
}

// data block
// data was feed site by site
// individual marginal posterior model for each site
data {
  int<lower=1> T;                                  // Total SNP size
  int<lower=1> C;                                  // Total number of subclones
  int<lower=1> S;                                  // Total number of segments
  //
  array[T] int<lower=1, upper=S> segment_id;       // SNP to segment id assignment
  array[S] int<lower=1, upper=C> subclone_id;      // segment to subclone id assignment
  array[S] int n;                                  // SNP density (counts) for each segment
  //
  array[T] int A;                                  // alternative allele count for each SNP [T]
  array[T] int N;                                  // Total allele count for each SNP [T]
  array[T] real<lower=0, upper=1> b;               // correpsonding allele ratio from WGS
  //
  array[C] real<lower=0> v_t;                      // empirical AF variance for each subclone (genome subclone based)
  array[S] real<lower=0> v_j;                      // empirical AF variance for each segment (genome segment based)
  real<lower=0> lambda;                            // Gamma scale, default is 0.2
  // hyper-parameters
  real<lower=0, upper=1> mu0;                      // Global mean of allele ratio (uniform prior/mapping bias)
  real<lower=0> sigma0;                            // Global sd of allele ratio
  real<lower=0> s;                                 // slope for scaled sigmoid, default is 0.5
  int<lower=0> c;                                  // median of SNP density in segments
}


// latent variables
parameters {
  vector<lower=1e-5, upper=1-1e-5>[T] a;           // selection strength avoid pure 0 and 1
  vector<lower=1e-5, upper=1-1e-5>[T] u;           // segmental SNP allele ratio mean
  vector<lower=1e-5, upper=1-1e-5>[T] v;           // SNP allele ratio variance
  vector<lower=1e-5, upper=1-1e-5>[T] mu;          // subclone allele ratio mean
}

transformed parameters {
  vector<lower=0, upper=1>[T] theta;
  for (id in 1:T) {
    theta[id] = f_a_b(a[id], b[id]);               // re-parameterization of theta and a
  }
}

model {
  // prior on subclone level
  // borrow global mapping bias and overdispersion as prior for each subclone
  mu ~ normal(mu0, sigma0);
  // segmental levels priors
  vector[S] w;
  vector[S] alpha_r;
  for (ids in 1:S) {
    // subclone index
    int sc_idx = subclone_id[ids];
    // logistic scaled sigmoid for weights calculating
    // this consider the sparsity of SNP density on different segments
    w[ids] = inv_logit(s * (n[ids] - c));
    // shape parameter for gamma modeling overdispersion
    // borrowing subclone overdispersion prior if sparsity is high in segment
    alpha_r[ids] = w[ids] * v_j[ids] + (1 - w[ids]) * v_t[sc_idx];
  }
  // SNP level
  for (i in 1:T) {
    // segmental id for SNP
    int sg_idx = segment_id[i];
    // hierarchical prior of subclone mean on segmental level
    u[i] ~ beta(mu[i] * v_t[subclone_id[sg_idx]], (1 - mu[i]) * v_t[subclone_id[sg_idx]]);
    v[i] ~ gamma(alpha_r[sg_idx]+1e-10, lambda); // avoid shape parameter to 0
    // did not integrate theta out here
    // segmental information helping modeling feature allele ratio
    // Likelihood
    theta[i] ~ beta(u[i] * v[i], (1 - u[i]) * v[i]);
    A[i] ~ binomial(N[i], theta[i]);
  }
}


// model evaluation block
generated quantities {
  vector[T] log_lik;         // For model comparison
  vector[T] theta_rep;       // Posterior predictive
  vector[T] overdispersion;  // Effective overdispersion at feature level
  vector[T] var_theta;       // Variance of theta_f from beta distribution
  //real sub_clone_u;     // subclone average allele ratio gap

  for (id in 1:T) {
    real alpha_theta = u[id] * v[id];
    real beta_theta = (1 - u[id]) * v[id];
    // Posterior predictive sample of theta
    theta_rep[id] = beta_rng(alpha_theta, beta_theta);
    // Overdispersion: inversely proportional to v
    // beta prior evaluation
    overdispersion[id] = 1 / (1 + v[id]);  // optional transformation
    // Variance of theta ~ Beta(alpha, beta)
    // theoretical theta around observed u
    var_theta[id] = (alpha_theta * beta_theta) /
      ((alpha_theta + beta_theta)^2 * (alpha_theta + beta_theta + 1));

    // Log-likelihood for loo/cv
    log_lik[id] = binomial_lpmf(A[id] | N[id], theta[id]);
  }
}
