
# from the GP chapter

set.seed(1234)

library(tidyverse)

l = 1            # for l, sigma_f, sigma_n, see note at covariance function
sigma_f = 1      
sigma_n = .25 
k_eps   = 1e-8   # see note at Kstarstar
n_prior = 5      # number of prior draws
n_post_pred = 5  # number of posterior predictive draws


N = 25
X_train = 15 * (runif(N) - .5)  
n_train = length(X_train)

# kept sine function for comparison to noise free result
y_train = sin(X_train) + rnorm(n = n_train, sd = .1)  

X_test = seq(-7.5, 7.5, length = 200)
n_test = length(X_test)

gp_mu <- function(x) {
  map_dbl(x, function(x) x = 0)
}

gp_K <- function(
  x,
  y = NULL,
  l = 1,
  sigma_f = 1,
  sigma_n = .5
) {
  
  if(!is.null(y)){
    sigma_f * exp( -(1/(2 * l^2)) * as.matrix(dist(x, upper = TRUE, diag = TRUE) ^ 2) ) +
      sigma_n*diag(length(x))    
  }
  else{
    sigma_f * exp( -(1/(2 * l^2)) * as.matrix(dist(x, upper = TRUE, diag = TRUE) ^ 2) )
  }  
}

Ky = gp_K(
  x = X_train,
  y = y_train,
  l = l,
  sigma_f = sigma_f,
  sigma_n = sigma_n
)

# initial matrix
K_ = gp_K(
  c(X_train, X_test),
  l = l,
  sigma_f = sigma_f,
  sigma_n = sigma_n
)

Kstar  = K_[1:n_train, (n_train+1):ncol(K_)]                    # dim = N x N*
tKstar = t(Kstar)                                               # dim = N* x N
Kstarstar = K_[(n_train+1):nrow(K_), (n_train+1):ncol(K_)] +    # dim = N* x N*
  k_eps*diag(n_test)      # the k_eps part is for positive definiteness

Kyinv = solve(Ky)

post_mu = gp_mu(X_test) + tKstar %*% Kyinv %*% (y_train - gp_mu(X_train))
post_K  = Kstarstar - tKstar %*% Kyinv %*% Kstar
s2 = diag(post_K)


y_pp = data.frame(t(MASS::mvrnorm(n_post_pred, mu = post_mu, Sigma = post_K)))

pp_data = data.frame(
  x = X_test,
  y = y_pp,
  fmean = post_mu, 
  se_lower = post_mu - 2 * sqrt(s2),
  se_upper = post_mu + 2 * sqrt(s2)
) %>% 
  pivot_longer(starts_with('y'), names_to = 'variable')


ggplot(aes(x = x, y = value), data = pp_data) +
  geom_ribbon(aes(ymin = se_lower, ymax = se_upper, group = variable),
              fill = 'gray95') +
  geom_line(aes(group = variable), color = '#ff550040', size = .5) +
  geom_point(
    aes(x = X_train, y = y_train),
    data = data.frame(X_train, y_train),
    color = '#00aaff80',
    size = 2
  ) +
  labs(x = '', y = '') +
  theme_void()

ggsave('img/title_plot.svg',  bg = "transparent")
ggsave('img/title_plot.png',  bg = "transparent")
