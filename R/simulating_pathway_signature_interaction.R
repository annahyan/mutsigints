library(here)
library(tidyverse)

p.paths = 0.5
N = 1000

simulate_pathways = function(N, prob, seed = 1435347456) {
    set.seed(seed)
    return(rbinom(N, size = 1, prob))
}


simulate_signatures = function(N, mu, zero.infl.p = 0.8, seed = 1435347456) {

    set.seed(seed)
    out = rnbinom(n = N, size = 1, mu = mu)
    
    set.seed(seed + 1)
    zero.pos = rbinom(n = N, size = 1, prob = zero.infl.p)
    
    out[zero.pos == 0] = 0
    return(out)
}

pathway.values = simulate_pathways(N, p.paths)

#' simulate signature and pathway values
#' @param N total number of samples
#' @param r the ratio of WT and MT pathways
#' @param mu1 mean signature intensity in WT
#' @param mu2 mean signature intensity in MT
#' @param zero.p1 zero-inflation probability in WT
#' @param zero.p2 zero-inflation probability in MT
#' @return 
    
simulate_exp = function(N, r, mu1, mu2, zero.p1, zero.p2) {
    
    pathway.values = simulate_pathways(N, r)
    
    wt.sigs = simulate_signatures(N = sum(pathway.values == 0),
                                  mu = mu1, zero.infl.p = zero.p1)
    mt.sigs = simulate_signatures(N = sum(pathway.values == 1),
                                  mu = mu2, zero.infl.p = zero.p2)
    
    out.df = data.frame(pathway = pathway.values, signature = 0)
    out.df[ out.df$pathway == 1, ]$signature = mt.sigs
    out.df[ out.df$pathway == 0, ]$signature = wt.sigs
    
    return(out.df)
}

N = 100
r = 0.3

mu1 = 100
mu2 = 300
zero.p1 = 0.6
zero.p2 = 0.6

sim.out = simulate_exp(N, r, mu1, mu2, zero.p1, zero.p2)


lm.robust.out = get_sig_path_lms(sim.out[, "signature", drop = FALSE],
                                 sim.out[, "pathway", drop = FALSE],
                                 sig.log = FALSE,
                                 robust = TRUE,
                                 path.to.sig = TRUE)
lm.robust.out

lm.out = get_sig_path_lms(sim.out[, "signature", drop = FALSE],
                                 sim.out[, "pathway", drop = FALSE],
                                 sig.log = FALSE,
                                 robust = FALSE,
                                 path.to.sig = TRUE)
lm.out


log.robust.out = get_sig_path_lms(sim.out[, "signature", drop = FALSE],
                                  sim.out[, "pathway", drop = FALSE],
                                  sig.log = FALSE,
                                  robust = TRUE,
                                  path.to.sig = FALSE)

log.out = get_sig_path_lms(sim.out[, "signature", drop = FALSE],
                                  sim.out[, "pathway", drop = FALSE],
                                  sig.log = FALSE,
                                  robust = FALSE,
                                  path.to.sig = FALSE)


lm.sim = lm(signature ~ pathway, data = sim.out)
lm.sim.summary = summary(lm.sim)

lm.sim.summary$residuals %>% 
    enframe %>% 
    ggplot(aes(x = value) ) + geom_histogram()


log.sim = glm(pathway ~ signature, data = sim.out, family = binomial)
log.sim.summary = summary(log.sim)

log.sim.summary$

x = runif(100)
y = 100 * x + 10 * rnorm(100)

y.x = lm(y ~ x)
summary(y.x)$coefficients

x.y = lm(x ~ y)
summary(x.y)$coefficients


