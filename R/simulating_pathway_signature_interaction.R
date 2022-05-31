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
    
}

pathway.values = simulate_pathways(N, p.paths)

