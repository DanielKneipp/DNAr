library(deSolve)

get_first_part <- function(react_str) {
    return(sub('->.*', '', react_str))
}

get_second_part <- function(react_str) {
    return(sub('.*->', '', react_str))
}

is_bimolecular <- function(react_str) {
    first_part <- get_first_part(react_str)
    return(grepl('\\+', first_part))
}

get_species_count <- function(species, reaction_part) {
    m <- gregexpr(paste('[1-9]*', species, sep = ''), reaction_part)
    matches <- regmatches(reaction_part, m)
    nums <- array(0, length(matches[[1]]))
    if(length(nums) > 0) {
        for(i in 1:length(nums)){
            nums[i] <- as.numeric(sub(species, '', matches[[1]][i]))
            if(is.na(nums[i])) {
                nums[i] <- 1
            }
        }
    }
    return(nums)
}

get_stoichiometry <- function(species, reaction) {
    f_p <- get_first_part(reaction)
    s_p <- get_second_part(reaction)
    
    f_p_n <- get_species_count(species, f_p)
    s_p_n <- get_species_count(species, s_p)
    
    r <- list(left_sto = sum(f_p_n), right_sto = sum(s_p_n))
    return(r)
}

reactants_in_reaction <- function(species, reaction) {
    f_p <- get_first_part(reaction)
    words <- strsplit(f_p, '[^a-zA-Z]')[[1]]
    r <- c()
    for(i in 1:length(species)) {
        if(any(words == species[i])) {
            r <- c(r, i)
        }
    }
    return(r)
}

# Reaction must be unimolecular or bimolecular
#  E.g.: A + B -> C
#        2A -> B
react <- function(species, ci, reactions, ki, t) {
    products <- matrix(data = 0, nrow = length(reactions), ncol = length(species))
    reactants <- matrix(data = 0, nrow = length(reactions), ncol = length(species))
    
    for(i in 1:length(reactions)) {
        for(j in 1:length(species)) {
            stoc <- get_stoichiometry(species[j], reactions[i])
            products[i,j] <- stoc$right_sto
            reactants[i,j] <- stoc$left_sto
        }
    }
    
    M <- products - reactants
    Mt <- t(M)
    
    fx <- function(t, y, parms) {
        dy <- numeric(length(y))
        
        v <- matrix(data = 0, nrow = length(reactions), ncol = 1)
        for(i in 1:length(reactions)) {
            s_is_reac <- reactants_in_reaction(species, reactions[i])
            v[i,1] <- ki[i] * prod(y[s_is_reac])
        }
        
        dy <- Mt %*% v
        
        return(list(dy))
    }
    
    result <- ode(times = t, y = ci, func = fx, parms = NULL)
    return(result)
}

behaviour <- react(
    species   = c('A', 'B'),
    ci        = c(20, 13), 
    reactions = c('A + B -> 2B',
                  'A -> 2A',
                  'B -> 0'),
    ki        = c(1e-4,
                  5e-3,
                  5e-3),
    t         = seq(0, 12000, 1)
)

matplot(x = behaviour[,1], y = behaviour[,2:dim(behaviour)[2]])
