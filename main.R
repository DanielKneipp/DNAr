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

react_4domain <- function(species, ci, reactions, ki, t) {
    
}

plot_behaviour <- function(behaviour) {
    matplot(x = behaviour[,1], y = behaviour[,2:dim(behaviour)[2]])
}

run_ApBeC <- function() {
    behaviour <- react(
        species   = c('A', 'B', 'C'),
        ci        = c(1e3, 1e3, 0), 
        reactions = c('A + B -> C'),
        ki        = c(1e-7),
        t         = seq(0, 72000, 1)
    )
    return(behaviour)
}

run_lotka <- function() {
    behaviour <- react(
        species   = c('A', 'B'),
        ci        = c(2, 1), 
        reactions = c('A + B -> 2B',
                      'A -> 2A',
                      'B -> 0'),
        ki        = c(1.5,
                      1,
                      1),
        t         = seq(0, 45, 0.1)
    )
    return(behaviour)
}

run_lotka_scaled <- function() {
    behaviour <- react(
        species   = c('A', 'B'),
        ci        = c(20e-9, 10e-9), 
        reactions = c('A + B -> 2B',
                      'A -> 2A',
                      'B -> 0'),
        ki        = c(5e5,
                      1/300,
                      1/300),
        t         = seq(0, 12600, 1)
    )
    return(behaviour)
}

run_ApBeC_4domain <- function() {
    behaviour <- react(
        species   = c('A', 'B', 'C', 'L', 'H', 'W', 'O', 'T'),
        ci        = c(1e3, 1e3, 0.0, 1e5, 0.0, 1e5, 0.0, 1e5), 
        reactions = c('A + L -> H + W',
                      'H + W -> A + L',
                      'B + H -> O',
                      'O + T -> C'),
        ki        = c(1e-7,
                      1e-3,
                      1e-3,
                      1e-3),
        t         = seq(0, 72000, 1)
    )
    return(behaviour)
}

run_origonator <- function() {
    # It is not working properly yet
    behaviour <- react(
        species   = c('A', 'B', 'C'),
        ci        = c(0.8e-9, 4.5e-9, 0.8e-9), 
        reactions = c('B -> A',
                      'B + A -> 0',
                      'A -> 2A + C',
                      '2A -> 0',
                      'C -> B',
                      'C -> 0'),
        ki        = c(3.889e-7,
                      5e5,
                      0.00232,
                      12500,
                      0.00198,
                      1.195e-5),
        t         = seq(0, 360000, 1)
    )
    return(behaviour)
}

run_rossler <- function() {
    # It is not working properly yet
    behaviour <- react(
        species   = c('A', 'B', 'C'),
        ci        = c(1.8e-9, 1.8e-9, 1.8e-9), 
        reactions = c('A -> 2A',
                      '2A -> A',
                      'B + A -> 2B',
                      'B -> 0',
                      'A + C -> 0',
                      'C -> 2C',
                      '2C -> C'),
        ki        = c(0.002,
                      200000,
                      400000,
                      0.000667,
                      400000,
                      0.001,
                      200000),
        t         = seq(0, 144000, 1)
    )
    return(behaviour)
}

consensus <- function() {
    behaviour <- react(
        species   = c('X', 'Y', 'B'),
        ci        = c(0.7 * 80e-9, 0.3 * 80e-9, 0.0),
        reactions = c('X + Y -> 2B',
                      'B + X -> 2X',
                      'B + Y -> 2Y'),
        ki        = c(2e3, 2e3, 2e3),
        t         = seq(0, 54000, 5)
    )
}

#behaviour <- run_ApBeC()
#behaviour <- run_lotka()
#behaviour <- run_lotka_scaled()
#behaviour <- run_ApBeC_4domain()
#behaviour <- run_origonator()
#behaviour <- run_rossler()
behaviour <- consensus()

plot_behaviour(behaviour)
