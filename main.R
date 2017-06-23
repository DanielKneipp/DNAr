library(deSolve)
library(stringr)

get_first_part <- function(react_str) {
    return(sub('->.*', '', react_str))
}

get_second_part <- function(react_str) {
    return(sub('.*->', '', react_str))
}

is_bimolecular <- function(react_str) {
    first_part <- get_first_part(react_str)
    return(get_stoichiometry_part(first_part) == 2)
}

get_onespecies_count <- function(one_species, reaction_part) {
    m <- gregexpr(paste('[1-9]*', one_species, sep = ''), reaction_part)
    matches <- regmatches(reaction_part, m)
    nums <- array(0, length(matches[[1]]))
    if(length(nums) > 0) {
        for(i in 1:length(nums)){
            nums[i] <- as.numeric(sub(one_species, '', matches[[1]][i]))
            if(is.na(nums[i])) {
                nums[i] <- 1
            }
        }
    }
    return(nums)
}

get_stoichiometry_onespecies <- function(one_species, reaction) {
    f_p <- get_first_part(reaction)
    s_p <- get_second_part(reaction)
    
    f_p_n <- get_onespecies_count(one_species, f_p)
    s_p_n <- get_onespecies_count(one_species, s_p)
    
    r <- list(left_sto = sum(f_p_n), right_sto = sum(s_p_n))
    return(r)
}

get_stoichiometry_part <- function(reaction_part) {
    matches <- str_match_all(reaction_part, '([1-9])*([a-zA-Z0-9_]+)')[[1]]
    count <- 0
    for(i in 1:dim(matches)[1]) {
        n <- as.numeric(matches[i,2])
        if(is.na(n)) {
            n <- 1
        }
        count <- count + n
    }
    return(count)
}

get_stoichiometry_all <- function(reaction) {
    f_p <- get_first_part(react_str)
    s_p <- get_second_part(reaction)
    r <- list(left_sto = get_stoichiometry_part(f_p), 
              right_sto = get_stoichiometry_part(s_p))
    return(r)
}

remove_stoichiometry <- function(species) {
    no_sto_spec <- c()
    for(i in 1:length(species)) {
        s <- sub('^[1-9]*', '', species[i])
        no_sto_spec <- c(no_sto_spec, s)
    }
    return(no_sto_spec)
}

get_species <- function(reaction_part) {
    specs <- strsplit(reaction_part, '[^a-zA-Z0-9_]')[[1]]
    specs <- specs[specs != '']
    specs <- remove_stoichiometry(specs)
    return(specs)
}

reactants_in_reaction <- function(species, reaction) {
    f_p <- get_first_part(reaction)
    words <- get_species(f_p)
    r <- c()
    for(i in 1:length(species)) {
        if(any(words == species[i])) {
            r <- c(r, i)
        }
    }
    return(r)
}

#
# Reaction must be unimolecular or bimolecular
#  E.g.: A + B -> C
#        2A -> B
#
# Limitations:
#  - It only supports uni or bimolecular reactions
#  - Bidiretional reactions is not supported (two separated reactions are needed)
#
react <- function(species, ci, reactions, ki, t) {
    products <- matrix(data = 0, nrow = length(reactions), ncol = length(species))
    reactants <- matrix(data = 0, nrow = length(reactions), ncol = length(species))
    
    for(i in 1:length(reactions)) {
        for(j in 1:length(species)) {
            stoc <- get_stoichiometry_onespecies(species[j], reactions[i])
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

plot_behaviour <- function(behaviour) {
    matplot(x = behaviour[,1], y = behaviour[,2:dim(behaviour)[2]])
}

#
# 4-domain DNA approach
#

# Only accepts uni or bimolecular reactions 
check_reaction_4domain <- function(reaction) {
    left_part <- get_first_part(reaction)
    sto <- get_stoichiometry_part(left_part)
    return(sto < 3)
}

react_4domain <- function(species, ci, reactions, qi, qmax, cmax, t) {
    for(r in reactions) {
        if(!check_reaction_4domain(r)) {
            stop(paste('Failed to processo reaction', r))
        }
    }
    
    new_reactions <- c()
    new_species <- c(species)
    new_ks <- c()
    new_cis <- c(ci)
    
    uni_aux_species <- c('G', 'O', 'T')
    bi_aux_species <- c('L', 'H', 'W', 'O', 'T')
    
    for(i in 1:length(reactions)) {
        new_reactions_for_i <- c()
        new_species_for_i <- c() 
        new_ks_for_i <- c()
        new_cis_for_i <- c()
        
        if(is_bimolecular(reactions[i])) {
            aux <- paste(bi_aux_species, as.character(i), sep = '')
            
            left_part <- get_first_part(reactions[i])
            l_p_specs <- get_species(left_part)
            
            # For a species, get its count on left.
            # If the count is 2 (2A), transform it into 2 molecules (A + A).
            # If there is only one species and the reaction is bimolecular, 
            # this unique molecule needs to be duplicated
            if(length(l_p_specs) == 1) {
                l_p_specs <- c(l_p_specs, l_p_specs)
            }
            
            # Get the products string.
            # This time we don't need to expand 2A to A + A
            right_part <- get_second_part(reactions[i])

            # Combine the molecules into reactions
            new_reactions_for_i <- c(paste(l_p_specs[1],' +', aux[1], '-->', aux[2], '+', aux[3]),
                                     paste(aux[2],' +', aux[3], '-->', l_p_specs[1], '+', aux[1]),
                                     paste(l_p_specs[2],' +', aux[2], '-->', aux[4]),
                                     paste(aux[4],' +', aux[5], '-->', right_part))
            
            # Set the new species, initial concentrations and ks
            new_species_for_i <- c(aux)
            #                     L      H     B       O     T
            new_cis_for_i <- c(cmax[i], 0.0, cmax[i], 0.0, cmax[i])
            new_ks_for_i <- c(qi[i], qmax[i], qmax[i], qmax[i])
        } else {
            aux <- paste(uni_aux_species, as.character(i), sep = '')
            
            left_part <- get_first_part(reactions[i])
            l_p_specs <- get_species(left_part)
            right_part <- get_second_part(reactions[i])
            
            new_reactions_for_i <- c(paste(l_p_specs,' +', aux[1], '-->', aux[2]),
                                     paste(aux[2],' +', aux[3], '-->', right_part, '+', aux[1]))
            new_species_for_i <- c(aux)
            #                     G      O     T
            new_cis_for_i <- c(cmax[i], 0.0, cmax[i])
            new_ks_for_i <- c(qi[i], qmax[i])
        }
        
        new_species <- c(new_species, new_species_for_i)
        new_reactions <- c(new_reactions, new_reactions_for_i)
        new_ks <- c(new_ks, new_ks_for_i)
        new_cis <- c(new_cis, new_cis_for_i)
    }
    
    # Run the reaction
    b <- react(
        species   = new_species,
        ci        = new_cis,
        reactions = new_reactions,
        ki        = new_ks,
        t         = t
    )
    
    # Return the behaviour of the input species, without the auxyliar ones
    return(b[,1:(length(species) + 1)])
}

#
# Examples
#

run_ApBeC <- function() {
    behaviour <- react(
        species   = c('A2', 'B', 'C'),
        ci        = c(1e3, 1e3, 0), 
        reactions = c('A2 + B -> C'),
        ki        = c(1e-7),
        t         = seq(0, 72000, 10)
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
#behaviour <- consensus()

#plot_behaviour(behaviour)


#
# 4-domain Examples
#

run_ApBeC_4d <- function() {
    behaviour <- react_4domain(
        species   = c('A', 'B', 'C'),
        ci        = c(1e3, 1e3, 0),
        reactions = c('A + B -> C'),
        qi        = c(1e-7),
        qmax      = c(1e-3),
        cmax      = c(1e5),
        t         = seq(0, 72000, 10)
    )
}

run_AeB_4d <- function() {
    behaviour <- react_4domain(
        species   = c('A', 'B'),
        ci        = c(1e4, 0),
        reactions = c('A -> B'),
        qi        = c(5e-5 / 1e5),
        qmax      = c(1e-5),
        cmax      = c(1e5),
        t         = seq(0, 72000, 10)
    )
}

#behaviour <- run_ApBeC_4d()
behaviour <- run_AeB_4d()

plot_behaviour(behaviour)
