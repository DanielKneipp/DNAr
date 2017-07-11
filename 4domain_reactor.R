#
# 4-domain DNA approach
#

source('crn_reactor.R')

# Only accepts uni or bimolecular reactions 
check_reaction_4domain <- function(reaction) {
    left_part <- get_first_part(reaction)
    sto <- get_stoichiometry_part(left_part)
    return(sto < 3)
}

get_buff_modules <- function(reactions, ki, qmax, cmax) {
    sigmas <- list()
    bff_aux_species <- c('LS', 'HS', 'WS')
    global_cmax <- max(cmax)
    global_qmax <- max(qmax)
    #                 LS      HS      WS
    bff_cis <- c(global_cmax, 0, global_cmax)
    
    uni_count <- 0
    for(i in 1:length(reactions)) {
        reactants <- get_reactants(reactions[i])
        first_reactant <- reactants[[1]]
        
        if(is_bimolecular(reactions[[i]])) {
            if(is.null(sigmas[[first_reactant]])) {
                sigmas[first_reactant] <- ki[[i]]
            } else {
                sigmas[first_reactant] <- sigmas[[first_reactant]] + ki[[i]]
            }
        } else {
            uni_count <- uni_count + 1
            if(is.null(sigmas[[first_reactant]])) {
                sigmas[first_reactant] <- 0
            }
        }
    }
    
    # There is no bimolecular reactions
    if(uni_count == length(reactions)) {
        return(NULL)
    }
    
    reaction_sigma <- max(unlist(sigmas))
    lambda_1 <- global_qmax / (global_qmax - reaction_sigma)
    
    new_ks <- c()
    new_bff_reactions <- c()
    new_species <- c()
    new_cis <- c()
    for(i in 1:length(sigmas)) {
        if(!sigmas[[i]] == reaction_sigma) {
            qs <- lambda_1 * (reaction_sigma - sigmas[[i]])
            aux_specs <- paste(bff_aux_species, as.character(i), sep = '')
            input_spec <- names(sigmas)[i]
            
            forward_reaction <- paste(input_spec, '+', aux_specs[1], '-->', aux_specs[2], '+', aux_specs[3])
            backward_reaction <- paste(aux_specs[2], '+', aux_specs[3], '-->', input_spec, '+', aux_specs[1])
            
            new_bff_reactions <- c(new_bff_reactions, forward_reaction, backward_reaction)
            new_ks <- c(new_ks, qs, global_qmax)
            new_species <- c(new_species, aux_specs[1], aux_specs[2], aux_specs[3])
            new_cis <- c(new_cis, bff_cis)
        }
    }
    
    ret <- list(
        lambda_1 = lambda_1,
        new_species = new_species,
        new_cis = new_cis,
        new_reactions = new_bff_reactions,
        new_ks = new_ks
    )
}

react_4domain <- function(species, ci, reactions, ki, qmax, cmax, t) {
    for(r in reactions) {
        if(!check_reaction_4domain(r)) {
            stop(paste('Failed to process reaction', r))
        }
    }
    
    new_reactions <- c()
    new_species <- c(species)
    new_ks <- c()
    new_cis <- c(ci)
    
    uni_aux_species <- c('G', 'O', 'T')
    bi_aux_species <- c('L', 'H', 'W', 'O', 'T')
    
    # Get buffer modules
    buffer_stuff <- get_buff_modules(reactions, ki, qmax, cmax)
    
    # Change cis according to the lambda^{-1} factor
    if(!is.null(buffer_stuff)) {
        for(i in 1:length(new_cis)) {
            new_cis[i] <- new_cis[[i]] * buffer_stuff$lambda_1
        }
    }
    
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
            new_reactions_for_i <- c(paste(l_p_specs[1], '+', aux[1], '-->', aux[2], '+', aux[3]),
                                     paste(aux[2], '+', aux[3], '-->', l_p_specs[1], '+', aux[1]),
                                     paste(l_p_specs[2], '+', aux[2], '-->', aux[4]))
            
            # Set the new species, initial concentrations and ks
            new_species_for_i <- c(aux[1], aux[2], aux[3], aux[4])
            #                     L      H     B       O 
            new_cis_for_i <- c(cmax[i], 0.0, cmax[i], 0.0)
            
            # Recalculate qis according to the buffer module theory
            qi_with_buff <- ki[i]
            if(!is.null(buffer_stuff)) {
                qi_with_buff <- qi_with_buff * buffer_stuff$lambda_1
            }
            new_ks_for_i <- c(qi_with_buff, qmax[i], qmax[i])
            
            # If the right part is only a 0, do not add the last reaction, otherwise:
            if(!isempty_part(right_part)) {
                new_reactions_for_i <- c(new_reactions_for_i,
                                         paste(aux[4], '+', aux[5], '-->', right_part))
                new_species_for_i <- c(new_species_for_i, aux[5])
                #                                   T
                new_cis_for_i <- c(new_cis_for_i, cmax[i])
                new_ks_for_i <- c(new_ks_for_i, qmax[i])
            }
        } else {
            aux <- paste(uni_aux_species, as.character(i), sep = '')
            
            left_part <- get_first_part(reactions[i])
            l_p_specs <- get_species(left_part)
            right_part <- get_second_part(reactions[i])
            
            new_reactions_for_i <- c(paste(l_p_specs, '+', aux[1], '-->', aux[2]))
            new_species_for_i <- c(aux[1], aux[2])
            #                     G      O
            new_cis_for_i <- c(cmax[i], 0.0)
            qi_with_buff <- ki[i] / cmax[i]
            if(!is.null(buffer_stuff)) {
                qi_with_buff <- qi_with_buff * buffer_stuff$lambda_1
            }
            new_ks_for_i <- c(qi_with_buff)
            
            if(!isempty_part(right_part)) {
                new_reactions_for_i <- c(new_reactions_for_i,
                                         paste(aux[2], '+', aux[3], '-->', right_part))
                new_species_for_i <- c(new_species_for_i, aux[3])
                #                                   T
                new_cis_for_i <- c(new_cis_for_i, cmax[i])
                new_ks_for_i <- c(new_ks_for_i, qmax[i])
            }
        }
        
        new_species <- c(new_species, new_species_for_i)
        new_reactions <- c(new_reactions, new_reactions_for_i)
        new_ks <- c(new_ks, new_ks_for_i)
        new_cis <- c(new_cis, new_cis_for_i)
    }
    
    # Add buffer stuff
    if(!is.null(buffer_stuff)) {
        new_species <- c(new_species, buffer_stuff$new_species)
        new_reactions <- c(new_reactions, buffer_stuff$new_reactions)
        new_ks <- c(new_ks, buffer_stuff$new_ks)
        new_cis <- c(new_cis, buffer_stuff$new_cis)
    }
    
    # Run the reaction
    b <- react(
        species   = new_species,
        ci        = new_cis,
        reactions = new_reactions,
        ki        = new_ks,
        t         = t
    )
    
    # Return the behaviour of the input species, without the auxiliar ones
    return(b[,1:(length(species) + 1)])
}
