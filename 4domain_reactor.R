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
