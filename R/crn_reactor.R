
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

isempty_part <- function(react_part) {
    return(!is.na(suppressWarnings(as.numeric(react_part))) && as.numeric(react_part) == 0)
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

get_reactants <- function(reaction) {
    f_p <- get_first_part(reaction)
    reactants <- get_species(f_p)
    return(reactants)
}

reactants_in_reaction <- function(species, reaction) {
    words <- get_reactants(reaction)
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
