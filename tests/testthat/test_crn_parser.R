library(DNAr)
context('CRN parser tools')

test_that('is_bimolecular checks if a reaction is bimolecular', {
    expect_equal(DNAr:::is_bimolecular('2A -> B'), TRUE)
    expect_equal(DNAr:::is_bimolecular('A + B -> C'), TRUE)
    expect_equal(DNAr:::is_bimolecular('A -> B'), FALSE)
    expect_equal(DNAr:::is_bimolecular('0 -> B'), FALSE)

    # Since 0 is considered a name of a species (for formation and degradation
    # reactions), this case will return TRUE. This reaction string is a invalid
    # construction
    expect_equal(DNAr:::is_bimolecular('0 + A -> B'), TRUE)
})

test_that('is_unimolecular checks if a reaction is unimolecular', {
    expect_equal(DNAr:::is_unimolecular('2A -> B'), FALSE)
    expect_equal(DNAr:::is_unimolecular('A + B -> C'), FALSE)
    expect_equal(DNAr:::is_unimolecular('A -> B'), TRUE)

    # This case will return TRUE because 0 is considered a special species
    expect_equal(DNAr:::is_unimolecular('0 -> B'), TRUE)
})

test_that('is_degradation checks if a reaction is of degradation', {
    expect_equal(DNAr:::is_degradation('A -> B'), FALSE)
    expect_equal(DNAr:::is_degradation('A + B -> B + C'), FALSE)
    expect_equal(DNAr:::is_degradation('2A -> 2B'), FALSE)
    expect_equal(DNAr:::is_degradation('A + B -> 0'), TRUE)
    expect_equal(DNAr:::is_degradation('A -> 0'), TRUE)
    expect_equal(DNAr:::is_degradation('2A -> 0'), TRUE)
    expect_equal(DNAr:::is_degradation('0 -> A + B'), FALSE)
})

test_that('is_formation checks if a reaction is of degradation', {
    expect_equal(DNAr:::is_formation('A -> B'), FALSE)
    expect_equal(DNAr:::is_formation('A + B -> B + C'), FALSE)
    expect_equal(DNAr:::is_formation('2A -> 2B'), FALSE)
    expect_equal(DNAr:::is_formation('0 -> A'), TRUE)
    expect_equal(DNAr:::is_formation('0 -> 2A'), TRUE)
    expect_equal(DNAr:::is_formation('0 -> A + B'), TRUE)
    expect_equal(DNAr:::is_formation('A + B -> 0'), FALSE)
})

test_that('get_species returns the species of normal and aspecial reactions', {
    expect_equal(DNAr:::get_species('A + 2B'), c('A', 'B'))
    expect_equal(DNAr:::get_species('A + A'), c('A', 'A'))
    expect_equal(DNAr:::get_species('2A'), c('A'))
    expect_equal(DNAr:::get_species('0'), c('0'))
    expect_equal(DNAr:::get_species('0 + 2A'), c('0', 'A'))
    expect_equal(DNAr:::get_species('20'), c('0'))

    test_warn <- function(test, exp_res) {
        expect_warning(r <- test())
        expect_equal(r, exp_res)
    }

    test_warn(function() {DNAr:::get_species('2')}, character(0))
    test_warn(function() {DNAr:::get_species('2 + 12')}, character(0))

    # 0 will be a 'species' with stoichiometry 2
    test_warn(function() {DNAr:::get_species('4 + 20')}, c('0'))
})

test_that('check_fix_reaction checks and fix some problematic reactions', {
    # Ok reactions
    test_ok <- function(reaction) {
        expect_equal(DNAr:::check_fix_reaction(reaction), reaction)
    }
    test_ok('A + B -> C')
    test_ok('2A2 -> 0')
    test_ok('0 -> A')
    test_ok('2A -> 0')

    # Fixable reactions
    test_warn <- function(bad_reaction, ok_reaction) {
        expect_warning(r <- DNAr:::check_fix_reaction(bad_reaction))
        expect_equal(r, ok_reaction)
    }
    test_warn('A + 0 -> B', 'A -> B')
    test_warn('0 + B -> 0', ' B -> 0')
    test_warn('0 -> 2A + 0', '0 -> 2A ')

    # Bad reactions
    test_bad <- function(reaction) {
        expect_error(DNAr:::check_fix_reaction(reaction))
    }
    test_bad('0 -> 0')
    test_bad('A + B ->')
    test_bad('A + B > C')
    test_bad('-> A + B')

    # TODO: Detect this errors
    # - `A ->> B` (In this case it will think that the
    #              last > is a species name)
    # - `A - B -> C` (Relation operator between species)
    # - `A -> B -> C` (Two or more relation operators between parts)
})

test_that('get_onespecies_count returns the correct value when
          two species have the same prefix name', {
    expect_equal(DNAr:::get_onespecies_count('A_i', 'A_i + 2A ')[[1]], 1)
    expect_equal(DNAr:::get_onespecies_count('A_i', 'A_i + A_i '), array(c(1, 1)))
    expect_equal(DNAr:::get_onespecies_count('A', 'A_i + 2A ')[[1]], 2)
})

test_that('check_crn throws an error when the number of species and initial
          concentrations are not the same', {
    # CRN with more species than cis
    crn1 <- list(
        species   = c('A', 'B'),
        ci        = c(1),
        reactions = c('A -> B',
                      'B -> 2A'),
        ki        = c(2, 2),
        t         = 0:10
    )

    expect_error(
        do.call(DNAr:::check_crn, crn1),
        'The length of species and ci are not equal',
        fixed = TRUE
    )

    # CRN with more cis than species
    crn2 <- list(
        species   = c('A'),
        ci        = c(1, 2),
        reactions = c('A -> B',
                      'B -> 2A'),
        ki        = c(2, 2),
        t         = 0:10
    )

    expect_error(
        do.call(DNAr:::check_crn, crn2),
        'The length of species and ci are not equal',
        fixed = TRUE
    )
})


test_that('check_crn throws an error when the number of reactions and ks
          are not the same', {
    # CRN with more reactions than kis
    crn1 <- list(
        species   = c('A', 'B'),
        ci        = c(1, 1),
        reactions = c('A -> B',
                      'B -> 2A'),
        ki        = c(2),
        t         = 0:10
    )

    expect_error(
      do.call(DNAr:::check_crn, crn1),
      'The length of reactions and ki are not equal',
      fixed = TRUE
    )

    # CRN with more kis than reactions
    crn2 <- list(
        species   = c('A', 'B'),
        ci        = c(1, 1),
        reactions = c('A -> B'),
        ki        = c(2, 2),
        t         = 0:10
    )

    expect_error(
        do.call(DNAr:::check_crn, crn2),
        'The length of reactions and ki are not equal',
        fixed = TRUE
    )
})
