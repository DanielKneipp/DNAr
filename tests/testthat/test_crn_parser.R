library(DNAr)
context('CRN parser tools')

test_that('is_bimolecular checks if a reaction is bimolecular', {
    expect_equal(DNAr:::is_bimolecular('2A -> B'), TRUE)
    expect_equal(DNAr:::is_bimolecular('A + B -> C'), TRUE)
    expect_equal(DNAr:::is_bimolecular('A -> B'), FALSE)
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
