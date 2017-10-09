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
