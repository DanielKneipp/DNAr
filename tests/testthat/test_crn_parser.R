library(DNAr)
context('CRN parser tools')

test_that('is_bimolecular checks if a reaction is bimolecular', {
    expect_equal(DNAr:::is_bimolecular('2A -> B'), TRUE)
    expect_equal(DNAr:::is_bimolecular('A + B -> C'), TRUE)
    expect_equal(DNAr:::is_bimolecular('A -> B'), FALSE)
})
