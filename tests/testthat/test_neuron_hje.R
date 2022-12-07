library(DNAr)
context('Hjelmfelt neuron')

run_reaction <- function(parms, filename, save = FALSE) {
    behavior <- do.call(react, parms)
    if(save) {
        save_behavior_csv(behavior, filename)
    } else {
        behavior_ok <- load_behavior_csv(filename)
        return(list(behavior, behavior_ok))
    }
}

test_that(
    'get_neuron_hje gives a neuron with the correct reaction behavior',
    {
        # Activated neuron
        n <- DNAr::get_neuron_hje('n', 1.1)
        n$t <- seq(0, 1, length.out = 100)
        n$name <- NULL
        b <- run_reaction(n, 'data/neuron_hje/neuron_hje_activated')
        expect_equal(b[[1]], b[[2]], tolerance = 1e-1)

        # Deactivated neuron
        n <- DNAr::get_neuron_hje('n', 0.9)
        n$t <- seq(0, 1, length.out = 100)
        n$name <- NULL
        b <- run_reaction(n, 'data/neuron_hje/neuron_hje_deactivated')
        expect_equal(b[[1]], b[[2]], tolerance = 1e-1)
    }
)

test_that(
    'get_neuron_AND_gate_hje gives a gate with the correct reaction behavior',
    {
        get_gate_crn <- function(input) {
            g <- DNAr::get_neuron_AND_gate_hje('_AND', input)
            and_crn <- DNAr::get_crn_from_neuron_gate_hje(g)
            and_crn$t <- seq(0, 1, length.out = 100)
            return(and_crn)
        }

        # Input 0 0
        and_crn <- get_gate_crn(c(0.5, 0.5))
        b <- run_reaction(and_crn, 'data/neuron_hje/neuron_hje_AND_00')
        expect_equal(b[[1]], b[[2]], tolerance = 1e-1)

        # Input 0 1
        and_crn <- get_gate_crn(c(0.5, 1.5))
        b <- run_reaction(and_crn, 'data/neuron_hje/neuron_hje_AND_01')
        expect_equal(b[[1]], b[[2]], tolerance = 1e-1)

        # Input 1 0
        and_crn <- get_gate_crn(c(1.5, 0.5))
        b <- run_reaction(and_crn, 'data/neuron_hje/neuron_hje_AND_10')
        expect_equal(b[[1]], b[[2]], tolerance = 1e-1)

        # Input 1 1
        and_crn <- get_gate_crn(c(1.5, 1.5))
        b <- run_reaction(and_crn, 'data/neuron_hje/neuron_hje_AND_11')
        expect_equal(b[[1]], b[[2]], tolerance = 1e-1)
    }
)

test_that(
    'if a circuit example with two AND gates behaves correctly',
    {
        get_circ_crn <- function(ins_g_1, ins_g_2) {
            g1    <- DNAr::get_neuron_AND_gate_hje('_AND1_', ins_g_1)
            g2    <- DNAr::get_neuron_AND_gate_hje('_AND2_', ins_g_2)
            g2    <- DNAr::define_neuron_gate_binding_hje(g1, g2, 1)
            circ  <- DNAr::get_circuit_from_neuron_gates_hje(list(g1, g2))
            crn   <- DNAr::get_crn_from_neuron_circuit_hje(circ)
            crn$t <- seq(0, 1, length.out = 100)
            return(crn)
        }
        # Input  0 0 0
        crn <- get_circ_crn(c(0.5, 0.5), c(0, 0.5))
        b <- run_reaction(crn, 'data/neuron_hje/neuron_hje_AND_AND_000')
        expect_equal(b[[1]], b[[2]], tolerance = 1e-1)

        # Input  0 0 1
        crn <- get_circ_crn(c(0.5, 1.5), c(0, 0.5))
        b <- run_reaction(crn, 'data/neuron_hje/neuron_hje_AND_AND_001')
        expect_equal(b[[1]], b[[2]], tolerance = 1e-1)

        # Input  0 1 1
        crn <- get_circ_crn(c(1.5, 1.5), c(0, 0.5))
        b <- run_reaction(crn, 'data/neuron_hje/neuron_hje_AND_AND_011')
        expect_equal(b[[1]], b[[2]], tolerance = 1e-1)

        # Input  1 1 1
        crn <- get_circ_crn(c(1.5, 1.5), c(0, 1.5))
        b <- run_reaction(crn, 'data/neuron_hje/neuron_hje_AND_AND_111')
        expect_equal(b[[1]], b[[2]], tolerance = 1e-1)
    }
)

test_that(
    'if a majority gate behaves correctly',
    {
        # 3 inputs
        g <- get_neuron_majority_gate('_MAJ_', 3, c(1.5, 1.5, 0))
        crn <- DNAr::get_crn_from_neuron_gate_hje(g)
        crn$t <- seq(0, 10, length.out = 100)
        b <- run_reaction(crn, 'data/neuron_hje/neuron_hje_MAJ_3i')
        expect_equal(b[[1]], b[[2]], tolerance = 1e-1)

        # 5 inputs
        g <- get_neuron_majority_gate('_MAJ_', 5, c(1.5, 1.5, 1.5, 0, 0))
        circ <- DNAr::get_circuit_from_neuron_gates_hje(list(g))
        crn <- DNAr::get_crn_from_neuron_circuit_hje(circ)
        crn$t <- seq(0, 10, length.out = 100)
        b <- run_reaction(crn, 'data/neuron_hje/neuron_hje_MAJ_5i')
        expect_equal(b[[1]], b[[2]], tolerance = 1e-1)

        # 7 inputs
        g <- get_neuron_majority_gate('_MAJ_', 7, c(1.5, 1.5, 1.5, 1.5, 0, 0, 0))
        crn <- DNAr::get_crn_from_neuron_gate_hje(g)
        crn$t <- seq(0, 10, length.out = 100)
        b <- run_reaction(crn, 'data/neuron_hje/neuron_hje_MAJ_7i')
        expect_equal(b[[1]], b[[2]], tolerance = 1e-1)
    }
)
