import gillespie1mR.dynamic as gillsim
import gillespie1mR.occ_dist as mrna
import gillespie1mR.free_pools as fp
import gillespie1mR.compl_transl as ct
import numpy as np

# 101 values describe the state of mRNA codons + RBS
initial_mrna = np.zeros(101)

free_ribo = 150
free_asrna = 150

# time to presumably reach steady state:
window_start_time = 2000

# final time to simulate until:
window_end_time = 3000

# 2d array (n x 4) of kinetic constants where each row is: ribo bind, ribo unbind, asRNA bind, asRNA unbind:
sweep_kinetic_const = np.array(
    [

    [0.001,0.001,0,0],
    [0.001,0.001,0.15,0.1]

    ])

# sweep over different rows of kinetic constants when simulating model:
for i in range(sweep_kinetic_const.shape[0]):

    # simulate until steady state is reached
    (trace_mrna, trace_ribo, trace_asrna, trace_time, steady_mrna, steady_ribo, steady_asrna, steady_time,
     kinetic_rates) = gillsim.dyn_pool(initial_mrna, free_ribo, free_asrna, 0, window_start_time,sweep_kinetic_const[i][0],sweep_kinetic_const[i][1],sweep_kinetic_const[i][2],sweep_kinetic_const[i][3])

    # simulate from steady state onwards
    (steadytrace_mrna, steadytrace_ribo, steadytrace_asrna, steadytrace_time, final_mrna, final_ribo, final_asrna,
     final_time, steadykinetic_rates) = gillsim.dyn_pool(steady_mrna, steady_ribo, steady_asrna, steady_time,
                                                          window_end_time,sweep_kinetic_const[i][0],sweep_kinetic_const[i][1],sweep_kinetic_const[i][2],sweep_kinetic_const[i][3])

    # generate and save plots
    mrna.occ_dist(steadykinetic_rates,final_mrna,steadytrace_time,steadytrace_mrna)
    fp.free_pools(steadytrace_time,steadytrace_ribo,steadytrace_asrna,steadykinetic_rates)
    ct.transl_count(steadytrace_time,steadytrace_mrna,final_mrna,steadykinetic_rates)


