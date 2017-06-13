import numpy as np

# ARBITRARILY INCREASED asRNA binding
as_on_multiplier = 2

####################################
# LOAD BASELINE QUEUED PROTEIN PRODUCTION
baseline = np.load('fin-asrna.npz')
phi0=baseline['phi0']

#####################################

num_param_iterations = 10

# iterate over parameter sweep
for i in range(0,num_param_iterations):

    #load each iteration of the parameter sweep
    param_iteration = np.load('%diter-%dasr-fin-sweep.npz' %(i,as_on_multiplier))
    steadytrace_mrna=param_iteration['steadytrace_mrna']
    steadytrace_ribo = param_iteration['steadytrace_ribo']
    steadytrace_asrna = param_iteration['steadytrace_asrna']
    steadytrace_time = param_iteration['steadytrace_time']
    final_mrna = param_iteration['final_mrna']
    final_ribo = param_iteration['final_ribo']
    final_time = param_iteration['final_time']
    steadykinetic_rates = param_iteration['steadykinetic_rates']
    window_start_time = param_iteration['window_start_time']
    window_end_time = param_iteration['window_end_time']
    transl_rate = param_iteration['transl_rate']
    fin_transl = param_iteration['fin_transl']

    ######################################################
    # EXTRACTING SEQUESTERED RIBOSOMES OVER TIME

    seq_ribo = np.array([])

    # count sequestered number of ribosomes over time
    for j in range(0, len(steadytrace_time)):

        # extract each state's rna occupancy
        state_of_rna = steadytrace_mrna[
                       (j * len(final_mrna)): (j + 1) * len(final_mrna)]

        # initialize count of total ribosomes for the state
        count_ribo = 0

        # loop over state of rna
        for k in range(len(state_of_rna)):

            # total the number of ribosomes on the strand
            if state_of_rna[k] == 1:
                count_ribo = count_ribo + 1

        # log the number of ribosomes sequestered on the mRNA for this time interval
        seq_ribo= np.append(seq_ribo, count_ribo)

    #####################################################
    # EXTRACTING COMPLETED TRANSLATIONS OVER TIME

    # keeps track of number of completed translations at each time stamp
    count_transl = np.zeros(1)

    # extract count of completed translation
    for j in range(1, len(steadytrace_time)):

        # extract the state's stop codon
        state_stop_codon = steadytrace_mrna[((j + 1) * len(final_mrna)) - 1]

        # extract previous state's stop codon
        prev_state_stop_codon = steadytrace_mrna[(j * len(final_mrna)) - 1]

        # if a ribosome has left the codon
        if state_stop_codon == 0 and prev_state_stop_codon == 1:

            count_transl = np.append(count_transl, count_transl[len(count_transl)-1] + 1)

        else:

            count_transl = np.append(count_transl, count_transl[len(count_transl) - 1])

    #####################################################
    # CALCULATE QUEUE LENGTH

    # calculate time-weighted average of number of sequestered ribosomes (ie. queue length)
    queue_len = 0

    for x in range(1, len(steadytrace_time)):
        time_weight = (steadytrace_time[x] - steadytrace_time[x - 1]) / (window_end_time - window_start_time)
        queue_len = queue_len + seq_ribo[x] * time_weight

    ########################################
    # CALCULATE RIBOSOME FLUX OUT (ie. IN)

    # essentially the rate at which translations of mRNA are completed
    ribo_flux = (count_transl[len(count_transl)-1] - count_transl[0])/(window_end_time - window_start_time)

    # if protein production matches or exceeds baseline
    if ribo_flux>=phi0:
        np.savez(str('%dasr-on-table-entry' %(as_on_multiplier)), queue_len=queue_len,ribo_flux=ribo_flux,as_bind=steadykinetic_rates[2],as_unbind=steadykinetic_rates[3])
        break

