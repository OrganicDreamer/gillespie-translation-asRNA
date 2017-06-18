import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# array to store steady state queue length over each iteration:
plot_queue_len = np.array([])

# store average flux of ribosomes in/out of queue at steady state:
plot_ribo_flux = np.array([])

# store deterministic prediction of ribosome flux in/out at steady state:
plot_det_flux = np.array([])

num_param_iterations = 50

# iterate over parameter sweep
for i in range(0,num_param_iterations):

    #load each iteration of the parameter sweep
    param_iteration = np.load('%dinf-sweep.npz' %(i))
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
    sweep_fin_transl = param_iteration['sweep_fin_transl']

    #####################################################
    # EXTRACTING PLOTS OF OCCUPANCY DISTRIBUTION

    # vector to keep track to total time site on mRNA is occupied
    rna_occ_rec = np.zeros((len(final_mrna)))

    # vector of all positions along mRNA
    rna_loc = np.arange(1, len(rna_occ_rec) + 1, 1)

    # mRNA site occupancy by ribosome
    # iterate over all state starting from 2nd one
    for j in range(1, len(steadytrace_time)):

        # time spent in previous state
        time_in_state = steadytrace_time[j] - steadytrace_time[j - 1]

        # extract the previous state's rna occupancy
        state_of_rna = steadytrace_mrna[((j - 1) * len(final_mrna)): j * len(final_mrna)]

        # record time that was spent with rbs/codons on the rna occupied
        for k in range(len(state_of_rna)):

            if state_of_rna[k] == 1:
                rna_occ_rec[k] = rna_occ_rec[k] + time_in_state

    # save time rna occupied by ribosome
    ribo_occ = np.copy(rna_occ_rec)

    # vector cleared to keep track to total time site on mRNA is occupied
    rna_occ_rec = np.zeros((len(final_mrna)))

    # mRNA site occupancy by asRNA
    for j in range(1, len(steadytrace_time)):

        # time spent in previous state
        time_in_state = steadytrace_time[j] - steadytrace_time[j - 1]

        # extract the previous state's rna occupancy
        state_of_rna = steadytrace_mrna[((j - 1) * len(final_mrna)): j * len(final_mrna)]

        # record time that was spent with rbs occupied by asRNA
        for k in range(len(state_of_rna)):

            if state_of_rna[k] == -1:
                rna_occ_rec[k] = rna_occ_rec[k] + time_in_state

    # save time rna occupied by asrna over all trials
    asrna_occ = np.copy(rna_occ_rec)

    # plot time spent at each spot on the rna
    plt.figure()
    plt.subplot(2,1,1)
    bar_width = 0.25
    plt.bar(rna_loc, ribo_occ, bar_width, hold=True)
    plt.xlabel('Position along mRNA')
    plt.ylabel('Time spent occupied (s)')
    plt.title('Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Ribosome occupancy along mRNA' %(steadykinetic_rates[0],steadykinetic_rates[1],steadykinetic_rates[2],steadykinetic_rates[3]))

    # plot time spent at each spot on the rna
    plt.subplot(2,1,2)
    bar_width = 0.25
    plt.bar(rna_loc, asrna_occ, bar_width, hold=True)
    plt.xlabel('Position along mRNA')
    plt.ylabel('Time spent occupied (s)')
    plt.title('asRNA occupancy along mRNA')

    # save figure
    plt.tight_layout()
    plt.savefig('Occupancy distribution, inf iter%d, %fs, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' %(i,window_end_time,steadykinetic_rates[0],steadykinetic_rates[1],steadykinetic_rates[2],steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure

    ######################################################

    # EXTRACTING PLOT OF SEQUESTERED RIBOSOMES OVER TIME

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

    # plot sequestered ribosomes over time
    plt.figure()
    plt.step(steadytrace_time, seq_ribo)
    plt.xlabel('Time (s)')
    plt.ylabel('Number of sequestered ribosomes on mRNA strand')
    plt.title('Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Sequestered ribosomes over time' % (
        steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2], steadykinetic_rates[3]))

    plt.tight_layout()
    plt.savefig(
        'Sequestered ribosomes, inf iter%d, %fs, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' %(i,window_end_time,steadykinetic_rates[0],steadykinetic_rates[1],steadykinetic_rates[2],steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure

    ##########################################################
    #EXTRACTING GRAPH OF QUEUE LENGTH DISTRIBUTION

    # find number of unique queue lengths
    uni_queue_len = np.unique(seq_ribo)

    # vector to hold time spent by strand with specific size of queue
    queue_len_occ = np.zeros(len(uni_queue_len))

    # find accumulated time spent with strand having these queue lengths
    for j in range(1,len(steadytrace_time)):

        # time spent in previous state
        time_in_state = steadytrace_time[j] - steadytrace_time[j - 1]

        # queue length during time spent in previous state
        seq_ribo_in_state = seq_ribo[j-1]

        for k in range(len(uni_queue_len)):

            if seq_ribo_in_state == uni_queue_len[k]:

                queue_len_occ[k] = queue_len_occ[k] + time_in_state

    # plot accumulated time over steady state period that mRNA has
    plt.figure()
    bar_width = 0.75
    plt.bar(uni_queue_len,queue_len_occ,bar_width)
    plt.xlabel('Queue length (ribosomes)')
    plt.ylabel('Accumulated time over steady state period (s)')
    plt.title('Distribution of queue length over steady state')

    # save figure
    plt.tight_layout()
    plt.savefig(
        'Queue Length Distribution inf iter%d, %fs, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' % (i, window_end_time, steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2],
    steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure

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

    # plot number of completed translations over time
    plt.figure()
    plt.step(steadytrace_time, count_transl)
    plt.ylabel('Number of completed mRNA translations')
    plt.xlabel('Time (s)')
    plt.title('Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Translations completed over time interval' % (
        steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2], steadykinetic_rates[3]))

    # save figure
    plt.tight_layout()
    plt.savefig('Completed Translations, inf iter%d, %fs, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' %(i,window_end_time,steadykinetic_rates[0],steadykinetic_rates[1],steadykinetic_rates[2],steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure

    #####################################################
    # CALCULATE QUEUE LENGTH

    # calculate time-weighted average of number of sequestered ribosomes (ie. queue length)
    queue_len = 0

    for x in range(1, len(steadytrace_time)):
        time_weight = (steadytrace_time[x] - steadytrace_time[x - 1]) / (window_end_time - window_start_time)
        queue_len = queue_len + seq_ribo[x] * time_weight

    # store result for this iteration
    plot_queue_len = np.append(plot_queue_len, queue_len)

    ########################################
    # CALCULATE RIBOSOME FLUX OUT (ie. IN)

    # essentially the rate at which translations of mRNA are completed
    ribo_flux = (count_transl[len(count_transl)-1] - count_transl[0])/(window_end_time - window_start_time)

    # store result for this iteration
    plot_ribo_flux = np.append(plot_ribo_flux,ribo_flux)

    # deterministic ribosome flux out:
    det_flux = sweep_fin_transl[i]*(transl_rate/(sweep_fin_transl[i]+transl_rate))

    # store result for this iteration
    plot_det_flux = np.append(plot_det_flux,det_flux)

######################################################################
# GENERATE AND SAVE PARAMETER SWEEP PLOTS
plt.figure()
plt.plot(sweep_fin_transl,plot_queue_len)
plt.ylabel('Number of ribosomes in queue')
plt.xlabel('Slow codon translation rate')
plt.title('Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Queue length vs. slow codon rate' % (
steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2], steadykinetic_rates[3]))

plt.tight_layout()
plt.savefig(
'Queue length sweep, inf, %fs, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' % (
window_end_time, steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2],
steadykinetic_rates[3]))  # save the figure to file
plt.close()  # close the figure


plt.figure()
plt.plot(sweep_fin_transl,plot_ribo_flux, 'b',label='Simulated flux', hold=True)
plt.plot(sweep_fin_transl,plot_det_flux,'r',label='Deterministic flux', hold=True)
plt.legend(loc='upper right')
plt.ylabel('Ribosome flux out')
plt.xlabel('Slow codon translation rate')
plt.title('Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Flux vs. slow codon rate' % (
steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2], steadykinetic_rates[3]))

plt.tight_layout()
plt.savefig(
'Flux sweep, inf, %fs, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' % (window_end_time, steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2],
steadykinetic_rates[3]))  # save the figure to file
plt.close()  # close the figure