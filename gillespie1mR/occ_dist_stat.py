import numpy as np
import matplotlib.pyplot as plt

def occ_dist(steadykinetic_rates,final_mrna,steadytrace_time,steadytrace_mrna):

    # EXTRACTING PLOTS OF OCCUPANCY DISTRIBUTION

    # vector to keep track to total time site on mRNA is occupied
    rna_occ_rec = np.zeros((len(final_mrna)))

    # vector of all positions along mRNA
    rna_loc = np.arange(1, len(rna_occ_rec) + 1, 1)

    # mRNA site occupancy by ribosome
    for j in range(1, len(steadytrace_time)):

        # time spent in state
        time_in_state = steadytrace_time[j] - steadytrace_time[j - 1]

        # extract the state's rna occupancy
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

        # time spent in state
        time_in_state = steadytrace_time[j] - steadytrace_time[j - 1]

        # extract the state's rna occupancy
        state_of_rna = steadytrace_mrna[((j - 1) * len(final_mrna)): j * len(final_mrna)]

        # record time that was spent with rbs occupied by asRNA
        for k in range(len(state_of_rna)):

            if state_of_rna[k] == -1:
                rna_occ_rec[k] = rna_occ_rec[k] + time_in_state

    # save time rna occupied by asrna over all trials
    asrna_occ = np.copy(rna_occ_rec)

    # plot time spent at each spot on the rna
    plt.figure()
    plt.subplot(2, 1, 1)
    bar_width = 1 / 1.5
    plt.bar(rna_loc, ribo_occ, bar_width, hold=True)
    plt.xlabel('Position along mRNA')
    plt.ylabel('Time spent occupied (s)')
    plt.title('Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Ribosome occupancy along mRNA' % (
    steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2], steadykinetic_rates[3]))

    # plot time spent at each spot on the rna
    plt.subplot(2, 1, 2)
    bar_width = 0.75
    plt.bar(rna_loc, asrna_occ, bar_width, hold=True)
    plt.xlabel('Position along mRNA')
    plt.ylabel('Time spent occupied (s)')
    plt.title('asRNA occupancy along mRNA')

    # save figure
    plt.tight_layout()
    plt.savefig(
        'C:/Users/fujitsu/Desktop/Y4 Project/1mR_stat/Occupancy distribution, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' % (
        steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2],
        steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure
