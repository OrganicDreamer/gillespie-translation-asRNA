import numpy as np
import matplotlib.pyplot as plt

def ribo_seq(steadytrace_time,steadytrace_mrna,final_mrna,steadykinetic_rates):

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
        seq_ribo = np.append(seq_ribo, count_ribo)

    # plot sequestered ribosomes over time
    plt.figure()
    plt.step(steadytrace_time, seq_ribo)
    plt.xlabel('Time (s)')
    plt.ylabel('Number of sequestered ribosomes on mRNA strand')
    plt.title('Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Sequestered ribosomes over time' % (
        steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2], steadykinetic_rates[3]))

    plt.savefig(
        'C:/Users/fujitsu/Desktop/Y4 Project/1mR_stat/Sequestered ribosomes, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' % (
        steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2],
        steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure
