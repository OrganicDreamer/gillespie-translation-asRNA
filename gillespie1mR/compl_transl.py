import matplotlib.pyplot as plt
import numpy as np

def transl_count(steadytrace_time,steadytrace_mrna,final_mrna,steadykinetic_rates):

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

            count_transl = np.append(count_transl, count_transl[len(count_transl) - 1] + 1)

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
    plt.savefig(
        'C:/Users/fujitsu/Desktop/Y4 Project/1mR_dyn/Completed Translations, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' % (
        steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2],
        steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure