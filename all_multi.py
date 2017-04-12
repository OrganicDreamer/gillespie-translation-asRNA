import numpy as np
import random, math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#####################################
#INITIALIZE ALL PARAMETERS HERE#####

# set specific directory to save output graphs to if desired (remember to include final '/' in directory path
save_dir = ''

# define size of ribosome and asrna pools
free_ribo = 150
free_asrna = 150

# time to presumably reach steady state:
window_start_time = 2000

# final time to simulate until:
window_end_time = 3000

# 2d array (n x 4) of kinetic constants where each row is: ribo bind, ribo unbind, asRNA bind, asRNA unbind:
sweep_kinetic_const = np.array(
    [

    [1,0,1,0.1]

    ])

#number of mRNA strands
num_strands = 1

# 101 values describe the state of mRNA codons + RBS
length_mrna = 101

# each column of 101 values describes the state of an mRNA's codons + RBS
initial_mrna = np.zeros((length_mrna,num_strands))

# Translation initiation, elongation, completion rates
start_transl = 1
transl_rate = 1
fin_transl = 0.01

# Array of values for translation elongation rates of specific codons on mRNA, not including stop codon or RBS:
codon_transl_rates = np.full((length_mrna - 2), transl_rate)

# Set specific rate for specific codon(s), in this case 50th codon:
codon_transl_rates[49] = 1

##########################################################################
# sweep over different rows of kinetic constants when simulating model:
for i in range(sweep_kinetic_const.shape[0]):

    def gillespie_multmR(mrna_state,free_ribo_pool,free_asrna_pool,start_time,sim_time, start_transl, fin_transl, codon_transl_rates, i = i):

        # SETTING UP SIMULATION AND MODEL PARAMETERS

        # array of values describe the state of mRNA codons + RBS
        current_mrna = np.copy(mrna_state)

        # Number of ribosomes:
        current_ribo = free_ribo_pool

        # Number of asRNA:
        current_asrna = free_asrna_pool

        # Initial time:
        current_time = start_time

        # Simulate until time:
        simulation_time = sim_time

        # Define rate constants:
        k_rib_bind = sweep_kinetic_const[i][0]
        rib_unbind_rate = sweep_kinetic_const[i][1]

        k_asrna_bind = sweep_kinetic_const[i][2]
        asrna_unbind_rate = sweep_kinetic_const[i][3]

        # Hold trajectory data:
        trace_asrna = np.array([current_asrna])
        trace_ribo = np.array([current_ribo])
        trace_mrna = np.array([current_mrna])
        trace_time = np.array([current_time])

        # RUNNING SIMULATION OF MODEL
        while current_time < simulation_time:

            # define rates of state transitions
            free_rib_bind_rate = k_rib_bind * current_ribo
            free_asrna_bind_rate = k_asrna_bind * current_asrna

            # initialize vector to record possible next states' mRNA, first entry is junk
            next_mrna = np.array([current_mrna])

            # initialize vector to record possible next states' ribosome pool
            next_ribo = np.array([])

            # initialize vector to record possible next states' asRNA pool
            next_asrna = np.array([])

            # initialize vector to record possible transition rates
            next_rates = np.array([])

            # initialize vector to explore possible state changes of mRNA
            possible_mrna = np.copy(current_mrna)

            # initialize variable to explore possible state changes of ribosome pool
            possible_ribo = current_ribo

            #initialize variable to explore possible state changes of asRNA pool
            possible_asrna = current_asrna

            # rules for next possible transitions and rates based on current state of mRNA strands:
            for i in range(np.shape(current_mrna)[0]):

                for n in range(np.shape(current_mrna)[1]):

                    # empty RBS
                    if (i == 0) and (current_mrna[i,n] == 0):

                        # Free ribosomes available to bind to RBS
                        if current_ribo >= 1:

                            # Free ribosome binds to RBS
                            possible_mrna[0,n] = 1
                            possible_ribo = possible_ribo - 1
                            possible_asrna = possible_asrna

                            # rate of the possible state change
                            possible_rate = free_rib_bind_rate

                            # store the rate of the new state transition
                            next_rates = np.append(next_rates,possible_rate)

                            # store the 101 values describing state's mRNA
                            next_mrna = np.append(next_mrna,[possible_mrna],0)

                            # store the value describing state's ribosome pool
                            next_ribo = np.append(next_ribo,possible_ribo)

                            # store the value describing state's asRNA pool
                            next_asrna = np.append(next_asrna,possible_asrna)

                            # revert possible state to explore other possibilities from current state
                            possible_mrna = np.copy(current_mrna)
                            possible_asrna = current_asrna
                            possible_ribo = current_ribo

                        # Free asRNA available to bind to RBS
                        if current_asrna >= 1:

                            # Free asRNA binds to RBS
                            possible_asrna = possible_asrna - 1
                            possible_ribo = possible_ribo
                            possible_mrna[0,n] = -1

                            # Rate of possible state change
                            possible_rate = free_asrna_bind_rate

                            # store the rate of the new state transition
                            next_rates = np.append(next_rates, possible_rate)

                            # store the 101 values describing the state's mRNA
                            next_mrna = np.append(next_mrna, [possible_mrna],0)

                            # store the value describing state's ribosome pool
                            next_ribo = np.append(next_ribo, possible_ribo)

                            # store the value describing state's asRNA pool
                            next_asrna = np.append(next_asrna, possible_asrna)

                            # revert possible state to explore other possibilities from current state
                            possible_mrna = np.copy(current_mrna)
                            possible_asrna = current_asrna
                            possible_ribo = current_ribo

                    # RBS occupied
                    if (i == 0) and current_mrna[0,n] != 0:

                        # By ribosome
                        if current_mrna[0,n] == 1:

                            # ribosome unbind from RBS
                            possible_mrna[0,n] = 0
                            possible_ribo = possible_ribo + 1
                            possible_asrna = possible_asrna

                            # unbind from rbs: rate of possible state change
                            possible_rate = rib_unbind_rate

                            # store the rate of the new state transition
                            next_rates = np.append(next_rates, possible_rate)

                            # store the values describing the state's mRNA, ribosome pool and asRNA pool
                            next_mrna = np.append(next_mrna, [possible_mrna],0)
                            next_ribo = np.append(next_ribo, possible_ribo)
                            next_asrna = np.append(next_asrna, possible_asrna)

                            # revert possible state to explore other possibilities from current state
                            possible_mrna = np.copy(current_mrna)
                            possible_asrna = current_asrna
                            possible_ribo = current_ribo

                            # start codon free
                            if current_mrna[1,n] == 0:

                                # ribosome moves from rbs to start codon
                                possible_mrna[0,n] = 0
                                possible_mrna[1,n] = 1
                                possible_ribo = possible_ribo
                                possible_asrna = possible_asrna

                                # rate of possible state change
                                possible_rate = start_transl

                                # store the rate of the new state transition
                                next_rates = np.append(next_rates, possible_rate)

                                # store the values describing the state's mRNA, ribosome pool and asRNA pool
                                next_mrna = np.append(next_mrna, [possible_mrna],0)
                                next_ribo = np.append(next_ribo, possible_ribo)
                                next_asrna = np.append(next_asrna, possible_asrna)

                                # revert possible state to explore other possibilities from current state
                                possible_mrna = np.copy(current_mrna)
                                possible_asrna = current_asrna
                                possible_ribo = current_ribo

                        # By asRNA
                        if current_mrna[0,n] == -1:

                            # asRNA unbind from RBS
                            possible_mrna[0,n] = 0
                            possible_asrna = possible_asrna + 1
                            possible_ribo = possible_ribo

                            # unbind from rbs: rate of possible state change
                            possible_rate = asrna_unbind_rate

                            # store the rate of the new state transition
                            next_rates = np.append(next_rates, possible_rate)

                            # store the values describing the state's mRNA, ribosome pool and asRNA pool
                            next_mrna = np.append(next_mrna, [possible_mrna],0)
                            next_ribo = np.append(next_ribo, possible_ribo)
                            next_asrna = np.append(next_asrna, possible_asrna)

                            # revert possible state to explore other possibilities from current state
                            possible_mrna = np.copy(current_mrna)
                            possible_asrna = current_asrna
                            possible_ribo = current_ribo

                    # codon occupied by ribosome
                    if i >= 1 and current_mrna[i,n] == 1:

                        # all codons except stop, next codon empty
                        if i < np.shape(current_mrna)[0]-1 and current_mrna[i + 1,n] == 0:

                            # ribosome moves to next codon
                            possible_mrna[i,n] = 0
                            possible_mrna[i + 1,n] = 1
                            possible_asrna = possible_asrna
                            possible_ribo = possible_ribo

                            # next codon is empty: rate of possible state change
                            possible_rate = codon_transl_rates[i-1]

                            # store the rate of the new state transition
                            next_rates = np.append(next_rates, possible_rate)

                            # store the values describing the state's mRNA, ribosome pool and asRNA pool
                            next_mrna = np.append(next_mrna, [possible_mrna],0)
                            next_ribo = np.append(next_ribo, possible_ribo)
                            next_asrna = np.append(next_asrna, possible_asrna)

                            # revert possible state to explore other possibilities from current state
                            possible_mrna = np.copy(current_mrna)
                            possible_asrna = current_asrna
                            possible_ribo = current_ribo

                        # on stop codon
                        if i == np.shape(current_mrna)[0]-1:

                            # stop codon reached, release rbs: possible state change
                            possible_mrna[i,n] = 0
                            possible_ribo = possible_ribo + 1
                            possible_asrna = possible_asrna

                            #rate of possible state change
                            possible_rate = fin_transl

                            # store the rate of the new state transition
                            next_rates = np.append(next_rates, possible_rate)

                            # store the values describing the state's mRNA, ribosome pool and asRNA pool
                            next_mrna = np.append(next_mrna, [possible_mrna],0)
                            next_ribo = np.append(next_ribo, possible_ribo)
                            next_asrna = np.append(next_asrna, possible_asrna)

                            # revert possible state to explore other possibilities from current state
                            possible_mrna = np.copy(current_mrna)
                            possible_asrna = current_asrna
                            possible_ribo = current_ribo

            # remove the first junk entry for the list of next mRNA states
            next_mrna = np.copy(next_mrna[1:])

            # total transition rate to other states
            esc_total = np.sum(next_rates)

            # Determine when next reaction occurs (ie. time interval):
            dt = float(- (math.log(1 - random.random())) / esc_total)

            # next event occurs after the end of simulation time
            if (current_time + dt) > simulation_time:

                current_time = simulation_time
                current_mrna = np.copy(current_mrna)
                current_ribo = current_ribo
                current_asrna = current_asrna

                trace_mrna = np.append(trace_mrna, [current_mrna],0)
                trace_ribo = np.append(trace_ribo, current_ribo)
                trace_asrna = np.append(trace_asrna, current_asrna)
                trace_time = np.append(trace_time, current_time)

            # next event occurs within the simulation time
            else:

                # Determine what reaction occurs next:
                # initialize random number threshold and scale to total transition rate
                threshold = random.random() * esc_total

                sum_trans_rates = 0

                for n in range(len(next_rates)):

                    # perform summation
                    sum_trans_rates = sum_trans_rates + next_rates[n]

                    # larger rates have proportionally higher chance to push sum above
                    if sum_trans_rates > threshold:

                        # define next state
                        new_mrna = next_mrna[n]
                        new_ribo = next_ribo[n]
                        new_asrna = next_asrna[n]
                        break

                current_time = current_time + dt
                current_mrna = np.copy(new_mrna)
                current_ribo = new_ribo
                current_asrna = new_asrna

                trace_mrna = np.append(trace_mrna, [current_mrna],0)
                trace_ribo = np.append(trace_ribo, current_ribo)
                trace_asrna = np.append(trace_asrna, current_asrna)
                trace_time = np.append(trace_time, current_time)

        return np.copy(trace_mrna), np.copy(trace_ribo), np.copy(trace_asrna), np.copy(trace_time), np.copy(
                    current_mrna), current_ribo, current_asrna, current_time, np.array(
                    [k_rib_bind, rib_unbind_rate, k_asrna_bind, asrna_unbind_rate])

    #####################################################
    # RUNNING SIMULATION

    # simulate until steady state is reached
    (trace_mrna, trace_ribo, trace_asrna, trace_time, steady_mrna, steady_ribo, steady_asrna, steady_time,
     kinetic_rates) = gillespie_multmR(initial_mrna, free_ribo, free_asrna, 0, window_start_time, start_transl, fin_transl, codon_transl_rates)

    # simulate from steady state onwards
    (steadytrace_mrna, steadytrace_ribo, steadytrace_asrna, steadytrace_time, final_mrna, final_ribo, final_asrna,
     final_time, steadykinetic_rates) = gillespie_multmR(steady_mrna, steady_ribo, steady_asrna, steady_time,
                                                          window_end_time, start_transl, fin_transl, codon_transl_rates)

    #####################################################
    # EXTRACTING PLOTS OF OCCUPANCY DISTRIBUTION

    # vector of all positions along mRNA
    rna_sites = np.arange(1, len(final_mrna) + 1, 1)

    # matrices to track total time sites on strands are occupied
    ribo_occ = np.zeros(np.shape(final_mrna))
    asrna_occ = np.zeros(np.shape(final_mrna))

    # mRNA site occupancy by ribosome
    # iterate over all states
    for j in range(1, len(steadytrace_time)):

        # time spent in state
        time_in_state = steadytrace_time[j] - steadytrace_time[j - 1]

        #iterate over all strands
        for i in range(0, num_strands):

            # extract each strand's occupancy in the state
            state_of_strand = steadytrace_mrna[j,:,i]

            # record time that was spent with rbs/codons on the rna occupied
            for k in range(len(state_of_strand)):

                if state_of_strand[k] == 1:
                    ribo_occ[k,i] = ribo_occ[k,i] + time_in_state

                if state_of_strand[k] == -1:
                    asrna_occ[k,i] = asrna_occ[k,i] + time_in_state

    avg_ribo_occ = np.mean(ribo_occ,axis=1)
    avg_asrna_occ = np.mean(asrna_occ,axis=1)

    std_asrna_occ = np.std(asrna_occ,axis=1)
    std_ribo_occ = np.std(ribo_occ,axis=1)

    if num_strands == 1:

        avg_ribo_occ = np.copy(ribo_occ)
        avg_asrna_occ = np.copy(asrna_occ)

    # plot time spent at each spot on the rna
    plt.figure()
    plt.subplot(2, 1, 1)
    bar_width = 1 / 1.5
    plt.bar(rna_sites, avg_ribo_occ, bar_width,yerr=std_ribo_occ,ecolor='r')
    plt.xlabel('Position along mRNA')
    plt.ylabel('Average time spent occupied (s)')
    plt.title('Number of mRNA strands: %d \n Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Mean ribosome occupancy along mRNA (Standard Deviation error)' % (num_strands,
        steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2], steadykinetic_rates[3]))

    plt.subplot(2, 1, 2)
    plt.bar(rna_sites, avg_asrna_occ, bar_width, yerr=std_asrna_occ, ecolor='r')
    plt.xlabel('Position along mRNA')
    plt.ylabel('Time spent occupied (s)')
    plt.title('Mean asRNA occupancy along mRNA (Standard Deviation error)')

    # save figure
    plt.tight_layout()
    plt.savefig('%sOccupancy distribution,num strands %d, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' % (
    save_dir,num_strands, steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2],
    steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure

    ######################################################
    # EXTRACTING PLOTS OF FREE POOLS

    # plot number of free ribosomes over time window
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.step(steadytrace_time, steadytrace_ribo, hold=True)
    plt.xlabel('Time (s)')
    plt.ylabel('Number of Free Ribosomes')
    plt.title('Number of mRNA strands: %d \n Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Free Ribosomes over time' %(num_strands,steadykinetic_rates[0],steadykinetic_rates[1],steadykinetic_rates[2],steadykinetic_rates[3]))

    # plot number of free asRNA over time window
    plt.subplot(2, 1, 2)
    plt.step(steadytrace_time,steadytrace_asrna,hold=True)
    plt.xlabel('Time (s)')
    plt.ylabel('Number of Free asRNA')
    plt.title('Free asRNA over time')

    # save figure
    plt.tight_layout()
    plt.savefig('%sFree pools, num strands %d, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' %(save_dir,num_strands,steadykinetic_rates[0],steadykinetic_rates[1],steadykinetic_rates[2],steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure

    #####################################################
    # EXTRACTING COMPLETED TRANSLATIONS OVER TIME

    # keeps track of number of completed translations of each strand at each time stamp
    count_transl = np.zeros(num_strands)

    # keeps track of total number of completed translations
    total_count_transl = np.zeros(1)

    # iterate over all states the mrna pool is in, starting from the 2nd state:
    for j in range(1, len(steadytrace_time)):

        # for each state iterate over all mrna strands
        for k in range(0,num_strands):

            # extract the strand's stop codon
            state_stop_codon = steadytrace_mrna[j,length_mrna-1,k]

            # extract previous state's stop codon
            prev_state_stop_codon = steadytrace_mrna[j-1,length_mrna-1,k]

            # if a ribosome has left the codon
            if state_stop_codon == 0 and prev_state_stop_codon == 1:

                count_transl[k] = count_transl[k] + 1

        # sum total translations completed by all strands together at each time stamp
        total_count_transl = np.append(total_count_transl, np.sum(count_transl))

    # plot number of completed translations over time
    plt.figure()
    plt.step(steadytrace_time, total_count_transl)
    plt.ylabel('Total number of completed mRNA translations')
    plt.xlabel('Time (s)')
    plt.title('Number of mRNA strands: %d \n Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Translations completed over time interval' % (num_strands,
        steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2], steadykinetic_rates[3]))

    # save figure
    plt.tight_layout()
    plt.savefig('%sTotal Completed Translations, num strands %d, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' % (
    save_dir, num_strands, steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2],
    steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure

    # plot distribution of translations completed by each mrna
    plt.figure()
    plt.hist(count_transl, bins = 125)
    plt.ylabel('Number of mRNA strands')
    plt.xlabel('Completed translations')
    plt.title('Total number of mRNA: %d \n Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Distribution of completed translations' % (num_strands,
        steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2], steadykinetic_rates[3]))

    # save figure
    plt.tight_layout()
    plt.savefig('%sDistribution of Completed Translations, num strands %d, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' % (
        save_dir,num_strands, steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2],
        steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure

    #####################################################
    # EXTRACTING AVERAGE SEQUESTERED RIBOSOMES ON MRNA

    # initialize count
    avg_seq_ribo = np.zeros(num_strands)

    # iterate over all states
    for j in range(1, len(steadytrace_time)):

        # time spent in state
        time_in_state = steadytrace_time[j] - steadytrace_time[j - 1]

        # iterate over all strands in each state
        for n in range(0,num_strands):

            # extract the each strand
            strand_of_rna = steadytrace_mrna[j,:,n]

            # initialize count of total ribosomes for the state
            count_ribo = 0

            # loop over strand of rna
            for k in range(0,len(strand_of_rna)):

                # total the number of ribosomes on the strand
                if strand_of_rna[k] == 1:
                    count_ribo = count_ribo + 1

            # calculate number of counted ribosomes weighted by time spent in the state during the simulation
            count_ribo = count_ribo * (time_in_state / (steadytrace_time[len(steadytrace_time) - 1] - steadytrace_time[0]))

            avg_seq_ribo[n] = avg_seq_ribo[n] + count_ribo

    # plot distribution of translations completed by each mrna
    plt.figure()
    plt.hist(avg_seq_ribo, bins=125)
    plt.ylabel('Number of mRNA strands')
    plt.xlabel('Average Sequestered ribosomes')
    plt.title(
        'Total Number of mRNA: %d \n Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Distribution of Sequestered Ribosomes' % (num_strands,
            steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2], steadykinetic_rates[3]))

    # save figure
    plt.tight_layout()
    plt.savefig(
        '%sDist. Sequestered Ribosomes, num strands %d, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' % (
            save_dir, num_strands, steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2],
            steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure