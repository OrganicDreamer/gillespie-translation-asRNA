import numpy as np
import random, math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

#####################################
#INITIALIZE ALL PARAMETERS HERE#####

# set specific directory to save output graphs to if desired (remember to include final '/' in directory path
save_dir = ''

# 101 values describe the state of mRNA codons + RBS
length_mrna = 267
initial_mrna = np.zeros(length_mrna)

free_ribo = 25
free_asrna = 25

# time to presumably reach steady state:
window_start_time = 2000

# final time to simulate until:
window_end_time = 3000

# 2d array (n x 4) of kinetic constants where each row is: ribo bind, ribo unbind, asRNA bind, asRNA unbind:
sweep_kinetic_const = np.array(
    [

    [1,1,0,0]

    ])

# Translation initiation, elongation, completion rates
start_transl = 0.0833333333333333
transl_rate = 15
fin_transl = range(transl_rate/500,(50*transl_rate)+transl_rate/10,transl_rate/10)

# Array of values for translation elongation rates of specific codons on mRNA, not including stop codon or RBS:
codon_transl_rates = np.full((length_mrna - 2), transl_rate)

# Set specific rate for specific codon(s), in this case 50th codon:
codon_transl_rates[49] = transl_rate

##########################################################################

# sweep over different rows of kinetic constants when simulating model:
for i in range(sweep_kinetic_const.shape[0]):

    # array to store steady state queue length over each iteration:
    plot_queue_len = np.array([])

    # store average flux of ribosomes in/out of queue at steady state:
    plot_ribo_flux = np.array([])

    # store deterministic prediction of ribosome flux in/out at steady state:
    plot_det_flux = np.array([])

    # parameter sweep slow codon rate:
    for n in range(len(fin_transl)):


        def gillespie_1mR_dyn(mrna_state,free_ribo_pool,free_asrna_pool,start_time,sim_time, start_transl, fin_transl, codon_transl_rates,i = i):

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

                free_rib_bind_rate = k_rib_bind * current_ribo
                free_asrna_bind_rate = k_asrna_bind * current_asrna

                # initialize vector to record possible next states' mRNA
                next_mrna = np.array([])

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

                # rules for next possible transitions and rates based on current state of mRNA:
                for i in range(len(current_mrna)):

                    # empty RBS
                    if (i == 0) and (current_mrna[i] == 0):

                        # Free ribosomes available to bind to RBS
                        if current_ribo >= 1:

                            # Free ribosome binds to RBS
                            possible_mrna[0] = 1
                            possible_ribo = possible_ribo - 1
                            possible_asrna = possible_asrna

                            # rate of the possible state change
                            possible_rate = free_rib_bind_rate

                            # store the rate of the new state transition
                            next_rates = np.append(next_rates,possible_rate)

                            # store the 101 values describing state's mRNA
                            next_mrna = np.append(next_mrna,possible_mrna)

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
                            possible_mrna[0] = -1

                            # Rate of possible state change
                            possible_rate = free_asrna_bind_rate

                            # store the rate of the new state transition
                            next_rates = np.append(next_rates, possible_rate)

                            # store the 101 values describing the state's mRNA
                            next_mrna = np.append(next_mrna, possible_mrna)

                            # store the value describing state's ribosome pool
                            next_ribo = np.append(next_ribo, possible_ribo)

                            # store the value describing state's asRNA pool
                            next_asrna = np.append(next_asrna, possible_asrna)

                            # revert possible state to explore other possibilities from current state
                            possible_mrna = np.copy(current_mrna)
                            possible_asrna = current_asrna
                            possible_ribo = current_ribo

                    # RBS occupied
                    if (i == 0) and current_mrna[0] != 0:

                        # By ribosome
                        if current_mrna[0] == 1:

                            # ribosome unbind from RBS
                            possible_mrna[0] = 0
                            possible_ribo = possible_ribo + 1
                            possible_asrna = possible_asrna

                            # unbind from rbs: rate of possible state change
                            possible_rate = rib_unbind_rate

                            # store the rate of the new state transition
                            next_rates = np.append(next_rates, possible_rate)

                            # store the values describing the state's mRNA, ribosome pool and asRNA pool
                            next_mrna = np.append(next_mrna, possible_mrna)
                            next_ribo = np.append(next_ribo, possible_ribo)
                            next_asrna = np.append(next_asrna, possible_asrna)

                            # revert possible state to explore other possibilities from current state
                            possible_mrna = np.copy(current_mrna)
                            possible_asrna = current_asrna
                            possible_ribo = current_ribo

                            # start codon free
                            if current_mrna[1] == 0:

                                # ribosome moves from rbs to start codon
                                possible_mrna[0] = 0
                                possible_mrna[1] = 1
                                possible_ribo = possible_ribo
                                possible_asrna = possible_asrna

                                # rate of possible state change
                                possible_rate = start_transl

                                # store the rate of the new state transition
                                next_rates = np.append(next_rates, possible_rate)

                                # store the values describing the state's mRNA, ribosome pool and asRNA pool
                                next_mrna = np.append(next_mrna, possible_mrna)
                                next_ribo = np.append(next_ribo, possible_ribo)
                                next_asrna = np.append(next_asrna, possible_asrna)

                                # revert possible state to explore other possibilities from current state
                                possible_mrna = np.copy(current_mrna)
                                possible_asrna = current_asrna
                                possible_ribo = current_ribo

                        # By asRNA
                        if current_mrna[0] == -1:

                            # asRNA unbind from RBS
                            possible_mrna[0] = 0
                            possible_asrna = possible_asrna + 1
                            possible_ribo = possible_ribo

                            # unbind from rbs: rate of possible state change
                            possible_rate = asrna_unbind_rate

                            # store the rate of the new state transition
                            next_rates = np.append(next_rates, possible_rate)

                            # store the values describing the state's mRNA, ribosome pool and asRNA pool
                            next_mrna = np.append(next_mrna, possible_mrna)
                            next_ribo = np.append(next_ribo, possible_ribo)
                            next_asrna = np.append(next_asrna, possible_asrna)

                            # revert possible state to explore other possibilities from current state
                            possible_mrna = np.copy(current_mrna)
                            possible_asrna = current_asrna
                            possible_ribo = current_ribo

                    # codon occupied by ribosome
                    if i >= 1 and current_mrna[i] == 1:

                        # all codons except stop, next codon empty
                        if i < len(current_mrna)-1 and current_mrna[i + 1] == 0:

                            # ribosome moves to next codon
                            possible_mrna[i] = 0
                            possible_mrna[i + 1] = 1
                            possible_asrna = possible_asrna
                            possible_ribo = possible_ribo

                            # next codon is empty: rate of possible state change
                            possible_rate = codon_transl_rates[i-1]

                            # store the rate of the new state transition
                            next_rates = np.append(next_rates, possible_rate)

                            # store the values describing the state's mRNA, ribosome pool and asRNA pool
                            next_mrna = np.append(next_mrna, possible_mrna)
                            next_ribo = np.append(next_ribo, possible_ribo)
                            next_asrna = np.append(next_asrna, possible_asrna)

                            # revert possible state to explore other possibilities from current state
                            possible_mrna = np.copy(current_mrna)
                            possible_asrna = current_asrna
                            possible_ribo = current_ribo

                        # on stop codon
                        if i == len(current_mrna)-1:

                            # stop codon reached, release rbs: possible state change
                            possible_mrna[i] = 0
                            possible_ribo = possible_ribo + 1
                            possible_asrna = possible_asrna

                            #rate of possible state change
                            possible_rate = fin_transl

                            # store the rate of the new state transition
                            next_rates = np.append(next_rates, possible_rate)

                            # store the values describing the state's mRNA, ribosome pool and asRNA pool
                            next_mrna = np.append(next_mrna, possible_mrna)
                            next_ribo = np.append(next_ribo, possible_ribo)
                            next_asrna = np.append(next_asrna, possible_asrna)

                            # revert possible state to explore other possibilities from current state
                            possible_mrna = np.copy(current_mrna)
                            possible_asrna = current_asrna
                            possible_ribo = current_ribo

                # total transition rate to other states
                esc_total = np.sum(next_rates)

                # Determine when next reaction occurs (ie. time interval):
                dt = - (math.log(1 - random.random())) / esc_total

                # next event occurs after the end of simulation time
                if (current_time + dt) > simulation_time:

                    current_time = simulation_time
                    current_mrna = np.copy(current_mrna)
                    current_ribo = current_ribo
                    current_asrna = current_asrna

                    trace_mrna = np.append(trace_mrna, current_mrna)
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
                            new_mrna = next_mrna[len(current_mrna)*n:len(current_mrna)*(n+1)]
                            new_ribo = next_ribo[n]
                            new_asrna = next_asrna[n]
                            break

                    current_time = current_time + dt
                    current_mrna = np.copy(new_mrna)
                    current_ribo = new_ribo
                    current_asrna = new_asrna

                    trace_mrna = np.append(trace_mrna, current_mrna)
                    trace_ribo = np.append(trace_ribo, current_ribo)
                    trace_asrna = np.append(trace_asrna, current_asrna)
                    trace_time = np.append(trace_time, current_time)

            return np.copy(trace_mrna),np.copy(trace_ribo),np.copy(trace_asrna),np.copy(trace_time), np.copy(current_mrna),current_ribo,current_asrna,current_time,np.array([k_rib_bind,rib_unbind_rate,k_asrna_bind,asrna_unbind_rate])

        #####################################################
        # RUNNING SIMULATION

        # simulate until steady state is reached
        (trace_mrna, trace_ribo, trace_asrna,trace_time, steady_mrna, steady_ribo, steady_asrna, steady_time,kinetic_rates) = gillespie_1mR_dyn(initial_mrna,free_ribo,free_asrna,0,window_start_time, start_transl, fin_transl, codon_transl_rates)

        # simulate from steady state onwards
        (steadytrace_mrna, steadytrace_ribo,steadytrace_asrna, steadytrace_time,final_mrna, final_ribo,final_asrna,final_time,steadykinetic_rates) = gillespie_1mR_dyn(steady_mrna,steady_ribo,steady_asrna,steady_time, window_end_time, start_transl, fin_transl, codon_transl_rates)

        #####################################################

        # EXTRACTING PLOTS OF OCCUPANCY DISTRIBUTION

        # vector to keep track to total time site on mRNA is occupied
        ribo_occ = np.zeros(np.shape(final_mrna))
        asrna_occ = np.zeros(np.shape(final_mrna))

        # vector of all positions along mRNA
        rna_loc = np.arange(1, len(final_mrna) + 1, 1)

        # mRNA site occupancy by ribosome
        #iterate over all states starting from 2nd state
        for j in range(1, len(steadytrace_time)):

            # time spent in previous state
            time_in_state = steadytrace_time[j] - steadytrace_time[j - 1]

            # extract the previous state's rna occupancy
            state_of_rna = steadytrace_mrna[((j - 1) * len(final_mrna)): j * len(final_mrna)]

            # record time that was spent with rbs/codons on the rna occupied
            for k in range(len(state_of_rna)):

                if state_of_rna[k] == 1:
                    ribo_occ[k] = ribo_occ[k] + time_in_state

                if state_of_rna[k] == -1:
                    asrna_occ[k] = asrna_occ[k] + time_in_state

        # plot time spent at each spot on the rna
        plt.figure()
        plt.subplot(2,1,1)
        bar_width = 1 / 1.5
        plt.bar(rna_loc, ribo_occ, bar_width, hold=True)
        plt.xlabel('Position along mRNA')
        plt.ylabel('Time spent occupied (s)')
        plt.title('Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Ribosome occupancy along mRNA' %(steadykinetic_rates[0],steadykinetic_rates[1],steadykinetic_rates[2],steadykinetic_rates[3]))

        # plot time spent at each spot on the rna
        plt.subplot(2,1,2)
        bar_width = 0.75
        plt.bar(rna_loc, asrna_occ, bar_width, hold=True)
        plt.xlabel('Position along mRNA')
        plt.ylabel('Time spent occupied (s)')
        plt.title('asRNA occupancy along mRNA')

        # save figure
        plt.tight_layout()
        plt.savefig('%sOccupancy distribution, finite, %fs, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' %(save_dir,window_end_time,steadykinetic_rates[0],steadykinetic_rates[1],steadykinetic_rates[2],steadykinetic_rates[3]))  # save the figure to file
        plt.close()  # close the figure
        ######################################################

        # EXTRACTIONG PLOTS OF FREE POOLS

        # plot number of free ribosomes over time window
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.step(steadytrace_time, steadytrace_ribo, hold=True)
        plt.xlabel('Time (s)')
        plt.ylabel('Number of Free Ribosomes')
        plt.title('Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Free Ribosomes over time' %(steadykinetic_rates[0],steadykinetic_rates[1],steadykinetic_rates[2],steadykinetic_rates[3]))

        # plot number of free asRNA over time window
        plt.subplot(2, 1, 2)
        plt.step(steadytrace_time,steadytrace_asrna,hold=True)
        plt.xlabel('Time (s)')
        plt.ylabel('Number of Free asRNA')
        plt.title('Free asRNA over time')

        # save figure
        plt.tight_layout()
        plt.savefig('%sFree pools, finite, %fs, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' %(save_dir,window_end_time,steadykinetic_rates[0],steadykinetic_rates[1],steadykinetic_rates[2],steadykinetic_rates[3]))  # save the figure to file
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
        plt.savefig('%sCompleted Translations, finite, %fs, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' %(save_dir,window_end_time,steadykinetic_rates[0],steadykinetic_rates[1],steadykinetic_rates[2],steadykinetic_rates[3]))  # save the figure to file
        plt.close()  # close the figure

        #####################################################
        # CALCULATE QUEUE LENGTH

        # find number of sequestered ribosomes on mRNA
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
        ribo_flux = (count_transl[len(count_transl) - 1] - count_transl[0]) / (window_end_time - window_start_time)

        # store result for this iteration
        plot_ribo_flux = np.append(plot_ribo_flux, ribo_flux)

        # deterministic ribosome flux out:
        det_flux = fin_transl * (transl_rate / (fin_transl + transl_rate))

        # store result for this iteration
        plot_det_flux = np.append(plot_det_flux, det_flux)

    ######################################################################
    # GENERATE AND SAVE PARAMETER SWEEP PLOTS
    plt.figure()
    plt.plot(fin_transl, plot_queue_len)
    plt.ylabel('Number of ribosomes in queue')
    plt.xlabel('Slow codon translation rate')
    plt.title('Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Queue length vs. slow codon rate' % (
        steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2], steadykinetic_rates[3]))

    plt.tight_layout()
    plt.savefig(
        '%sQueue length sweep, fin, %fs, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' % (
            save_dir, window_end_time, steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2],
            steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure

    plt.figure()
    plt.plot(fin_transl, plot_ribo_flux, 'b', label='Simulated flux', hold=True)
    plt.plot(fin_transl, plot_det_flux, 'r', label='Deterministic flux', hold=True)
    plt.legend(loc='upper right')
    plt.ylabel('Ribosome flux out')
    plt.xlabel('Slow codon translation rate')
    plt.title('Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Flux vs. slow codon rate' % (
        steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2], steadykinetic_rates[3]))

    plt.tight_layout()
    plt.savefig(
        '%sFlux sweep, fin, %fs, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' % (
            save_dir, window_end_time, steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2],
            steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure