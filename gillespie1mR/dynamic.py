import numpy as np
import random, math

def dyn_pool(mrna_state, free_ribo_pool, free_asrna_pool, start_time, sim_time,k_rib_bind, rib_unbind_rate,k_asrna_bind,asrna_unbind_rate):

    # SETTING UP SIMULATION AND MODEL PARAMETERS

    # array of values describe the state of mRNA codons + RBS
    current_mrna = mrna_state

    # Number of ribosomes:
    current_ribo = free_ribo_pool

    # Number of asRNA:
    current_asrna = free_asrna_pool

    # Initial time:
    current_time = start_time

    # Simulate until time:
    simulation_time = sim_time

    # Translation initiation, elongation, completion rates:
    start_transl = 1
    transl_rate = 1
    fin_transl = 0.01

    # Array of values for translation elongation rates of specific codons on mRNA, not including stop codon or RBS:
    codon_transl_rates = np.full((len(current_mrna) - 2), transl_rate)

    # Set specific rate for specific codon(s), in this case 50th codon:
    codon_transl_rates[49] = 1

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

        # initialize variable to explore possible state changes of asRNA pool
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
                    next_rates = np.append(next_rates, possible_rate)

                    # store the 101 values describing state's mRNA
                    next_mrna = np.append(next_mrna, possible_mrna)

                    # store the value describing state's ribosome pool
                    next_ribo = np.append(next_ribo, possible_ribo)

                    # store the value describing state's asRNA pool
                    next_asrna = np.append(next_asrna, possible_asrna)

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
                if i < len(current_mrna) - 1 and current_mrna[i + 1] == 0:
                    # ribosome moves to next codon
                    possible_mrna[i] = 0
                    possible_mrna[i + 1] = 1
                    possible_asrna = possible_asrna
                    possible_ribo = possible_ribo

                    # next codon is empty: rate of possible state change
                    possible_rate = codon_transl_rates[i - 1]

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
                if i == len(current_mrna) - 1:
                    # stop codon reached, release rbs: possible state change
                    possible_mrna[i] = 0
                    possible_ribo = possible_ribo + 1
                    possible_asrna = possible_asrna

                    # rate of possible state change
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
            current_mrna = current_mrna
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
                    new_mrna = next_mrna[len(current_mrna) * n:len(current_mrna) * (n + 1)]
                    new_ribo = next_ribo[n]
                    new_asrna = next_asrna[n]
                    break

            current_time = current_time + dt
            current_mrna = new_mrna
            current_ribo = new_ribo
            current_asrna = new_asrna

            trace_mrna = np.append(trace_mrna, current_mrna)
            trace_ribo = np.append(trace_ribo, current_ribo)
            trace_asrna = np.append(trace_asrna, current_asrna)
            trace_time = np.append(trace_time, current_time)

    return np.copy(trace_mrna), np.copy(trace_ribo), np.copy(trace_asrna), np.copy(trace_time), np.copy(
        current_mrna), current_ribo, current_asrna, current_time, np.array(
        [k_rib_bind, rib_unbind_rate, k_asrna_bind, asrna_unbind_rate])