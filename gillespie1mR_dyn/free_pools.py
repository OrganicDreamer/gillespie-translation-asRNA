import matplotlib.pyplot as plt

def free_pools(steadytrace_time,steadytrace_ribo,steadytrace_asrna,steadykinetic_rates):

    # EXTRACTIONG PLOTS OF FREE POOLS

    # plot number of free ribosomes over time window
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.step(steadytrace_time, steadytrace_ribo, hold=True)
    plt.xlabel('Time (s)')
    plt.ylabel('Number of Free Ribosomes')
    plt.title('Ribo bind: %f, unbind: %f \n asRNA bind: %f, unbind: %f \n Free Ribosomes over time' % (
    steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2], steadykinetic_rates[3]))

    # plot number of free asRNA over time window
    plt.subplot(2, 1, 2)
    plt.step(steadytrace_time, steadytrace_asrna, hold=True)
    plt.xlabel('Time (s)')
    plt.ylabel('Number of Free asRNA')
    plt.title('Free asRNA over time')

    # save figure
    plt.tight_layout()
    plt.savefig(
        'C:/Users/fujitsu/Desktop/Y4 Project/1mR_dyn/Free pools, ribo bind %f, ribo unbind %f, asR bind %f, asR unbind %f.png' % (
        steadykinetic_rates[0], steadykinetic_rates[1], steadykinetic_rates[2],
        steadykinetic_rates[3]))  # save the figure to file
    plt.close()  # close the figure