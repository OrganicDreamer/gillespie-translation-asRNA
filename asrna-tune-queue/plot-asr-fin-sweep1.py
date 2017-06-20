import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

as_on_mult_sweep = np.linspace(2,128,num=25)
num_as_on_rates = len(as_on_mult_sweep) + 1
results = np.empty([num_as_on_rates,4])

for i in range(num_as_on_rates):

    if i == 0:

        results_entry = np.load('fin-asrna.npz')
        results[i,0] = results_entry['as_bind']
        results[i,1] = results_entry['as_unbind']
        results[i,2] = results_entry['phi0']
        results[i,3] = results_entry['queue_len']

    else:

        results_entry = np.load('%dasr-on-table-entry.npz' %(as_on_mult_sweep[i-1]))
        results[i, 0] = as_on_mult_sweep[i-1]*results_entry['as_bind']
        results[i, 1] = results_entry['as_unbind']
        results[i, 2] = results_entry['ribo_flux']
        results[i, 3] = results_entry['queue_len']

plt.figure()
plt.scatter(results[:,0],results[:,2])
plt.ylabel('Ribosome flux out (ribosomes/(mRNA * s))')
plt.xlabel('asRNA binding rate (1/(asRNA mol * s))')
plt.title('Translation vs. asRNA binding')
plt.tight_layout()
plt.savefig(
'fin-final1.png')  # save the figure to file
plt.close()  # close the figure

plt.figure()
plt.scatter(results[:,0],results[:,3])
plt.ylabel('Queue length (ribosomes)')
plt.xlabel('asRNA binding rate (1/(asRNA mol * s))')
plt.title('Sequestration vs. asRNA binding')
plt.tight_layout()
plt.savefig(
'fin-final2.png')  # save the figure to file
plt.close()  # close the figure

plt.figure()
plt.scatter(results[:,1],results[:,2])
plt.ylabel('Ribosome flux out (ribosomes/(mRNA * s))')
plt.xlabel('asRNA unbinding rate (1/s)')
plt.title('Translation vs. asRNA unbinding')
plt.tight_layout()
plt.savefig(
'fin-final3.png')  # save the figure to file
plt.close()  # close the figure

plt.figure()
plt.scatter(results[:,1],results[:,3])
plt.ylabel('Queue length (ribosomes)')
plt.xlabel('asRNA unbinding rate (1/s)')
plt.title('Sequestration vs. asRNA unbinding')
plt.tight_layout()
plt.savefig(
'fin-final4.png')  # save the figure to file
plt.close()  # close the figure

plt.figure()
plt.scatter(results[:,1],results[:,0])
plt.ylabel('asRNA binding rate (1/(asRNA mol * s))')
plt.xlabel('asRNA unbinding rate (1/s)')
plt.title('At an equal or increased rate of protein production')
plt.tight_layout()
plt.savefig(
'fin-final5.png')  # save the figure to file
plt.close()  # close the figure