import numpy as np

as_on_mult_sweep = np.linspace(2,128,num=25)
baseline = np.load('fin-asrna.npz')
base_asunbind = baseline['as_unbind']
sweep_as_unbind = np.linspace(base_asunbind,65*base_asunbind,num=25)

params = np.empty([1,2])

for i in range(len(as_on_mult_sweep)):

    for j in range(len(sweep_as_unbind)):

        row = np.array([as_on_mult_sweep[i],sweep_as_unbind[j]])
        params=np.append(params,[row],axis=0)

params = params[1:len(params),:]

for x in range(len(params)):
    np.savez(str('%diter-params' %(x)),params=params[x])
