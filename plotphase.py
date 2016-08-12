import numpy as np
import pylab as pl
from brian2 import *

def frate(st, window_width=5, dt=0.01):
    '''
    st : spike train
    window_width : gaussian kernel width in ms
    dt : time step in ms
    '''
    width_dt = int(window_width/dt)
    window = exp(-arange(-4 * width_dt, 4 * width_dt + 1) ** 2 * 1. / (2 * (width_dt) ** 2))
    rate = zeros(int(500./(dt)))
    if len(st)>0:
        st = [int(x/dt) for x in st]
        rate[st] = 1000./dt
    return convolve(rate, window * 1. / sum(window), mode='same')

C = 250*pF              # Membrane Capacitance (pF)
gL = 15*nS              # Leak Conductance (nS)
DeltaT = 2*mV           # Spike Slope Factor (mV)
taum = C/gL             # Membrane Time Constant (pF/nS=ms)
EL = -55*mV             # Leak Reversal Potential (mV)
VT = -50*mV             # Spike Threshold (mV)
Vr = -45*mV             # Reset potential (mV)
Vr2 = -20*mV
Vcut = 0*mV             # Spiking Threshold (mV)
a = 4*nS                # Subthreshold Adaptation (nS)
b = 60*pA              # Spike-triggered Adaptation (nA)
tauw = 120*ms            # Adaptation Time Constant (ms)

inputIcl = 700*pA
vms = np.linspace(-120*mV, Vcut, 2000)

vnull =gL*(EL-vms)+gL*DeltaT*exp((vms-VT)/DeltaT)+inputIcl #v-nullcline
vnullbase =gL*(EL-vms)+gL*DeltaT*exp((vms-VT)/DeltaT) #v-nullcline I =0
wnull =a*(vms-EL) #w-nullcline

eqs="""
dvm/dt=(gL*(EL-vm)+gL*DeltaT*exp((vm-VT)/DeltaT)+I-w)/C : volt
dw/dt=(a*(vm-EL)-w)*(1./tauw)                           : amp
I                                                       : amp
"""

reset="""
vm=Vr;
w+=b
"""

reset2="""
vm=Vr2;
"""

defaultclock.dt = 0.01*ms

def ppp():

    NV = 12
    Nw = 11

    Neuron = NeuronGroup(NV*Nw, model=eqs, threshold='vm>Vcut', reset=reset2, method='euler')
    Neurontra = NeuronGroup(1, model=eqs, threshold='vm>Vcut', reset=reset, method='euler')

    inputI = 0*amp
    Neuron.I=inputI
    Neurontra.I=inputI

    run(100*ms)

    M = StateMonitor(Neuron, ('vm', 'w'), record=True)
    s = SpikeMonitor(Neuron)

    Mtra = StateMonitor(Neurontra, ('vm', 'w'), record=True)
    stra = SpikeMonitor(Neurontra)

    inputI = inputIcl

    Neuron.I=inputI
    Neurontra.I=inputI

    VS = np.linspace(-80, -25, NV)
    WS = np.linspace(-0.5, 1.5, Nw)

    #defaultclock.t=0*second
    #run(100*ms, report='text', report_period=1*second)
    for i in xrange(NV*Nw):
        wi, vi = divmod(i, NV)
        Neuron.vm[i] = VS[vi]*mV
        Neuron.w[i] = WS[wi]*nA

    Neurontra.vm[0] = EL
    Neurontra.w[0] = 0*amp

    print "stimulation start"
    run(50*ms, report='text', report_period=1*second)
    print "stimulation over"

    inputI = 0*nA
    Neuron.I=inputI
    Neurontra.I=inputI

    run(150*ms, report='text', report_period=1*second)


    pl.figure()
    Nvis = 100

    for i in xrange(NV*Nw):
        wi, vi = divmod(i, NV)
        pl.plot(M.vm[i][1:Nvis]/mV, M.w[i][1:Nvis]/nA, 'b', ms=2, lw=1)

    pl.plot(vms*1000, vnullbase*1e9, 'k--', lw=1)
    pl.plot(vms*1000, vnull*1e9, 'k--', lw=1)
    pl.plot(vms*1000, wnull*1e9, 'k--', lw=1)

    pl.plot(Mtra.vm[0]/mV, Mtra.w[0]/nA, 'k.', ms=0.5)

    pl.xlim(-80, -30)
    pl.ylim(-0.5, 1.5)
    pl.xlabel("Membrane Potential (mV)")
    pl.ylabel("Adaptive current (nA)")
    
    
    pl.figure()
    xstart = 80
    xend = 200
    pl.subplot(311)
    pl.plot(Mtra.t/ms, Mtra.vm[0]/mV)
    pl.xlim(xstart, xend)
    pl.ylabel("Membrane potential (mV)")
    
    pl.subplot(312)
    pl.plot(Mtra.t/ms, Mtra.w[0]/nA)
    pl.xlim(xstart, xend)
    pl.ylabel("Adaptive Current (nA)")

    pl.subplot(313)
    strain = stra.spike_trains()
    firing = frate(strain[0]/ms)
    pl.plot(np.linspace(0, 500, len(firing)), firing)
    pl.xlim(xstart, xend)
    pl.xlabel("Time (ms)")
    pl.ylabel("Spike density (spks/s)")
    
    pl.show()
    print stra.num_spikes
