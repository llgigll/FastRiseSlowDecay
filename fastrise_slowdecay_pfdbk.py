# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 14:35:55 2014

This script shifts the synaptic depression into two equations for each neuron
rather than computing it for every synapse. This is possible since every
excitatory neuron only has two discrete possibilities for STD.

@author: Grant
"""
#from Spike_Stats import *
import brian as bn
import time
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import seaborn as sns
from scipy.ndimage import filters
from scipy.signal import gaussian
import pickle

def rep_supp(exp_struct, net_struct,new_connectivity,save_file):
  bn.seed(int(time.time()))
  bn.reinit_default_clock()
#  bn.seed(1412958308+2)
  bn.defaultclock.dt = 0.5*bn.ms
  
  #==============================================================================
  # Define constants for the model.
  #==============================================================================
  print net_struct
  print new_connectivity
#  normalize_connections = 'population'
  normalize_connections = 'neuron'
#  normalize_connections = ''
  
  stim1 = exp_struct['stim1']
  ISI1 = exp_struct['ISI1']
  set1_begin = exp_struct['set1_begin']
  set1_end = exp_struct['set1_end']
  
  T = set1_end
  
  delta_q = exp_struct['delta_q']
  delta_tr= exp_struct['delta_tr']
  delta_u = exp_struct['delta_u']
  
  u_base = exp_struct['u_base']
  trec_base = exp_struct['trec_base']
  
  uff = exp_struct['uff']
  trecff = exp_struct['trecff']
  Jin = exp_struct['Jin']
  
  neuron_mult = exp_struct['neuron_mult']
  baseline_mult = exp_struct['baseline_mult']
  k = exp_struct['k']
  w = exp_struct['w']
  
  ro = exp_struct['ro']

  
  qee1 = 0.75+delta_q # Fraction of NMDA receptors for e to e connections
  qee2 = 0.25+delta_q
  qie1 = 0.75-delta_q # Fraction of NMDA receptors for e to i connections
  qie2 = 0.25-delta_q
  uee1 = u_base-delta_u
  uee2 = u_base+delta_u
  trec1 = trec_base-delta_tr
  trec2 = trec_base+delta_tr
  
  Ne   = 3200*neuron_mult # number of excitatory neurons
  Ni   = 800*neuron_mult # number of inhibitory neurons
  Nout = 100 # number of neurons in the output group
  Nin  = 5000 # number of neurons in the input group
  N = Ne + Ni + Nout
  
  pcon = 0.2 # probability of connection
  pcon_out = 0.025
  pcon_in = 0.2
  
  Jee = w/(Ne*pcon)
  Jie = 0
  Jii = 0
  Jei = 0

  Jout = 0.2/(Ne*pcon_out)
  Jin = Jin/(Nin*pcon_in)
  
  El = -60.0*bn.mV # leak reversal potential
  Vreset = -52.0*bn.mV # reversal potential
  Vthresh = -40.0*bn.mV # spiking threshold
  
  tref = 2.0*bn.ms # refractory period
  te = 20.0*bn.ms # membrane time constant of excitatory neurons
  ti = 20.0*bn.ms # membrane time constant of inhibitory neruons
  tee_ampa = 10.0*bn.ms # time const of ampa currents at excitatory neurons
  tee_nmda = 200.0*bn.ms # time const of nmda currents at excitatory neurons
  tie_ampa = 10.0*bn.ms  # time const of ampa currents at inhibitory neurons
  tie_nmda = 200.0*bn.ms # time const of nmda currents at inhibitory neurons
  tii_gaba = 10.0*bn.ms # time const of GABA currents at inhibitory neurons
  tei_gaba = 10.0*bn.ms # time const of GABA currents at excitatory neurons
  teo_input = 100.0*bn.ms
  
  #==============================================================================
  # Define model structure
  #==============================================================================
  
  model = '''
  dV/dt = (-(V-El)+I_ampa+I_nmda-I_gaba+I_input+eta)/tm : bn.volt
  dI_ampa/dt = -I_ampa/t_ampa : bn.volt
  dI_nmda/dt = -I_nmda/t_nmda : bn.volt
  dI_gaba/dt = -I_gaba/t_gaba : bn.volt
  dI_input/dt = (-I_input+mu)/t_input : bn.volt
  dx1/dt = (1-x1)/t1_rec : 1
  dx2/dt = (1-x2)/t2_rec : 1
  u1 : 1
  t1_rec : bn.second
  u2 : 1
  t2_rec : bn.second
  mu : bn.volt
  eta : bn.volt
  J_ampa : 1
  J_nmda : 1
  J_gaba : 1
  J_input : 1
  tm : bn.second
  t_ampa : bn.second
  t_nmda : bn.second
  t_gaba : bn.second
  t_input : bn.second
  '''
  
  P_reset = "V=-52*bn.mV;x1+=-u1*x1;x2+=-u2*x2"
  
  Se_model = '''
  we_ampa1 : bn.volt
  we_nmda1 : bn.volt
  we_ampa2 : bn.volt
  we_nmda2 : bn.volt
  '''
  
  Se_pre = ('I_ampa += x1_pre*we_ampa1','I_nmda += x1_pre*we_nmda1',
            'I_ampa += x2_pre*we_ampa2','I_nmda += x2_pre*we_nmda2')
    
  Si_model = '''
  wi_gaba : bn.volt
  '''
  
  Si_pre = 'I_gaba += wi_gaba'  
  
  Sout_model = '''
  wout_ampa : bn.volt
  '''
  
  Sout_pre = 'I_ampa += wout_ampa'
  
  Sin_model ='''
  x : 1
  u : 1
  trec : bn.second
  win_ampa : bn.volt
  win_nmda : bn.volt
  win_input : bn.volt
  '''
  
  Sin_pre = '''
  x=1+(x-1)*bn.exp(-(t-lastupdate)/trec)
  I_ampa+=win_ampa*x
  I_nmda+=win_nmda*x
  I_input+=win_input*x
  x+=-u*x
  '''
  
  #==============================================================================
  # Define populations
  #==============================================================================
  
  P = bn.NeuronGroup(N, model, threshold = Vthresh, reset = P_reset, refractory = tref)
  
  Pe = P[0:Ne]
  Pe.tm = te
  Pe.t_ampa = tee_ampa
  Pe.t_nmda = tee_nmda
  Pe.t_gaba = tei_gaba
  Pe.t_input = teo_input
  Pe.I_ampa = 0*bn.mV
  Pe.I_nmda = 0*bn.mV
  Pe.I_gaba = 0*bn.mV
  Pe.I_input = 0*bn.mV
  Pe.V = (np.random.rand(Pe.V.size)*12-52)*bn.mV
  
  Pe.x1 = 1.0
  Pe.x2 = 1.0
  Pe.u1 = uee1
  Pe.u2 = uee2
  Pe.t1_rec = trec1
  Pe.t2_rec = trec2
  
  
  Pi = P[Ne:(Ne+Ni)]
  Pi.tm = ti
  Pi.t_ampa = tie_ampa
  Pi.t_nmda = tie_nmda
  Pi.t_gaba = tii_gaba
  Pi.t_input = teo_input
  Pi.I_ampa = 0*bn.mV
  Pi.I_nmda = 0*bn.mV
  Pi.I_gaba = 0*bn.mV
  Pi.I_input = 0*bn.mV
  Pi.V = (np.random.rand(Pi.V.size)*12-52)*bn.mV
  
  Pi.x1 = 1.0
  Pi.x2 = 1.0
  Pi.u1 = 0.0
  Pi.u2 = 0.0
  Pi.t1_rec = 1.0
  Pi.t2_rec = 1.0
  
  Pout = P[(Ne+Ni):]
  Pout.tm = te
  Pout.t_ampa = tee_ampa
  Pout.t_nmda = tee_nmda
  Pout.t_gaba = tii_gaba
  Pout.t_input = teo_input
  Pout.I_ampa = 0*bn.mV
  Pout.I_nmda = 0*bn.mV
  Pout.I_gaba = 0*bn.mV
  Pout.I_input = 0*bn.mV
  Pout.V = (np.random.rand(Pout.V.size)*12-52)*bn.mV
  
  Pout.x1 = 1.0
  Pout.x2 = 1.0
  Pout.u1 = 0.0
  Pout.u2 = 0.0
  Pout.t1_rec = 1.0
  Pout.t2_rec = 1.0
  
  #==============================================================================
  # Define inputs
  #==============================================================================

  Background_eo = bn.PoissonInput(Pe, N=5000, rate=5.0*baseline_mult*bn.Hz, weight=0.2*bn.mV, state='I_ampa')
#  Background_io = bn.PoissonInput(Pi, N=5000, rate=5.1*baseline_mult*bn.Hz, weight=0.2*bn.mV, state='I_ampa')
#  Background_out = bn.PoissonInput(Pout, N=1000, rate=0.9*bn.Hz, weight=0.2*bn.mV, state='I_ampa')

  def firing_function(t):
    if (t > set1_begin) & (t < set1_end):
        if ((t-set1_begin) % (stim1+ISI1)) <= stim1:
            return ro
        else:
            return 0.0 * bn.Hz
    else:
        return 0.0 * bn.Hz

  Pe.mu = 0*bn.mV
  Pe.eta = 0*bn.mV#, dt=0.5*bn.ms)  
  Pi.mu =  0.0*bn.mV
  Pi.eta = 0*bn.mV#, dt=0.5*bn.ms)

  Pin = bn.PoissonGroup(Nin,rates=lambda t:firing_function(t))  
#  Pin = bn.PoissonGroup(Nin,rates=0.5*bn.Hz)  

  #==============================================================================
  # Define synapses  
  #==============================================================================
  
  See1 = bn.Synapses(Pe, Pe, model = Se_model, pre = Se_pre)
  See2 = bn.Synapses(Pe, Pe, model = Se_model, pre = Se_pre)
  Sie1 = bn.Synapses(Pe, Pi, model = Se_model, pre = Se_pre)
  Sie2 = bn.Synapses(Pe, Pi, model = Se_model, pre = Se_pre)
  
  Sei = bn.Synapses(Pi, Pe, model = Si_model, pre = Si_pre)
  Sii = bn.Synapses(Pi, Pi, model = Si_model, pre = Si_pre)
  
  Sin = bn.Synapses(Pin, Pe, model = Sin_model, pre = Sin_pre)  
  Sout = bn.Synapses(Pe, Pout, model = Sout_model, pre = Sout_pre)
  
  #==============================================================================
  # Define random connections
  #==============================================================================
  
  
  if new_connectivity:
    See1.connect_random(Pe,Pe,sparseness=pcon/2.0)
    See2.connect_random(Pe,Pe,sparseness=pcon/2.0)
    Sie1.connect_random(Pe,Pi,sparseness=pcon/2.0)  
    Sie2.connect_random(Pe,Pi,sparseness=pcon/2.0)
    Sii.connect_random(Pi,Pi,sparseness=pcon)
    Sei.connect_random(Pi,Pe,sparseness=pcon)
    Sin.connect_random(Pin,Pe,sparseness=pcon_in)
    Sout.connect_random(Pe,Pout,sparseness=pcon_out)
  
    print 'Saving'
    See1.save_connectivity('./See1_connections_fastrise'+str(net_struct))
    See2.save_connectivity('./See2_connections_fastrise'+str(net_struct))
    Sie1.save_connectivity('./Sie1_connections_fastrise'+str(net_struct))
    Sie2.save_connectivity('./Sie2_connections_fastrise'+str(net_struct))
    Sii.save_connectivity('./Sii_connections_fastrise'+str(net_struct))
    Sei.save_connectivity('./Sei_connections_fastrise'+str(net_struct))
    Sin.save_connectivity('./Sin_connections_fastrise'+str(net_struct))
    Sout.save_connectivity('./Sout_connections_fastrise'+str(net_struct))

  else:
    print 'Loading'
    See1.load_connectivity('./See1_connections_fastrise'+str(net_struct))
    See2.load_connectivity('./See2_connections_fastrise'+str(net_struct))
    Sie1.load_connectivity('./Sie1_connections_fastrise'+str(net_struct))
    Sie2.load_connectivity('./Sie2_connections_fastrise'+str(net_struct))
    Sii.load_connectivity('./Sii_connections_fastrise'+str(net_struct))
    Sei.load_connectivity('./Sei_connections_fastrise'+str(net_struct))
    Sin.load_connectivity('./Sin_connections_fastrise'+str(net_struct))
    Sout.load_connectivity('./Sout_connections_fastrise'+str(net_struct))
  
  
  if normalize_connections == 'population':
      See1_norm = (Ne*Ne*(pcon/2))/len(See1.postsynaptic)
      See2_norm = (Ne*Ne*(pcon/2))/len(See2.postsynaptic)
      Sie1_norm = (Ni*Ne*(pcon/2))/len(Sie1.postsynaptic)
      Sie2_norm = (Ni*Ne*(pcon/2))/len(Sie2.postsynaptic)
      Sii_norm  = (Ni*Ni*pcon)/len(Sii.postsynaptic)
      Sei_norm  = (Ne*Ni*pcon)/len(Sei.postsynaptic)
      
      Sin_norm  = (Nin*Ne*pcon_in)/len(Sin.postsynaptic)
      Sout_norm = (Ne*Nout*pcon_out)/len(Sout.postsynaptic)
  else:
      See1_norm = 1.0
      See2_norm = 1.0
      Sie1_norm = 1.0
      Sie2_norm = 1.0
      Sii_norm  = 1.0
      Sei_norm  = 1.0
      Sin_norm  = 1.0
      Sout_norm = 1.0
      
  See1.we_ampa1 = See1_norm*Jee*(1-qee1)*bn.mV/tee_ampa
  See1.we_nmda1 = See1_norm*Jee*qee1*bn.mV/tee_nmda
  See1.we_ampa2 = 0.0*bn.mV/tee_ampa
  See1.we_nmda2 = 0.0*bn.mV/tee_nmda
  
  
  See2.we_ampa1 = 0.0*bn.mV/tee_ampa
  See2.we_nmda1 = 0.0*bn.mV/tee_nmda
  See2.we_ampa2 = See2_norm*Jee*(1-qee2)*bn.mV/tee_ampa
  See2.we_nmda2 = See2_norm*Jee*qee2*bn.mV/tee_nmda
  
  
  Sie1.we_ampa1 = 0.0*bn.mV/tie_ampa
  Sie1.we_nmda1 = 0.0*bn.mV/tie_nmda
  Sie1.we_ampa2 = Sie1_norm*Jie*(1-qie1)*bn.mV/tie_ampa
  Sie1.we_nmda2 = Sie1_norm*Jie*qie1*bn.mV/tie_nmda
  
  
  Sie2.we_ampa1 = Sie2_norm*Jie*(1-qie2)*bn.mV/tie_ampa
  Sie2.we_nmda1 = Sie2_norm*Jie*qie2*bn.mV/tie_nmda
  Sie2.we_ampa2 = 0.0*bn.mV/tie_ampa
  Sie2.we_nmda2 = 0.0*bn.mV/tie_nmda
  
  
  Sei.wi_gaba = Sei_norm*Jei*bn.mV/tei_gaba
  Sii.wi_gaba = Sii_norm*Jii*bn.mV/tii_gaba
  
  Sin.x = 1.0
  Sin.u = uff
  Sin.trec = trecff
  Sin.win_ampa = 0.0*Sin_norm*Jin*bn.mV/tee_ampa
  Sin.win_nmda = 0.0*Sin_norm*Jin*bn.mV/tee_nmda
  Sin.win_input = 1.0*Sin_norm*Jin*bn.mV/teo_input
  
  Sout.wout_ampa = Sout_norm*Jout*bn.mV/tee_ampa
  
  if normalize_connections == 'neuron':
      for i in iter(See1.synapses_post):          
          norm = (Ne*pcon/2)/len(i)
          See1.we_ampa1[i[:]] = norm*See1.we_ampa1[i[0]]
          See1.we_nmda1[i[:]] = norm*See1.we_nmda1[i[0]]
          
      for i in iter(See2.synapses_post):          
          norm = (Ne*pcon/2)/len(i)
          See2.we_ampa2[i[:]] = norm*See2.we_ampa2[i[0]]
          See2.we_nmda2[i[:]] = norm*See2.we_nmda2[i[0]]
          
      for i in iter(Sie1.synapses_post):          
          norm = (Ne*pcon/2)/len(i)
          Sie1.we_ampa2[i[:]] = norm*Sie1.we_ampa2[i[0]]
          Sie1.we_nmda2[i[:]] = norm*Sie1.we_nmda2[i[0]]
          
      for i in iter(Sie2.synapses_post):          
          norm = (Ne*pcon/2)/len(i)
          Sie2.we_ampa1[i[:]] = norm*Sie2.we_ampa1[i[0]]
          Sie2.we_nmda1[i[:]] = norm*Sie2.we_nmda1[i[0]]
          
      for i in iter(Sei.synapses_post):          
          norm = (Ni*pcon)/len(i)
          Sei.wi_gaba[i[:]] = norm*Sei.wi_gaba[i[0]]
          
      for i in iter(Sii.synapses_post):          
          norm = (Ni*pcon)/len(i)
          Sii.wi_gaba[i[:]] = norm*Sii.wi_gaba[i[0]]
      
    
  #==============================================================================
  #  Define monitors
  #==============================================================================
    
  Pe_ratemon = bn.PopulationRateMonitor(Pe,bin=1.0*bn.ms)
  Pi_ratemon = bn.PopulationRateMonitor(Pi,bin=1.0*bn.ms)
 
  Pe_spikes = bn.SpikeMonitor(Pe)
  
  #==============================================================================
  # Run model
  #==============================================================================
  timer = 0*bn.second
  t_start = time.time()
  bn.run(T, report='graphical')
  timer = timer + T
  print '-------------------------------------------------------'
  print 'Time is ' + str(timer)+' seconds'
  t_end = time.time()
  print 'Time to compute last ' +str(T)+' seconds is: ' + \
        str(t_end - t_start) + ' seconds'
  print '-------------------------------------------------------\n'
  
  #==============================================================================
  # Save into a Matlab file 
  #==============================================================================
  
  holder = {'Pe_rate':Pe_ratemon.rate,'Pe_time':Pe_ratemon.times,'exp_struct':exp_struct}
  scipy.io.savemat(save_file+'_network'+str(net_struct), mdict=holder)

  
  return {'Pe_rate':Pe_ratemon.rate,'Pe_times':Pe_ratemon.times,
          'Pi_rate':Pi_ratemon.rate, 'Pi_times':Pi_ratemon.times,
          'Pe_spikes': Pe_spikes.spiketimes,'T_step':bn.defaultclock.dt}

  bn.clear(erase=True, all=True)



#==============================================================================
#  Compute the frequency response within a window (to,tend)
#==============================================================================
def freq_comp(spikes,to,tend,fo,fend):

    T_step = spikes['Pe_times'][1]-spikes['Pe_times'][0]
    rates = spikes['Pe_rate'][(spikes['Pe_times']>=to) & (spikes['Pe_times']<=tend)]
    
    freqs = np.squeeze(np.fft.fftfreq(rates.size,T_step))
    inds_int = (freqs >= fo) & (freqs <= fend)
    
    rate_fft = np.abs(np.fft.fft(rates))[inds_int]
    rate_area = np.trapz(rate_fft,dx = freqs[1]-freqs[0],axis=0)
    freqs_list = freqs[inds_int]
        
    return {'freqs':freqs_list,'rate_fft':rate_fft,'rate_area':rate_area}
    
#==============================================================================
# Smooth binned output with a Gaussian filter   
#==============================================================================
def psth_smooth(rate, pnts, std):
    filt = gaussian(pnts, std)
    return filters.convolve1d(rate, filt/filt.sum())
    


#==============================================================================
# Define the parameters of the network: treats each parameter change as a case
#==============================================================================
def define_network_structure():

  stim1 = 0.0 * bn.second
  ISI1 = 0.0 * bn.second
  stim_reps1 = 0
  set1_begin = 0.0 * bn.second
  set1_end = set1_begin + (stim1+ISI1)*stim_reps1
  
  trec_base = 1000.0 * bn.ms
  u_base = 0.15
  
  delta_q = 0.0
  delta_tr= 0.0 * bn.ms
  delta_u = 0.0
  uff = 0.0
  trecff = 0.0*bn.ms
  Jin = 0.0
  
  baseline_mult = 0.8
  
  neuron_mult = 1.0
  
  k = 0
  w = 10.0
  
  ro = 15.0*bn.Hz
  
  return {'stim1':stim1,'ISI1':ISI1,'set1_begin':set1_begin,
    'set1_end':set1_end,
    'trec_base':trec_base,'u_base':u_base,
    'delta_q':delta_q,'delta_u':delta_u,'delta_tr':delta_tr,     
    'uff':uff,'trecff':trecff,'Jin':Jin,'neuron_mult':neuron_mult,
    'baseline_mult':baseline_mult,'k':k,'w':w,'ro':ro}


if __name__ == '__main__':

    net_struct = 3
    num_reps = 1
    new_connectivity =False
    #
    exp_struct = define_network_structure()    
      
    stim1 = 2.0
    ISI1 = 2.0
    burnin = 0.5
    exp_struct['stim1'] = stim1 * bn.second
    exp_struct['ISI1'] = ISI1 * bn.second
    exp_struct['stim_reps1'] = 1
    exp_struct['set1_begin'] = burnin * bn.second
    #  exp_struct['set1_end'] = 1.7*bn.second
    exp_struct['set1_end'] = (exp_struct['set1_begin'] + 
                (exp_struct['stim1']+exp_struct['ISI1'])*exp_struct['stim_reps1'])

    exp_struct['trec_base'] = 500.0 * bn.ms

    exp_struct['delta_q'] = 0.0
    exp_struct['delta_u'] = 0.0
    exp_struct['delta_tr'] = 0.00 * bn.ms
    
    exp_struct['trecff'] = 500.0*bn.ms
    exp_struct['uff'] = 0.0
    
    exp_struct['ro'] = 20.0*bn.Hz

    exp_struct['k'] = 0.0
    exp_struct['neuron_mult'] = net_struct
    exp_struct['baseline_mult'] = 0.365
    
#    exp_struct['w'] = 0.146
#    u_hold = [0.0,0.05,0.1,0.2]
#    Jin_hold = [0.011318,0.049,0.069,0.094]
#    file_base = 'pfdbk_final4_u'
    
    exp_struct['w'] = 0.145
    u_hold = [0.0,0.05,0.1,0.2]
    Jin_hold = [0.014,0.049,0.069,0.094]
    file_base = 'pfdbk_final4_u'
    
#    exp_struct['w'] = 0.146
#    u_hold = [0.0]
#    Jin_hold = [0.011318]
#    file_base = 'pfdbk_testing_w146Jin011318_u'
    
    for i in range(len(u_hold)):
        local_u = u_hold[i]
        local_Jin = Jin_hold[i]
        exp_struct['u_base'] = local_u 
        exp_struct['Jin'] = local_Jin
        save_file = file_base+str(int(100*local_u))
        print 'u='+str(local_u)+'   Jin='+str(local_Jin)+'\n Save File: '+save_file
        hold = rep_supp(exp_struct, net_struct,new_connectivity,save_file)   
        
        Pe_rate = hold['Pe_rate']
        T = hold['Pe_times']
        plt.plot(T,psth_smooth(Pe_rate,79,5))      
        plt.savefig(save_file+'.png')
     

  
  
  
  

  

  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  