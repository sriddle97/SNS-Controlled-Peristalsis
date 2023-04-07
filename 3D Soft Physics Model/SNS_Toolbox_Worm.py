## SNS controller for peristaltic locomotion of a simulated worm robot
# Author: Shane Riddle
# Last edited: 03/12/2023
#########################################################################

# The motors use an on-board microcontroller so you can give them a speed
# command and they will run at that speed. Can also give a position target
# and they will run at max speed to reach said target. Using poistion
# control in this simulation for better accuracy

########################### Import packages #############################

# SNS-Toolbox Packages
from sns_toolbox.neurons import NonSpikingNeuron, NonSpikingNeuronWithPersistentSodiumChannel
from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.networks import Network
from sns_toolbox.renderer import render

# Basic Utility Packages
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

# Movie Making Packages  (USE mediapy FOR Mujoco stuff) or moviepy?????
import matplotlib.animation as animation


########################### Define functions ############################

## m and h steady state and time constant functions (z is a placeholder)
# Steady-state value of z, as a function of U.
def zInf(U, Kz, Sz, Ez):
    return 1/(1+Kz*np.exp(Sz*(Ez-U)))

# Time constant of z, as a function of U.
def tauz(U, tauzmax, Kz, Sz, Ez):
    return tauzmax*zInf(U, Kz, Sz, Ez)*np.sqrt(Kz*np.exp(Sz*(Ez-U)))

## Synaptic conductance function
def g(U, gsyn, R):
    return np.minimum(np.maximum(U/R, 0), 1)*gsyn



######################### Network Parameters ############################

## Neuron and synapse properties (Units are nF, uS, mV, ms, nA)
C = 5
Gm = 1
tauHmax = 300       # 300 is default (200 and 100 also used for faster action)
Er = -60            # resting potential, biologically typical value
delta = .01         # bifurcation thingy
R = 20
offset = 0.1        # R offset for CPG stopped states
Nseg = 3            # 6, number of segments

# Speed gain for controller (if using instead of max/min Speed)
k = 1.3             # 1.3 is min value when tauHmax is 300 (1.9 for 200, 3.6 for 100)

# Stability Notes: Period for tauHmax=300 and k=1.3 is 5250ms (measured by sensor signal). Already "stable" after the first cycle.

# makes sense biologically that excitatory delE higher than inhibitory
# Inhibitory Synapse Porperties (???? transmitter)
Einh = -100
delEsyn = Einh - Er

# Excitatory Synapse Properties (calcium transmitter)
Eex = 134
delEsyn_ex = Eex-Er

# Persistent sodium (NaP) channel
Ena = 50
S = .05                 # Slope of the sigmoid of hInf, mInf.
Kh = 0.5
Eh = 0
Km = 1
Em = R
delEna = Ena - Er       # 50 mV ??

# Solve for the conductance of the NaP channel to get U* = R, with no external current or synaptic inputs.
Gna = Gm*R/(zInf(R, Km, S, Em)*zInf(R, Kh, -S, Eh)*(delEna - R))

# Now we know that U* = R, and we can find h* and m* based on that.
Ustar = R
mStar = zInf(Ustar, Km, S, Em)
hStar = zInf(Ustar, Kh, -S, Eh)



############################## SNS Setup #################################

## Set up synapse parameters
# CPG synapses
Esyn1_2 = Einh
delEsyn1_2 = Esyn1_2-Er
gSyn_cpg = (-delta - delta*Gna*zInf(delta, Km, S, Em)*zInf(delta, Kh, -S, Eh) + Gna*zInf(delta, Km, S, Em)*zInf(delta, Kh, -S, Eh)*delEna)/(delta - delEsyn)
gSyn1_2 = gSyn_cpg

# Excitation synapses
# 3-2
Esyn32 = Eex
delEsyn32 = Esyn32-Er
gSyn32 = ((delEsyn1_2*gSyn1_2*(-offset))/((offset+R)*R)-1-gSyn1_2*(-offset)/R)/(1-delEsyn32/(R+offset))
# 2-4
Esyn24 = Eex
delEsyn24 = Esyn24-Er
gSyn24 = R/((delEsyn24/R-1)*(R+offset))

# Inhibition synapses
# 3-1
Esyn31 = Einh
delEsyn31 = Esyn31-Er
gSyn31 = ((delEsyn1_2*gSyn1_2*(R+offset))/(-offset*R)-1-gSyn1_2*(R+offset)/R)/(1+delEsyn31/(-offset))
# 3-4
Esyn34 = Einh
delEsyn34 = Esyn34-Er
gSyn34 = (-delEsyn24*gSyn24*(R+offset))/(delEsyn34*R)
# 4-3
Esyn43 = Einh
delEsyn43 = Esyn43-Er
gSyn43 = -R/delEsyn43



# Define neuron and synapse types
neuron = NonSpikingNeuron(membrane_capacitance = C, membrane_conductance = Gm)
cpg_neuron = NonSpikingNeuronWithPersistentSodiumChannel(membrane_capacitance = C, membrane_conductance = Gm,
                                                         g_ion = [Gna], e_ion = [delEna],
                                                         k_m = [Km], slope_m = [S], e_m = [Em],
                                                         k_h = [Kh], slope_h = [-S], e_h = [Eh], tau_max_h = [tauHmax])

synapse_32 = NonSpikingSynapse(max_conductance = gSyn32, reversal_potential = delEsyn32)
synapse_24 = NonSpikingSynapse(max_conductance = gSyn24, reversal_potential = delEsyn24)
synapse_31 = NonSpikingSynapse(max_conductance = gSyn31, reversal_potential = delEsyn31)
synapse_34 = NonSpikingSynapse(max_conductance = gSyn34, reversal_potential = delEsyn34)
synapse_43 = NonSpikingSynapse(max_conductance = gSyn43, reversal_potential = delEsyn43)
cpg_synapse = NonSpikingSynapse(max_conductance = gSyn1_2, reversal_potential = delEsyn1_2)


## Build SNS network
net = Network()
U1 = 'U1_'
U2 = 'U2_'
U3 = 'U3_'
U4 = 'U4_'
# if statements initialize neuron voltages for segment 1 contracted and the rest expanded
for i in range(Nseg):
    # CPG
    if i == 0:
        net.add_neuron(cpg_neuron, name = U1+str(i+1), color = 'green', initial_value = R)
    else:
        net.add_neuron(cpg_neuron, name = U1+str(i+1), color = 'green', initial_value = -offset)
    net.add_input(U1+str(i+1)) 
    net.add_output(U1+str(i+1))
    if i == 0:
        net.add_neuron(cpg_neuron, name = U2+str(i+1), color = 'green', initial_value = 0)
    else:
        net.add_neuron(cpg_neuron, name = U2+str(i+1), color = 'green', initial_value = R+offset)
    net.add_input(U2+str(i+1)) 
    net.add_output(U2+str(i+1))
    net.add_connection(cpg_synapse, U1+str(i+1), U2+str(i+1))
    net.add_connection(cpg_synapse, U2+str(i+1), U1+str(i+1))
    # Interneurons
    if i == 0:
        net.add_neuron(neuron, name = U3+str(i+1), color = 'blue', initial_value = 0)
    else:
        net.add_neuron(neuron, name = U3+str(i+1), color = 'blue', initial_value = R)
    net.add_input(U3+str(i+1), color = 'red')
    net.add_neuron(neuron, name = U4+str(i+1), color = 'blue', initial_value = 0)
    net.add_connection(synapse_32, U3+str(i+1), U2+str(i+1))
    net.add_connection(synapse_24, U2+str(i+1), U4+str(i+1))
    net.add_connection(synapse_31, U3+str(i+1), U1+str(i+1))
    net.add_connection(synapse_34, U3+str(i+1), U4+str(i+1))
    # Synapse between segments
    if i != 0:
        net.add_connection(synapse_43, U4+str(i), U3+str(i+1))
net.add_connection(synapse_43, U4+str(Nseg), U3+str(1))      # from last to first segment

# Simulation Time(ms) and Stimulation Stuff
dtSim = 1
tmax = 5000            # 10000
tSim = np.arange(0, tmax, dtSim)
numSteps = np.size(tSim)

# Perturbation current (to get the ball rolling on segment 1)
Ipert = np.zeros((np.size(tSim), Nseg))
Ipert[0,0] = 1

Iapp = np.zeros((np.size(tSim), Nseg))     # none for now

# Initialize sensor current signal vectors
Isens = np.zeros((np.size(tSim), Nseg))
Isens[0,:] = R      # this was designed
Isens[0,0] = 0      # Since seg 1 starts contracted

# Initialize U(t) for each neuron, can control starting voltage for each
Usim1 = np.zeros((np.size(tSim), Nseg))
Usim2 = np.zeros((np.size(tSim), Nseg))
# Usim3 = np.zeros((np.size(tSim), Nseg))
# Usim4 = np.zeros((np.size(tSim), Nseg))
Usim1[0,:] = -offset
Usim1[0,0] = R
Usim2[0,:] = R+offset
Usim2[0,0] = 0



###################### Kinematics Simulation Setup ############################

hmax = 11                                                       # height max, cm
hmin = 6.5                                                      # height min, cm
l0 = 7.3                                                        # rhombus side length, cm
lmax = math.sqrt(-1*(math.pow(hmin,2))+(math.pow(2*l0,2)))      # cm
lmin = math.sqrt(-1*(math.pow(hmax,2))+(math.pow(2*l0,2)))      # cm
# Speed_nl = 6.39 cm/s  at 12 V (no-load motor speed, linear actuation)   or 2.09?
# Speed = 2*pi*20/60;     # cm/s  From performance graph, using 0.6 N-m of Torque?
maxSpeed = 2*math.pi*20/60/300          # technically cm/ms but arbitrary for now

# initialize all segments in expanded state except first (need first two time steps for sim)
dh = np.zeros((np.size(tSim), Nseg))
h = np.zeros((np.size(tSim), Nseg))
l = np.zeros((np.size(tSim), Nseg))
for i in range(Nseg):
    # if i == 0:
    #     h[0,i] = hmin
    # else:
    #     h[0,i] = hmax


    # FIX THE REST OF THE INITIALIZATIONS FOR THE NETWORK????
    h[0,i] = hmax
    l[0,i] = math.sqrt(-1*(math.pow(h[0,i],2))+(math.pow(2*l0,2)))



##################### SNS Network and Mujoco Simulation ######################

model = net.compile(backend='numpy', dt=dtSim)      # compiles the network
data = np.zeros(net.get_num_outputs_actual())       # get_num_outputs_actual gives size of output vector

inputs_mat = np.zeros((np.size(tSim), Nseg, 3))     # 3 inputs (at U1, U2, and U3)
inputs_mat[:,:,0] = Ipert + Iapp                    # Goes into U1
inputs_mat[:,:,1] = -Iapp                            # Goes into U2
inputs_mat[:,:,2] = Isens                           # Goes into U3

# input vector generated in order inputs were defined when building network
# transform matrix to 1D vector in order
inputs = np.zeros(Nseg*3)
for j in range(Nseg):
    inputs[(3*j):3*(j+1)] = inputs_mat[0,j,:]

# initialize vertices matrix for 2d motion plotting
vert = np.zeros((5, 2, Nseg, np.size(tSim)))        # (5 points per rhombus, 2 dimensions (x and y), Nseg segments, timesteps)

# run network simulation
for i in range(numSteps):
    # data vector generated in order outputs were defined when building network
    # print(inputs)
    data = model(inputs)

    #print(data)

    # send data to Usim variables and do kinematics/feedback to define inputs for next time step
    for j in range(Nseg):

        if i != numSteps-1:
            Usim1[i+1,j] = data[(2*j)]
            Usim2[i+1,j] = data[(2*j)+1]

        
            ######## Robot Kinematics Calculations (Position Control) ########
            # Usim1 for diameter contraction and Usim2 for diameter expansion
            cont_sig = -Usim1[i+1,j]+Usim2[i+1,j]      # expand when +, contract when -
            
            #print(cont_sig)

            # Cut off target positions beyond robot physical limitations 
            cont_sig_corrected = np.minimum(np.maximum(cont_sig, -R), R)

            # mapping from (-R,R) signal range to (hmin,hmax) position range
            pos_target = cont_sig_corrected*(hmax-hmin)/(2*R)+(hmax+hmin)/2

            # if j == 0:
            #     print(pos_target)

            # Proportional control for Speed
            # k is the gain, defined in properties section up top
            Speed = k*(pos_target-h[i,j])*maxSpeed
            if h[i,j] <= hmin:
                Speed = np.maximum(0, Speed)
            elif h[i,j] >= hmax:
                Speed = np.minimum(0, Speed)

            # Segment geometry calculations
            dh[i+1,j] = Speed*dtSim
            h[i+1,j] = h[i,j]+dh[i,j]
            l[i+1,j] = math.sqrt(-1*(math.pow(h[i,j], 2))+(math.pow(2*l0, 2)))
            # vertices calculated below (left and bottom constrained)
            #vert = np.array([0, h[i+1,j]/2],
            #                [l[i+1,j]/2, 0],
            #                [l[i+1,j], h[i+1,j]/2],
            #                [l[i+1,j]/2, h[i+1,j]],
            #                [0, h[i+1,j]/2])

            # If statement approximation for stretch sensor current (may need to change)
            if l[i+1,j] <= lmin+0.01:
                Isens[i+1,j] = R
            else:
                Isens[i+1,j] = 0

            #inputs for next time step (Ipert and Iapp don't change in this network so only Isens needs redefined)
            inputs_mat[i+1,j,2] = Isens[i+1,j]
            inputs[(3*j):3*(j+1)] = inputs_mat[i+1,j,:]



############################## Plots ###############################

# Plot 1: CPG U1 neuron potentials for each time step
# Plot 2: CPG U2 neuron potentials for each time step
# Plots 3-Nseg: Change in segemnt height per time step (basically speed)
# color_mat = np.array([[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]])

color_mat = cm.rainbow(np.linspace(0, 1, Nseg))

fig = plt.figure()
Segment = "Segment"
for ct in range(1, Nseg+3):
    ax = fig.add_subplot(Nseg+2, 1, ct)
    if ct == 1:
        for j in range(Nseg):
            ax.plot(tSim, Usim1[:,j], color = color_mat[j,:])
        ax.set_ylabel('Voltage (mV)')
        ax.set_title('CPG U1 (Contraction)')
    elif ct == 2:
        for j in range(Nseg):
            ax.plot(tSim, Usim2[:,j]*1/100, color = color_mat[j,:])
        ax.set_ylabel('Voltage (mV)')
        ax.set_title('CPG U2 (Expansion)')
    else:
        # ax.plot(tSim, h[:,ct-3], color = color_mat[ct-3,:])
        # ax.set_ylabel('w (m)')
        # ax.set_title(Segment + str(ct-2))

        ax.plot(tSim, l[:,ct-3]*1/100, color = color_mat[ct-3,:])
        ax.set_ylabel('l (m)')
        ax.set_title(Segment + str(ct-2))
plt.xlabel('Time (ms)')
plt.show()