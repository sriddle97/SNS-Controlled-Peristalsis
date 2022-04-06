function P = Riddle_cycle_new(x0,Isens0,T)
%% Simulation turned into a function for stability analysis %%
%  Shane Riddle, Zhuojun Yu                                  %
%  Last edited 04/05/2022                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input x0: initial condition; Isens0: initial censor currents; T: period
% Output P: final point

%% Neuron and Synapse Properties
% Units are nF, uS, mV, ms, nA
C = 5;
Gm = 1;
tauHmax = 300;      %300 is default (200 and 100 also used for faster action)
Er = -60;           %resting potential, biologically typical value
delta = .01;        %bifurcation thingy
R = 20;
offset = 0.1;       %R offset for CPG stopped states
Nseg = 6;               %number of segments

% Speed gain for controller (if using instead of max/min Speed)
k=1.3;      %1.3 is min value when tauHmax is 300 (1.9 for 200, 3.6 for 100)
% Period for tauHmax=300 and k=1.3 is 5250ms (measured by sensor signal)
% Already "stable" after the first cycle

%makes sense biologically that excitatory delE higher than inhibitory
%Inhibitory Synapse Porperties (???? transmitter)
Einh = -100;
delEsyn = Einh - Er;

%Excitatory Synapse Properties (calcium transmitter)
Eex = 134;
delEsyn_ex = Eex-Er;

%Persistent sodium (NaP) channel
Ena = 50;

S = .05; %Slope of the sigmoid of hInf, mInf.
delEna = Ena - Er; %50 mV
Utest = -R:.1:3*R;
mInf = @(U) 1./(1 + exp(S*(R-U))); %Steady-state value of m, as a function of U.
hInf = @(U) 1./(1 + .5*exp(S*U)); %Steady-state value of h, as a function of U.
tauh = @(U) tauHmax*hInf(U).*sqrt(.5*exp(S*U)); %Time constant of h, as a function of U.

%Solve for the conductance of the NaP channel to get U* = R, with no
%external current or synaptic inputs.
Gna = Gm*R/(mInf(R)*hInf(R)*(delEna - R));

%Now we know that U* = R, and we can find h* and m* based on that.
Ustar = R;
mStar = mInf(Ustar);
hStar = hInf(Ustar);


%% Simulation Time(ms) and Stimulation Stuff
dtSim = 1;
tmax = T;
tSim = 0:dtSim:tmax;
numSteps = length(tSim);


%% SNS Simulation Setup
Isens = zeros(length(tSim),Nseg);
Isens(1,:) = Isens0;

%Initialize U(t) for each neuron, can control starting voltage for each
Usim1 = zeros(length(tSim),Nseg);
Usim2 = zeros(length(tSim),Nseg);
Usim3 = zeros(length(tSim),Nseg);
Usim4 = zeros(length(tSim),Nseg);
Usim1(1,:) = x0(1:6); Usim2(1,:) = x0(7:12); Usim3(1,:) = x0(13:18); Usim4(1,:) = x0(19:24);

hSim1 = zeros(length(tSim),Nseg);
hSim2 = zeros(length(tSim),Nseg);
hSim1(1,:) = x0(25:30); hSim2(1,:) = x0(31:36);

% Set up synapses and function to use later
%Calculate gSyn for the CPG neurons
gSyn = (-delta - delta*Gna*mInf(delta)*hInf(delta) + Gna*mInf(delta)*hInf(delta)*delEna)/(delta - delEsyn);
g = @(U,gsyn) min(max(U/R,0),1)*gsyn;

% CPG synapses
Esyn1 = Einh;
Esyn2 = Einh;
delEsyn1 = Esyn1-Er;
delEsyn2 = Esyn2-Er;
gSyn1 = gSyn;
gSyn2 = gSyn;

% Excitation synapses
% 3-2
Esyn32 = Eex;
delEsyn32 = Esyn32-Er;
gSyn32 = ((delEsyn1*gSyn1*(-offset))/((offset+R)*R)-1-gSyn1*(-offset)/R)/(1-delEsyn32/(R+offset));
% 2-4
Esyn24 = Eex;
delEsyn24 = Esyn24-Er;
gSyn24 = R/((delEsyn24/R-1)*(R+offset));

% Inhibition synapses
% 3-1
Esyn31 = Einh;
delEsyn31 = Esyn31-Er;
gSyn31 = ((delEsyn2*gSyn2*(R+offset))/(-offset*R)-1-gSyn2*(R+offset)/R)/(1+delEsyn31/(-offset));
% 3-4
Esyn34 = Einh;
delEsyn34 = Esyn34-Er;
gSyn34 = (-delEsyn24*gSyn24*(R+offset))/(delEsyn34*R);
% 4-3
Esyn43 = Einh;
delEsyn43 = Esyn43-Er;
gSyn43 = -R/delEsyn43;

% Initialize Synapses
Gsyn1 = zeros(length(tSim),Nseg);
Gsyn2 = zeros(length(tSim),Nseg);
Gsyn32 = zeros(length(tSim),Nseg);
Gsyn24 = zeros(length(tSim),Nseg);
Gsyn31 = zeros(length(tSim),Nseg);
Gsyn34 = zeros(length(tSim),Nseg);
Gsyn43 = zeros(length(tSim),Nseg);

for q=1:Nseg
    Gsyn1(1,q) = g(Usim1(1,q),gSyn1);
    Gsyn2(1,q) = g(Usim2(1,q),gSyn2);
    Gsyn32(1,q) = g(Usim3(1,q),gSyn32);
    Gsyn24(1,q) = g(Usim2(1,q),gSyn24);
    Gsyn31(1,q) = g(Usim3(1,q),gSyn31);
    Gsyn34(1,q) = g(Usim3(1,q),gSyn34);
    Gsyn43(1,q) = g(Usim4(1,q),gSyn43);
end


%% Kinematics Simulation Setup

hmax = 11;              %height max, cm
hmin = 6.5;             %height min, cm
l0 = 7.3;               %rhombus side length, cm
lmax = sqrt(-1*(hmin^2)+((2*l0)^2));        %cm
lmin = sqrt(-1*(hmax^2)+((2*l0)^2));        %cm
% Speed_nl = 6.39 cm/s  at 12 V (no-load motor speed, linear actuation)   or 2.09?
% Speed = 2*pi*20/60;     %cm/s  From performance graph, using 0.6 N-m of Torque?
maxSpeed = 2*pi*20/60/300;         % technically cm/ms but arbitrary for now

%initialize all segments in expanded state
h = zeros(length(tSim),Nseg);
h(1,:) =  x0(37:42);
l = zeros(length(tSim),Nseg);
l(1,:) = sqrt(-1*(h(1,:).^2)+((2*l0)^2));
dh = zeros(length(tSim),Nseg);


%% Simulation Execution
%can turn neurons on/off to debug
of1 = 1;
of2 = 1;
of3 = 1;
of31 = 1;
of32 = 1;
of34 = 1;
of4 = 1;

for i=2:numSteps
    for j=1:Nseg
        % Neural Network Forward Euler Calculations
        % U1
        % inhib from other CPG neuron (U2), inhib from U3
        Usim1(i,j) = Usim1(i-1,j) + dtSim/C*(-Gm*Usim1(i-1,j) + of2*Gsyn2(i-1,j)*(delEsyn-Usim1(i-1,j)) + of3*of31*Gsyn31(i-1,j)*(delEsyn31-Usim1(i-1,j)) + Gna*mInf(Usim1(i-1,j))*hSim1(i-1,j)*(delEna-Usim1(i-1,j)));
        hSim1(i,j) = hSim1(i-1,j) + dtSim/tauh(Usim1(i-1,j))*(hInf(Usim1(i-1,j)) - hSim1(i-1,j));
        Gsyn1(i,j) = g(Usim1(i,j),gSyn1);

        % U2
        % inhib from other CPG neuron (U1), excite from U3
        Usim2(i,j) = Usim2(i-1,j) + dtSim/C*(-Gm*Usim2(i-1,j) + of1*Gsyn1(i-1,j)*(delEsyn-Usim2(i-1,j)) + of3*of32*Gsyn32(i-1,j)*(delEsyn32-Usim2(i-1,j)) + Gna*mInf(Usim2(i-1,j))*hSim2(i-1,j)*(delEna-Usim2(i-1,j)));
        hSim2(i,j) = hSim2(i-1,j) + dtSim/tauh(Usim2(i-1,j))*(hInf(Usim2(i-1,j)) - hSim2(i-1,j));
        Gsyn2(i,j) = g(Usim2(i,j),gSyn2);
        Gsyn24(i,j) = g(Usim2(i,j),gSyn24);

        % U3 interneuron
        % stim current from sensor output (if statement), inhib from previous segment's U4 (j-1)
        % if statement to loop Nth segment 4-3 synapse back to 1st segment
        if j == 1
            Usim3(i,j) = Usim3(i-1,j) + dtSim/C*(-Gm*Usim3(i-1,j) + of4*Gsyn43(i-1,Nseg)*(delEsyn43-Usim3(i-1,j)) + Isens(i-1,j));
        else
            Usim3(i,j) = Usim3(i-1,j) + dtSim/C*(-Gm*Usim3(i-1,j) + of4*Gsyn43(i-1,j-1)*(delEsyn43-Usim3(i-1,j)) + Isens(i-1,j));
        end
        Gsyn31(i,j) = g(Usim3(i,j),gSyn31);
        Gsyn32(i,j) = g(Usim3(i,j),gSyn32);
        Gsyn34(i,j) = g(Usim3(i,j),gSyn34);

        % U4 interneuron
        % excite from expand neuron (U2), inhib from U3
        Usim4(i,j) = Usim4(i-1,j) + dtSim/C*(-Gm*Usim4(i-1,j) + of2*Gsyn24(i-1,j)*(delEsyn24-Usim4(i-1,j)) + of3*of34*Gsyn34(i-1,j)*(delEsyn34-Usim4(i-1,j)));
        Gsyn43(i,j) = g(Usim4(i,j),gSyn43);


        %%%%%%%% Robot Kinematics Calculations (Position Control) %%%%%%%%
        % Usim1 for diameter contraction and Usim2 for diameter expansion
        cont_sig = -Usim1(i,j)+Usim2(i,j);   %expand when +, contract when -

        % Cut off target positions beyond robot physical limitations 
        cont_sig_corrected = min(max(cont_sig,-R),R);

        % mapping from (-R,R) signal range to (hmin,hmax) position range
        pos_target = cont_sig_corrected*(hmax-hmin)/(2*R)+(hmax+hmin)/2;      
        
        % Position target sets motor speed (all or nothing)
        % This should also shut down motor when physical limits reached
%         if pos_target<=hmin && h(i-1,j)>hmin
%             Speed = -maxSpeed;
%         elseif pos_target>=hmax && h(i-1,j)<hmax
%             Speed = maxSpeed;
%         else
%             Speed = 0;
%         end

        % Proportional control for Speed
        % k is the gain, defined in properties section up top
        k=1.3;
        Speed = k*(pos_target-h(i-1,j))*maxSpeed;
        if h(i-1,j) <= hmin
            Speed = max(0,Speed);
        elseif h(i-1,j) >= hmax
            Speed = min(0,Speed);
        end
        
        
        % Segment geometry calculations
        dh(i,j) = Speed*dtSim;
        h(i,j) = h(i-1,j)+dh(i,j);
        l(i,j) = sqrt(-1*(h(i,j)^2)+((2*l0)^2));
        %vertices calculated below (left and bottom constrained)
        vert = [0 h(i,j)/2; l(i,j)/2 0; l(i,j) h(i,j)/2; l(i,j)/2 h(i,j); 0 h(i,j)/2];

        % If statement approximation for stretch sensor current (changed sensor threshold!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
        if l(i,j) <= lmin+0.01
            Isens(i,j) = R;
        else
            Isens(i,j) = 0;
        end
    end
    
end


P = [Usim1(end,:) Usim2(end,:) Usim3(end,:) Usim4(end,:) hSim1(end,:) hSim2(end,:) h(end,:)];

end