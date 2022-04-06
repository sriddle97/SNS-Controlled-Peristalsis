%%  Synthetic Nervous System Controller for Worm Robot Simulation  %%
%   Shane Riddle                                                    %
%   Last edited 04/05/2022                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Kinametaic Model (Speed proportional to voltage)
clear
clc
% The motors use an on-board microcontroller so you can give them a speed
% command and they will run at that speed. Can also give a position target
% and they will run at max speed to reach said target. Using poistion
% control in this simulation for better accuracy.


%% Neuron and Synapse Properties
% Units are nF, uS, mV, ms, nA
C = 5;
Gm = 1;
tauHmax = 300;      %300 is default (200 and 100 also used for faster action)
Er = -60;           %resting potential, biologically typical value
delta = .01;        %bifurcation thingy
R = 20;
offset = 0.1;       %R offset for CPG stopped states
Nseg = 3;               %6, number of segments

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
tmax = 5000;           %10000
tSim = 0:dtSim:tmax;
numSteps = length(tSim);

%Perturbation current (to get the ball rolling on segment 1)
Ipert = zeros(length(tSim),Nseg);
Ipert(1,1) = 1;

Iapp = zeros(length(tSim),Nseg);    %none for now
%%%%%%%%%%%%%% Current Pulse (FOR TESTING AND DEBUGGING) %%%%%%%%%%%%%%%%%
% I = 0.1;              % magnitude
% tStart = 2000;
% tEnd = 4000;
% Iapp = I*(tSim >= tStart & tSim <= tEnd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Stability Check (from provided code)
%Once we know all of the equilibrium states, we can compute the Jacobian
%matrix and its eigenvalues. This will tell us about the stability of each
%neuron in the network (although not about the system as a whole).
%Specifically, we want to avoid complex eigenvalues, because this will
%cause unwanted oscillations that will make our analysis fail.
dU = 1e-3;
dm_dU = @(U) (mInf(U+dU) - mInf(U))/dU;
dh_dU = @(U) (hInf(U+dU) - hInf(U))/dU;
dtauh_dU = @(U) (tauh(U+dU) - tauh(U))/dU;

dUdot_dU = 1/C*(-1 + dm_dU(Ustar)*Gna*hStar*(delEna - Ustar) - Gna*mStar*hStar);
dUdot_dh = 1/C*(Gna*mStar*(delEna - Ustar));
dhdot_dU = dh_dU(Ustar)/tauh(Ustar);
dhdot_dh = -1/tauh(Ustar);

J = [   dUdot_dU   dUdot_dh;...
        dhdot_dU   dhdot_dh];

%Print warning if the Jacobian has complex eigenvalues.
eigs(J)
if any(imag(eigs(J)))
    warning('on')
    warning('The active neuron has complex eigenvalues, which may complicate analysis.')
end


%% SNS Simulation Setup
Isens = zeros(length(tSim),Nseg);
Isens(1,:) = R;     % this was designed
Isens(1,1) = 0;     % Since seg 1 starts contracted

%Initialize U(t) for each neuron, can control starting voltage for each
Usim1 = zeros(length(tSim),Nseg);
Usim2 = zeros(length(tSim),Nseg);
Usim3 = zeros(length(tSim),Nseg);
Usim4 = zeros(length(tSim),Nseg);

Usim1(1,:) = -offset;
Usim1(1:1) = R;     % Since seg 1 starts contracted
Usim2(1,:) = R+offset;
Usim2(1:1) = 0;     % Since seg 1 starts contracted
Usim3(1,:) = R;
Usim3(1,1) = 0;     % Since seg 1 starts contracted
Usim4(1,:) = 0;

hSim1 = zeros(length(tSim),Nseg);
hSim1(1,:) = hInf(Usim1(1));
hSim2 = zeros(length(tSim),Nseg);
hSim2(1,:) = hInf(Usim2(1));

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
h(1,:) = hmax;
h(1,1) = hmin;
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

% Video capture stuff
myVideo = VideoWriter('Single_Wave.mp4'); %open video file
myVideo.FrameRate = 100;  %can adjust this
open(myVideo)

for i=2:numSteps
    for j=1:Nseg
        % Neural Network Forward Euler Calculations
        % U1
        % inhib from other CPG neuron (U2), inhib from U3
        Usim1(i,j) = Usim1(i-1,j) + dtSim/C*(-Gm*Usim1(i-1,j) + of2*Gsyn2(i-1,j)*(delEsyn-Usim1(i-1,j)) + of3*of31*Gsyn31(i-1,j)*(delEsyn31-Usim1(i-1,j)) + Gna*mInf(Usim1(i-1,j))*hSim1(i-1,j)*(delEna-Usim1(i-1,j)) + Iapp(i-1,j) + Ipert(i-1,j));
        hSim1(i,j) = hSim1(i-1,j) + dtSim/tauh(Usim1(i-1,j))*(hInf(Usim1(i-1,j)) - hSim1(i-1,j));
        Gsyn1(i,j) = g(Usim1(i,j),gSyn1);

        % U2
        % inhib from other CPG neuron (U1), excite from U3
        Usim2(i,j) = Usim2(i-1,j) + dtSim/C*(-Gm*Usim2(i-1,j) + of1*Gsyn1(i-1,j)*(delEsyn-Usim2(i-1,j)) + of3*of32*Gsyn32(i-1,j)*(delEsyn32-Usim2(i-1,j)) + Gna*mInf(Usim2(i-1,j))*hSim2(i-1,j)*(delEna-Usim2(i-1,j)) - Iapp(i-1,j));
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
    
%     % Uncomment this section and the two lines after the for loop end to
%     % make and play the simulation video. Warning, it's pretty slow.
%     for j = 1:Nseg
%         le = [sum(l(i,j+1:Nseg)) hmax/2];
%         bot = [sum(l(i,j+1:Nseg))+l(i,j)/2 hmax/2-h(i,j)/2];
%         ri = [sum(l(i,j:Nseg)) hmax/2];
%         top = [sum(l(i,j+1:Nseg))+l(i,j)/2 hmax/2+h(i,j)/2];
%         vert(:,:,j) = [le; bot; ri; top; le];
%         plot(vert(:,1,j),vert(:,2,j),'k')
%         hold on
%     end
%     daspect([1 1 1])
%     % You may need to change the next line to fit the plot on your screen
%     set(gcf,'Position',[1000 500 2000 500])
%     xlim([0 lmax*Nseg])
%     ylim([0 12])
%     drawnow
%     F(i) = getframe;
%     hold off
%     writeVideo(myVideo, F(i));
end
% close(myVideo)
% movie(F,1,50)

% Plot each segment height at each time step individually
figure
for j=1:Nseg
    subplot(Nseg,1,j)
    plot(tSim,h(:,j))
%     xlabel('time (ms)')
    ylabel('h (cm)')
    title(['Segment ' num2str(j)])
end
subplot(Nseg,1,Nseg)
xlabel('time (ms)')
% title('Segment Heights')

figure
for j=1:Nseg
    plot(tSim,h(:,j))
    hold on
    xlabel('time (ms)')
    ylabel('h (cm)')
    title(['Segment Heights'])
end
hold off

% Plot change in segemnt height per time step (basically speed)
figure
subplot(6,1,5)
for j=1:Nseg
    plot(tSim,dh(:,j))
    hold on
end
xlabel('time (ms)')
ylabel('dh (cm)')
title('dh (speed) plot')
% legend(h_dh_names)
hold off

% Plot segment height at each time step
subplot(6,1,6)
for j=1:Nseg
    plot(tSim,h(:,j))
    hold on
end
xlabel('time (ms)')
ylabel('h (cm)')
title('segment height plot')
% legend(h_dh_names)
hold off

% Plot CPG neuron potentials for each time step
subplot(6,1,1)
scale_factor = 100; %pick one so you can see dh overlayed on CPG signals
for j=1:Nseg
    % plot(tSim,Usim1(:,j),tSim,Usim2(:,j),tSim,dh(:,j)*scale_factor,'--k')
    plot(tSim,Usim1(:,j),tSim,Usim2(:,j))
%     plot(tSim,Usim1(:,j))
    hold on
end
xlabel('time (ms)')
ylabel('voltage (mV)')
title('CPG voltage plot')
% legend(Usim1_Usim2_names)
hold off

% Plot other neuron potentials for each time step
subplot(6,1,2)
for j=1:Nseg
    plot(tSim,Usim3(:,j))
    hold on
end
xlabel('time (ms)')
ylabel('voltage (mV)')
title('Neuron 3 voltage plot')
% legend(Usim1_Usim2_names)
hold off

subplot(6,1,3)
for j=1:Nseg
    plot(tSim,Usim4(:,j))
    hold on
end
xlabel('time (ms)')
ylabel('voltage (mV)')
title('Neuron 4 voltage plot')
% legend(Usim1_Usim2_names)
hold off

% Plot the sensor output
subplot(6,1,4)
for j=1:Nseg
    plot(tSim,Isens(:,j))
    hold on
end
xlabel('time (ms)')
ylabel('Current (nA)')
title('Sensor Current plot')
% legend(Usim1_Usim2_names)
hold off







% Plot change in segemnt height per time step (basically speed)
figure
% % subplot(7,1,6)
% for j=1:Nseg
%     plot(tSim,dh(:,j))
%     hold on
% end
% xlabel('time (ms)')
% ylabel('dh (cm)')
% title('dh (speed) plot')
% % legend(h_dh_names)
% hold off
% 
% % Plot segment height at each time step
% subplot(7,1,7)
% for j=1:Nseg
%     plot(tSim,h(:,j))
%     hold on
% end
% xlabel('time (ms)')
% ylabel('h (cm)')
% title('segment height plot')
% % legend(h_dh_names)
% hold off

% Plot CPG U1 neuron potentials for each time step
subplot(Nseg+2,1,1)
scale_factor = 100; %pick one so you can see dh overlayed on CPG signals
for j=1:Nseg
    % plot(tSim,Usim1(:,j),tSim,Usim2(:,j),tSim,dh(:,j)*scale_factor,'--k')
%     plot(tSim,Usim1(:,j),tSim,Usim2(:,j))
    plot(tSim,Usim1(:,j))
    hold on
end
%xlabel('time (ms)')
ylabel('voltage (mV)')
title('CPG U1')
% legend(Usim1_Usim2_names)
hold off

% Plot CPG U2 neuron potentials for each time step
subplot(Nseg+2,1,2)
scale_factor = 100; %pick one so you can see dh overlayed on CPG signals
for j=1:Nseg
    % plot(tSim,Usim1(:,j),tSim,Usim2(:,j),tSim,dh(:,j)*scale_factor,'--k')
%     plot(tSim,Usim1(:,j),tSim,Usim2(:,j))
    plot(tSim,Usim2(:,j))
    hold on
end
% xlabel('time (ms)')
ylabel('voltage (mV)')
title('CPG U2')
% legend(Usim1_Usim2_names)
hold off


%color_mat = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
for j=1:Nseg
    subplot(Nseg+2,1,j+2)
%     plot(tSim,h(:,j),'color',color_mat(j,:))
    plot(tSim,h(:,j))
%     xlabel('time (ms)')
    ylabel('h (cm)')
    title(['Segment ' num2str(j)])
end
subplot(Nseg+2,1,Nseg+2)
xlabel('time (ms)')
% title('Segment Heights')
% 
% % Plot other neuron potentials for each time step
% subplot(5,1,3)
% for j=1:Nseg
%     plot(tSim,Usim3(:,j))
%     hold on
% end
% xlabel('time (ms)')
% ylabel('voltage (mV)')
% title('Neuron 3 voltage plot')
% % legend(Usim1_Usim2_names)
% hold off
% 
% subplot(5,1,4)
% for j=1:Nseg
% 
%     plot(tSim,Usim4(:,j))
%     hold on
% end
% xlabel('time (ms)')
% ylabel('voltage (mV)')
% title('Neuron 4 voltage plot')
% % legend(Usim1_Usim2_names)
% hold off
% 
% % Plot the sensor output
% subplot(5,1,5)
% for j=1:Nseg
%     plot(tSim,Isens(:,j))
%     hold on
% end
% xlabel('time (ms)')
% ylabel('Current (nA)')
% title('Sensor Current plot')
% % legend(Usim1_Usim2_names)
% hold off
% 
% 
% 
% 
% 
% 
