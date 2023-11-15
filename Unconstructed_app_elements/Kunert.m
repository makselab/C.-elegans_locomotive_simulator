function neuron_state_evolution = Kunert(noisy_connectome,nc,norm_chem,norm_gap,x0,W_neg_chem,W_pos_chem,W_gap,frequencies,drivers,std_,chem_switch,gap_switch,noise_switch,total_time,dt,V_Osc,V_drive,V_thresh,synvarstd,synvarmean,seed)
% MODEL - Kunert et al., 2014
addpath functions
Nn=size(x0,1);

%% Introduce variations to the connections
if noisy_connectome==1
    W_pos_chem = W_pos_chem + (W_pos_chem>0).*normrnd(0, nc, [Nn, Nn]);
    W_neg_chem = W_neg_chem + (W_neg_chem>0).*normrnd(0, nc, [Nn, Nn]);
    temp = (W_gap>0).*normrnd(0, nc/(sqrt(2)), [Nn, Nn]);%careful with adding two normal distributions 
    W_gap = W_gap + (temp+temp');
end

%% Normalize the total in-degree of each neuron
if norm_chem==1
    W_pos_chem = bsxfun(@rdivide, W_pos_chem, sum((W_pos_chem),1)); W_pos_chem(isnan(W_pos_chem)) = 0;
    W_neg_chem = bsxfun(@rdivide, W_neg_chem, sum((W_neg_chem),1)); W_neg_chem(isnan(W_neg_chem)) = 0;
end
if norm_gap==1
    W_gap = bsxfun(@rdivide, W_gap, sum((W_gap),1)); W_gap(isnan(W_gap)) = 0;
end

%% CONSTANTS & MEMORY

ar = 1.0; % Synaptic activity's rise time
ad = 5.0; % Synaptic activity's decay time
seq = ar./(ar+(2*ad));
dt_sqrt = sqrt(dt);% sqaure root of time 
dt_frac = power(dt,3/4);
glS_c = 10;% time constant G^c/C; kunert value 100
Vr = -35;%V resting potential
chem_tc = 100; % chem membrane time constant g/C; kunert value 200
gap_tc = 100; % chem membrane time constant g/C; kunert value 10
Vs_act = 0;% V reversal potential for activating connections
Vs_rep = -48;% V reversal potential for inhibiting connections
gamma =  0.125;
time = 0;% inital time value
in_current_const = 1000;% 1pA/1pF
%drivers are the neurons driving the system through an input
%frequencies are the frequencies of these neurons driving the system

s0 = (normrnd(synvarmean, synvarstd, [Nn, 1]));
simtime = 0:dt:total_time-dt;
duration = length(simtime);
neuron_state_evolution = nan(Nn,duration);
%neuron_state_evolution = nan(Nn*3,duration); % kunert testing
neuron_state_evolution(:,1) = x0;
N_k4 = normrnd(0, std_, [Nn, 1]);% create initial gaussian noise

%% potential threshold calculation
if V_thresh==1
    delta = -24*ones([Nn 1]);%V threshold potential
elseif V_thresh==0
    A = eye(Nn) + (gap_tc/glS_c)*diag(gap_switch*sum(W_gap,1)) + (chem_tc/glS_c)*diag(chem_switch*seq.*sum((W_neg_chem+W_pos_chem),1)) - (gap_tc/glS_c)*(W_gap - diag(diag(W_gap)))*gap_switch;
    b = Vr*ones([Nn 1]) + chem_switch*Vs_act*(chem_tc/glS_c)*seq.*sum((W_pos_chem),1)' + chem_switch*Vs_rep*(chem_tc/glS_c)*seq.*sum((W_neg_chem),1)' + (in_current_const/glS_c)*V_drive*abs(drivers);
    delta = A\b;
end
save('kunert_jacobian','delta')

%s0 = (ar*(1./(1+exp(-gamma*(x0-delta)))))./((ar*(1./(1+exp(-gamma*(x0-delta))))) + ad);
%neuron_state_evolution(:,1) = [x0; s0; delta]; % kunert testing

%% FUNCTIONS
driver_update=@(time) in_current_const*(V_drive.*abs(drivers) + V_Osc*sin(frequencies*time).*drivers);%
noise_update=@(Noise) in_current_const*Noise.*abs(drivers);%holding noise contant for each loop (not the best concept)
synaptic_update=@(x0,s0) (ar*(1./(1+exp(-gamma*(x0-delta)))).*(1-s0))-(ad.*s0);
S=@(x0,s0) (x0-Vs_act).*((s0'*W_pos_chem))' + (x0-Vs_rep).*((s0'*W_neg_chem))';% chemical effects
Gap=@(x0) sum(W_gap.*((x0-x0')),2);% gap junction effects
neuron_update=@(x0,s0) -glS_c.*(x0-Vr) - chem_switch*S(x0,s0)*chem_tc - gap_switch*Gap(x0)*gap_tc;% neuron model

%% Random walk generator
if seed>0
    rng(seed,'philox')
end
if noise_switch~=0
    %RW = nan(Nn,duration);
    RW = nan(1,duration);
    %RW(:,1) = drivers.*normrnd(0, std_, [Nn, 1]);
    RW(1) = normrnd(0, std_);
    for i = 2:duration
        %RW(:,i) = RW(:,i-1)+(drivers.*normrnd(0, std_, [Nn, 1])) ;
%        RW(i) = RW(i-1)+normrnd(0, std_) ;
        if ((abs(RW(i-1))/100)>rand)==0
            RW(i) = RW(i-1)+normrnd(0, std_);
        else
            RW(i) = RW(i-1)-(abs(normrnd(0, std_))*sign(RW(i-1)));
        end
    end
    % normalize between specifiec strength
    RW = RW.*abs(drivers);
    RW = (((RW - (min(RW,[],2)))./(max(RW,[],2)-min(RW,[],2)))*(2*noise_switch))-(noise_switch);
    RW(isnan(RW)) = 0;
    N_k4=RW(:,1);
end

%% MAIN SIMULATION LOOP
% Runge-kutta-4 update scheme
for t = 1:duration-1
    %use two noise values to get the average of these two for
    %intermediate steps k2 and k3
    if noise_switch~=0
        N_k1 = N_k4;
        %N_k4 = normrnd(0, std_, [Nn, 1]);
        N_k4=RW(:,t+1);
        N_k2_k3 = (N_k4 + N_k1) / 2;%still considered white noise (dcov = 0)
    else
        N_k1 = zeros([Nn, 1]); N_k2_k3 = zeros([Nn, 1]); N_k4 = zeros([Nn, 1]);
    end

    %runge-kutta slopes
    %STOCHASTIC EFFECTS IN PHYSICAL SYSTEMS by
    %MAXI SAN-MIGUEL and RAUL TORAL equtation: 2.75 (stochastic ODEs)
    x0_k1 = neuron_update(  x0,s0)*dt + driver_update(time)*dt + noise_update(N_k1)*dt_sqrt;
    s0_k1 = synaptic_update(x0,s0)*dt;

    %
    x0_k2 = neuron_update(  x0 + x0_k1/2, s0 + s0_k1/2)*dt + driver_update(time + (dt/2))*dt + noise_update(N_k2_k3)*dt_frac;
    s0_k2 = synaptic_update(x0 + x0_k1/2, s0 + s0_k1/2)*dt;
    %

    x0_k3 = neuron_update(  x0 + x0_k2/2, s0 + s0_k2/2)*dt + driver_update(time + (dt/2))*dt + noise_update(N_k2_k3)*dt_frac;
    s0_k3 = synaptic_update(x0 + x0_k2/2, s0 + s0_k2/2)*dt;
    %

    x0_k4 = neuron_update(  x0 + x0_k3,   s0 + s0_k3)*dt + driver_update(time +  dt   )*dt + noise_update(N_k4)*dt;
    s0_k4 = synaptic_update(x0 + x0_k3,   s0 + s0_k3)*dt;
    %

    x0 = x0 + (x0_k1 + (2* x0_k2) + (2* x0_k3) + x0_k4)/6;
    s0 = s0 + (s0_k1 + (2* s0_k2) + (2* s0_k3) + s0_k4)/6;

    %Euler's  method
%     x0 = x0 + x0_k1;
%     s0 = s0 + s0_k1;
    neuron_state_evolution(:, t+1) = x0;
    %neuron_state_evolution(:, t+1) = [x0; s0; delta];% kunert testing
    time = time + dt;
end



end