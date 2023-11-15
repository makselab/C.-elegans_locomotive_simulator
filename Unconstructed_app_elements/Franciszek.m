function neuron_state_evolution = Franciszek(noisy_connectome,nc,norm_chem,norm_gap,x0,W_neg_chem,W_pos_chem,W_gap,frequencies,drivers,std_,chem_switch,gap_switch,noise_switch,total_time,dt,V_Osc,V_drive,V_thresh,seed)
% MODEL - Synaptic polarity of the interneuron circuit controlling C. elegans locomotion, Franciszek
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
dt_sqrt = sqrt(dt);% sqaure root of time
dt_frac = power(dt,3/4);
glS_c = 10;% time constant g_l*S/C*S [1/s]
Vr = -35;%V resting potential
chem_tc = 100; % chem membrane time constant q_s/C = 0.1nS/1.5pF 
gap_tc = 100; % chem membrane time constant q_e/C  = 0.1nS/1.5pF 
Vs_act = 0;% V reversal potential for activating connections
Vs_rep = -48;% V reversal potential for inhibiting connections
gamma = 0.125; %1/0.15e-3V slope of potential 150.5
time = 0;% inital time value
in_current_const = 1000;% 1pA/1.5pF
%drivers are the neurons driving the system through an input
%frequencies are the frequencies of these neurons driving the system

simtime = 0:dt:total_time-dt;
duration = length(simtime);
neuron_state_evolution = nan(Nn,duration);
%neuron_state_evolution = nan(Nn*3,duration); % kunert testing
neuron_state_evolution(:,1) = x0;
N_k4 = normrnd(0, std_, [Nn, 1]);% create initial gaussian noise


if V_thresh==1
    delta = -24*ones([Nn 1]);%V threshold potential
elseif V_thresh==0
%     A = eye(Nn) + (gap_tc/glS_c)*diag(gap_switch*sum(W_gap,1)) + (1/2)*(chem_tc/glS_c)*diag(chem_switch*sum((W_neg_chem+W_pos_chem),1)) - (gap_tc/glS_c)*(W_gap - diag(diag(W_gap)))*gap_switch;
%     b = Vr*ones([Nn 1]) + (1/2)*chem_switch*Vs_act*(chem_tc/glS_c)*sum((W_pos_chem),1)' + (1/2)*chem_switch*Vs_rep*(chem_tc/glS_c)*sum((W_neg_chem),1)' + (in_current_const/glS_c)*V_drive*abs(drivers);
%     delta = A\b;
    A = eye(Nn) + (gap_tc/glS_c)*diag(gap_switch*sum(W_gap,1)) + (chem_tc/glS_c)*diag(chem_switch*(1/2).*sum((W_neg_chem+W_pos_chem),1)) - (gap_tc/glS_c)*(W_gap - diag(diag(W_gap)))*gap_switch;
    b = Vr*ones([Nn 1]) + chem_switch*Vs_act*(chem_tc/glS_c)*(1/2).*sum((W_pos_chem),1)' + chem_switch*Vs_rep*(chem_tc/glS_c)*(1/2).*sum((W_neg_chem),1)' + (in_current_const/glS_c)*V_drive*abs(drivers);
    delta = A\b;
end
save('delta','delta')
%neuron_state_evolution(:,1) = [x0; s0; delta]; % kunert testing
%% FUNCTIONS
%S=@(x0) (x0-Vs_act).*(sum((1./(1+exp(-gamma*(x0-delta)))).*W_pos_chem,1))' + (x0-Vs_rep).*(((1./(1+exp(-gamma*(x0-delta))))'*W_neg_chem))';% chemical effects
S=@(x0) (x0-Vs_act).*(((1./(1+exp(-gamma*(x0-delta))))'*W_pos_chem))' + (x0-Vs_rep).*(((1./(1+exp(-gamma*(x0-delta))))'*W_neg_chem))';% chemical effects
Gap=@(x0) sum(W_gap.*((x0-x0')),2);% gap junction effects
neuron_update=@(x0) -glS_c.*(x0-Vr) - chem_switch*S(x0)*chem_tc - gap_switch*Gap(x0)*gap_tc;% model based on Franciszek
driver_update_osc=@(time) in_current_const*(V_Osc*sin(frequencies*time).*drivers);% driving force %negative in_current produces osciallations
driver_update_con = in_current_const*(V_drive.*abs(drivers));
noise_update=@(Noise) in_current_const*Noise.*abs(drivers);%noise

%% Random walk generator
if seed>0
    rng(seed,'philox')
end
if noise_switch~=0
    %RW = nan(Nn,duration);%
    %RW(:,1) = drivers.*normrnd(0, std_, [Nn, 1]);%
    RW = nan(1,duration);
    RW(1) = normrnd(0, std_);
    for i = 2:duration
        %RW(:,i) = RW(:,i-1)+(drivers.*normrnd(0, std_, [Nn, 1])) ;
        %RW(i) = RW(i-1)+normrnd(0, std_) ;
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
    % MAXI SAN-MIGUEL and RAUL TORAL EQ: 2.75 (stochastic ODEs)
%     if t<15000
%     x0_k1 = neuron_update(x0)*dt + driver_update_con*dt;
%     x0_k2 = neuron_update(x0 + x0_k1/2 )*dt + driver_update_con*dt;
%     x0_k3 = neuron_update(x0 + x0_k2/2 )*dt + driver_update_con*dt;
%     x0_k4 = neuron_update(x0 + x0_k3   )*dt + driver_update_con*dt;
%     else
    x0_k1 = neuron_update(x0)*dt + driver_update_con*dt + driver_update_osc(time)*dt + noise_update(N_k1)*dt_sqrt;
    x0_k2 = neuron_update(x0 + x0_k1/2 )*dt + driver_update_con*dt + driver_update_osc(time + (dt/2))*dt + noise_update(N_k2_k3)*dt_frac;
    x0_k3 = neuron_update(x0 + x0_k2/2 )*dt + driver_update_con*dt + driver_update_osc(time + (dt/2))*dt + noise_update(N_k2_k3)*dt_frac;
    x0_k4 = neuron_update(x0 + x0_k3   )*dt + driver_update_con*dt + driver_update_osc(time +  dt   )*dt + noise_update(N_k4)*dt;
%     end
    x0 = x0 + (x0_k1 + (2* x0_k2) + (2* x0_k3) + x0_k4)/6;
    
    %neuron_state_evolution(:,t+1) = [x0; s0; delta]; % kunert testing
    neuron_state_evolution(:, t+1) = x0;
    time = time + dt;
end

end