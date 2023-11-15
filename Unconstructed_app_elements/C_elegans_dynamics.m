function activities = C_elegans_dynamics(noisy_connect,ns,norm_indegree_chem,norm_indegree_gap,expression_type,diagonal_symmetry,off_diagonal_symmetry,intra_connections,binary_system,random_removal,noise,chem_strength,gap_strength,V_Osc,V_drive,voltages,SynapVar,rand_initial,V_thresh,V_mean,V_std,total_time,synvarstd,synvarmean,synch_measure,seed)
addpath functions

%% CONSTANTS & INITIAL CONDITIONS
worm = Organizer(expression_type,diagonal_symmetry,off_diagonal_symmetry);
Nodes = vertcat(worm.NeuronNames.NodesF,worm.NeuronNames.NodesB);
%% EXPRESSION TYPE
if expression_type=="Chem"
    input_neurons = sign(voltages(1,[34,36,37,35,31,32,33,38,39,11,12,13,15,16,53,52,51,50,55,54,14,30,10,29,17,18,19,20,21,1,2,3,4,5,8,9,22,27,24,25,28,26,6,7,49,48,43,42,45,44,47,46,41,40,23]));%chem+plus+neurons+with+no+conn
    frequencies = abs(voltages(1,[34,36,37,35,31,32,33,38,39,11,12,13,15,16,53,52,51,50,55,54,14,30,10,29,17,18,19,20,21,1,2,3,4,5,8,9,22,27,24,25,28,26,6,7,49,48,43,42,45,44,47,46,41,40,23]));%chem+plus+neurons+with+no+conn [Hz]
    Nn = size(worm.Matrices.W_ff_chem,1)+size(worm.Matrices.W_bb_chem,1); % total population size
    Nf = size(worm.Matrices.W_ff_chem,1); % forward population size
elseif expression_type=="Gap"
    input_neurons = sign(voltages(1,[14,15,16,19,20,13,30,10,32,12,33,11,29,34,35,31,36,37,55,54,51,50,53,52,18,19,22,23,24,25,26,27,28,3,6,7,5,9,8,4,1,17,2,20,21,45,44,47,46,49,48,41,40,43,42]));%gap+plus+neurons+with+no+conn
    frequencies = abs(voltages(1,[14,15,16,19,20,13,30,10,32,12,33,11,29,34,35,31,36,37,55,54,51,50,53,52,18,19,22,23,24,25,26,27,28,3,6,7,5,9,8,4,1,17,2,20,21,45,44,47,46,49,48,41,40,43,42]));%gap+plus+neurons+with+no+conn [Hz]
    Nn = size(worm.Matrices.W_ff_gap,1)+size(worm.Matrices.W_bb_gap,1); % total population size
    Nf = size(worm.Matrices.W_ff_gap,1); % forward population size
end

%% Parameters
omega = 2*pi*frequencies; % angular frequency

% set trial duration
dt = 0.00005;%time step size

transient_steps = 1; %number of steps from witch to begin plotting
simtime = 0:dt:total_time-dt;
if rand_initial == 1
    x0 = normrnd(V_mean, V_std, [Nn, 1]);%random dis around inactive voltage
elseif rand_initial == 0
    x0 = ones([Nn,1])*V_mean;%all inital volages 
end

% noise
std_ = 1;

%% RANDOM REMOVAL OF LINKS
if random_removal>0
    [row,col] = find(worm.Matrices.W_ff_chem);%forward
    removals = randperm(size(row,1),random_removal);
    for i = 1:random_removal
        worm.Matrices.W_ff_chem(row(removals(i)),col(removals(i)))=0;
        disp (string(worm.NeuronNames.NodesF(row(removals(i))))+"➜"+string(worm.NeuronNames.NodesF(col(removals(i)))))
    end
    [row,col] = find(worm.Matrices.W_bb_chem);%forward
    removals = randperm(size(row,1),random_removal);
    for i = 1:random_removal
        worm.Matrices.W_bb_chem(row(removals(i)),col(removals(i)))=0;
        disp (string(worm.NeuronNames.NodesB(row(removals(i))))+"➜"+string(worm.NeuronNames.NodesB(col(removals(i)))))
    end
    [row,col] = find(worm.Matrices.W_ff_gap);%forward
    removals = randperm(size(row,1),random_removal);
    for i = 1:random_removal
        worm.Matrices.W_ff_gap(row(removals(i)),col(removals(i)))=0;
        worm.Matrices.W_ff_gap(col(removals(i)),row(removals(i)))=0;
        disp (string(worm.NeuronNames.NodesF(row(removals(i))))+"⇆"+string(worm.NeuronNames.NodesF(col(removals(i)))))
    end
    [row,col] = find(worm.Matrices.W_bb_gap);%forward
    removals = randperm(size(row,1),random_removal);
    for i = 1:random_removal
        worm.Matrices.W_bb_gap(row(removals(i)),col(removals(i)))=0;
        worm.Matrices.W_bb_gap(col(removals(i)),row(removals(i)))=0;
        disp (string(worm.NeuronNames.NodesB(row(removals(i))))+"⇆"+string(worm.NeuronNames.NodesB(col(removals(i)))))
    end
end

%% GAP CONNECTIONS
W_gap = zeros(max(size(Nodes)),max(size(Nodes)));
W_gap(1:Nf,1:Nf) = worm.Matrices.W_ff_gap;
W_gap(Nf+1:Nn,Nf+1:Nn) = worm.Matrices.W_bb_gap;
if intra_connections==1
    W_gap(1:Nf,Nf+1:Nn) = worm.Matrices.W_fb_gap;%GAP INTRA CONNECTIONS
    W_gap(Nf+1:Nn,1:Nf) = worm.Matrices.W_bf_gap;%GAP INTRA CONNECTIONS
elseif intra_connections==0
    W_gap(1:Nf,Nf+1:Nn) = worm.Matrices.W_fb_gap*0;%GAP NO INTRA CONNECTIONS
    W_gap(Nf+1:Nn,1:Nf) = worm.Matrices.W_bf_gap*0;%GAP NO INTRA CONNECTIONS
end

if binary_system==1
    W_gap = double(W_gap>0);
end

%% CHEMICAL CONNECTIONS
W_pos_chem = zeros(max(size(Nodes)),max(size(Nodes)));
W_pos_chem(1:Nf,1:Nf) = worm.Matrices.W_ff_chem;
W_pos_chem(Nf+1:Nn,Nf+1:Nn) = worm.Matrices.W_bb_chem;
pos_chem_mask = double(sum(W_pos_chem,2)>0);

W_neg_chem = zeros(max(size(Nodes)),max(size(Nodes)));
W_neg_chem(Nf+1:Nn,1:Nf) = worm.Matrices.W_bf_chem; %INTRA CONNECTIONS
W_neg_chem(1:Nf,Nf+1:Nn) = worm.Matrices.W_fb_chem; %INTRA CONNECTIONS
neg_chem_mask = double(sum(W_neg_chem,2)>0);

if binary_system==1
    W_neg_chem = double(W_neg_chem>0);
    W_pos_chem = double(W_pos_chem>0);
end
if intra_connections==0
    W_neg_chem = W_neg_chem*0;
end
%% MASKS needed for intra-connections dynamics
joint_regulation_mask = neg_chem_mask.*pos_chem_mask;
inverse_joint_regulation_mask = double(~(joint_regulation_mask>0));
active_regulation_mask = pos_chem_mask.*inverse_joint_regulation_mask;
repressive_regulation_mask = neg_chem_mask.*inverse_joint_regulation_mask;

%% Color codes for lines based on expected simulation data
if expression_type=="Chem"
        range_f = [1:20]; Titulo1_f = "Forward Chemical"; nombres_f = worm.NeuronNames.NodesF_chem;
        color_codes_f = ['#2E9072'; '#2E9072'; '#FF5C81'; '#FF5C81'; '#C686E9'; '#C686E9'; '#C686E9'; '#C686E9'; '#C686E9'; '#C686E9'; '#C686E9'; '#C686E9'; '#C686E9'; '#C686E9'; '#5FC613'; '#5FC613'; '#00CAFF'; '#00CAFF'; '#EA8615'; '#EA8615'];

        range_b = [25:54]; Titulo1_b = "Backward Chemical"; nombres_b = worm.NeuronNames.NodesB_chem;
        color_codes_b = ['#A95A8D'; '#54D70C'; '#54D70C'; '#54D70C'; '#54D70C';'#00DCFF'; '#00DCFF'; '#00DCFF'; '#00DCFF'; '#9CA4FF'; '#9CA4FF';'#9CA4FF'; '#9CA4FF'; '#9CA4FF'; '#D6C6A1'; '#D6C6A1'; '#FF791E';'#FF791E'; '#FF791E'; '#FF791E'; '#005E6D'; '#005E6D'; '#FF79FF'; '#FF79FF'; '#000000'; '#000000'; '#D3A900'; '#D3A900'; '#00CC96'; '#00CC96'];

elseif expression_type=="Gap"
        range_f = [1:22]; Titulo1_f = "Forward Gap"; nombres_f = worm.NeuronNames.NodesF_gap;
        color_codes_f = ['#DF89FF'; '#DF89FF'; '#DF89FF'; '#DF89FF'; '#DF89FF'; '#75E8E6'; '#ff5584'; '#ff5584'; '#4c463e'; '#4c463e'; '#00BD94'; '#00BD94'; '#FF8805'; '#FF8805'; '#c5a9a7'; '#c5a9a7'; '#6ba613'; '#6ba613'; '#6ba613'; '#6ba613'; '#00c4ff'; '#00c4ff'];

        range_b = [25:53]; Titulo1_b = "Backward Gap"; nombres_b = worm.NeuronNames.NodesB_gap;
        color_codes_b = ['#E28BFF'; '#E28BFF'; '#E28BFF'; '#E28BFF'; '#E28BFF';'#E28BFF'; '#E28BFF'; '#E28BFF'; '#E28BFF'; '#E28BFF'; '#E28BFF';'#E28BFF'; '#d3b3b0'; '#73c000'; '#73c000'; '#73c000'; '#73c000';'#73c000'; '#73c000'; '#ff5584'; '#ff5584'; '#4c463e'; '#4c463e'; '#3ab7dd'; '#3ab7dd'; '#cfb42a'; '#cfb42a'; '#ff8805'; '#ff8805'];   
end

%% MAIN SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
if SynapVar==0
    activities = Franciszek(noisy_connect,ns,norm_indegree_chem,norm_indegree_gap,x0,W_neg_chem,W_pos_chem,W_gap,omega',input_neurons',std_,chem_strength,gap_strength,noise,total_time,dt,V_Osc,V_drive,V_thresh,seed);
else
    activities = Kunert(noisy_connect,ns,norm_indegree_chem,norm_indegree_gap,x0,W_neg_chem,W_pos_chem,W_gap,omega',input_neurons',std_,chem_strength,gap_strength,noise,total_time,dt,V_Osc,V_drive,V_thresh,synvarstd,synvarmean,seed);
end
save('activities','activities')

%activities is returned in mV

%% PLOTS
R = activities;
taxis = simtime(transient_steps:end)-(transient_steps-1)*dt;
idx = (1:Nn)';
linestyle={':','-.'};

%FORWARD DYNAMICS FIGURE OUTPUT (SAVED AND CALLED)
f1 = figure('visible','off');
for i=range_f
   lind = mod(i,2);
   plot(taxis(1,1:end),R(idx(i),:),linestyle{lind+1},'linew',3);
   hold on
   xlabel('time [Seconds]'); ylabel('Neuron state [mV]'); set(gca,'fontsize',17); 
end
colororder(color_codes_f)%colors of lines in plot
str = strings(size(nombres_f,1),1); j = 1;
for i = range_f
    if input_neurons(i)~=0
        str(j,1) = "\bf ";
    end
    j = j + 1;
end
legend(str+nombres_f,'fontsize',14,'Location','eastoutside');
if binary_system==1
    Titulo2 = split(diagonal_symmetry); Titulo2 = strjoin(Titulo2(1:end-1)," "); Titulo2 = Titulo2+"- Binarized";
elseif binary_system==0
    Titulo2 = diagonal_symmetry;
end
if norm_indegree_chem==1
    Titulo2 = Titulo2+" - Chem Normalized";
elseif norm_indegree_gap==1
    Titulo2 = Titulo2+" - Gap Normalized";
end

title({Titulo1_f,Titulo2,"Neruon States"}) 
saveas(f1,"forward_dynamics")

%BACKWARD DYNAMICS FIGURE OUTPUT (SAVED AND CALLED)
f2 = figure('visible','off');
for i=range_b
   lind = mod(i,2);
   plot(taxis(1,1:end),R(idx(i),:),linestyle{lind+1},'linew',3);
   hold on
   xlabel('time [Seconds]'); ylabel('Neuron state [mV]'); set(gca,'fontsize',17); 
end
colororder(color_codes_b)%colors of lines in plot
str = strings(size(nombres_b,1),1); j = 1;
for i = range_b
    if input_neurons(i)~=0
        str(j,1) = "\bf ";
    end
    j = j + 1;
end
legend(str+nombres_b,'fontsize',14,'Location','eastoutside');

if binary_system==1
    Titulo2 = split(diagonal_symmetry); Titulo2 = strjoin(Titulo2(1:end-1)," "); Titulo2 = Titulo2+" - Binarized";
elseif binary_system==0
    Titulo2 = diagonal_symmetry;
end
if norm_indegree_chem==1
    Titulo2 = Titulo2+" - Chem Normalized";
elseif norm_indegree_gap==1
    Titulo2 = Titulo2+" - Gap Normalized";
end

title({Titulo1_b,Titulo2,"Neruon States"}) 
saveas(f2,"backward_dynamics")

la_B_gap = [1 1 1 1 1 1 1 1 1 1 1 1 2 3 3 3 3 3 3 4 4 5 5 6 6 7 7 8 8];
la_F_gap = [1 1 1 1 1 2 3 3 4 4 5 5 6 6 7 7 8 8 8 8 9 9];
la_B_chem = [1 2 2 2 2 3 3 3 3 4 4 4 4 4 5 5 6 6 6 6 7 7 8 8 9 9 9 9 10 10];
la_F_chem = [1 1 2 2 3 3 3 3 3 3 3 3 3 3 4 4 5 5 6 6]; 

if expression_type=="Chem"
    last_assignement = la_F_chem;
elseif expression_type=="Gap"
    last_assignement = la_F_gap;
end

%FORWARD SYNCHS FIGURE OUTPUT (SAVED AND CALLED)
f3 = figure('visible','off');


time_window = 1;%sec
distance_matrix=zeros(Nn,Nn);
for i=1:Nn
        for j=1:Nn
            if ((i <= Nf) && (j <= Nf)) || ((i > Nf) && (j > Nf)) %no interest in cross synchrony
                if i==j
                    distance_matrix(i,j) = 1;
                else
                    distance_matrix(i,j) = measuremnet_types(activities(i,end-round(time_window/dt):end),activities(j,end-round(time_window/dt):end),synch_measure);
                end
            end
        end
end
synch_measure = distance_matrix;
%remove comment to get idealized LoS
% synch_measure(synch_measure>0.99)=1;
% synch_measure(synch_measure<0.99)=0;

data = synch_measure(range_f,range_f);
imagesc(data,'alphadata', ~(isinf(data)))
colorbar
hold on
size_commu = zeros(max(last_assignement)+1,2);

xticks(0:100:length(data));
yticks(0:100:length(data));
title({Titulo1_f,"Synchronicity Measure"}) 
ylabel('Neurons'); xlabel('Neurons'); set(gca,'fontsize',17);

aR=gca;
aR.XTick = 1:Nn;
aR.XTickLabel = Nodes(range_f)';
aR.XTickLabelRotation = 90;
aR.YTick = 1:Nn;
aR.YTickLabel = Nodes(range_f)';
aR.YTickLabelRotation = 0;
aR.XAxis.FontSize = 12;
aR.YAxis.FontSize = 12;
aR.GridAlpha=0.8;
aR.MinorGridAlpha = 1;
aR.CLim = [0 1];

idx_communities = (1:max(last_assignement))';

for j = 1:length(idx_communities)

for i = 1:size(data,1) % from 1 to max number of nodes
    
    comm = idx_communities(j);
    
    if last_assignement(i) == comm
        size_commu(j,1) = j;
        size_commu(j,2) = size_commu(j,2) + 1;
    else

    end
end
end

%%%%%%%%  BLOCKS   %%%%%%%%

hold on
orig = 0.5;
edge = size_commu(1,2);

for i = 1:(length(size_commu)-1)
    
    h(i) = rectangle('position',[orig orig edge edge]);
    set(h(i),'EdgeColor',[1 0 0],'linewidth',2);
    yline(orig);
    orig = orig + size_commu(i,2);
    edge = size_commu(i+1,2);
    
end
saveas(f3,"forward_synchs")

if expression_type=="Chem"
    last_assignement = la_B_chem;
elseif expression_type=="Gap"
    last_assignement = la_B_gap;
end

%BACKWARD SYNCHS FIGURE OUTPUT (SAVED AND CALLED)
f4 = figure('Visible','off');
data = synch_measure(range_b,range_b);
imagesc(data,'alphadata', ~(isinf(data)))
colorbar
hold on
size_commu = zeros(max(last_assignement)+1,2);

xticks(0:100:length(data));
yticks(0:100:length(data));
title({Titulo1_b,"Synchronicity Measure"})
ylabel('Neurons'); xlabel('Neurons'); set(gca,'fontsize',17);

aR=gca;
aR.XTick = 1:Nn;
aR.XTickLabel = Nodes(range_b)';
aR.XTickLabelRotation = 90;
aR.YTick = 1:Nn;
aR.YTickLabel = Nodes(range_b)';
aR.YTickLabelRotation = 0;
aR.XAxis.FontSize = 12;
aR.YAxis.FontSize = 12;
aR.GridAlpha=0.8;
aR.MinorGridAlpha = 1;
aR.CLim = [0 1];

idx_communities = (1:max(last_assignement))';

for j = 1:length(idx_communities)

for i = 1:size(data,1) % from 1 to max number of nodes
    
    comm = idx_communities(j);
    
    if last_assignement(i) == comm
        size_commu(j,1) = j;
        size_commu(j,2) = size_commu(j,2) + 1;
    else

    end
end
end

%%%%%%%%  BLOCKS   %%%%%%%%

hold on
orig = 0.5;
edge = size_commu(1,2);

for i = 1:(length(size_commu)-1)
    
    h(i) = rectangle('position',[orig orig edge edge]);
    set(h(i),'EdgeColor',[1 0 0],'linewidth',2);
    yline(orig);
    orig = orig + size_commu(i,2);
    edge = size_commu(i+1,2);
    
end
saveas(f4,"backward_synchs")
% f5 = figure('Visible','off');
% heatmap(synch_measure,'XData',Nodes,'YData',Nodes);
% saveas(f5,"heat")
load('delta');
sum(delta-activities(:,end),'all')
end