%{
NOTE:
This is a script to simulate time traces for dCov_nonlinear analysis
%}

clear; close all
cd('~/workspace/DDC/')

%--------------------------------------------
% 3 NEURON SIMULATION
%--------------------------------------------

savedir = '~/workspace/DDC/out/';
NTrial = 1;

%% Simulation parameters
% linear simulation
deltaT = 0.01;
Ttotal = 10; %default 1000 
TR = deltaT;
params.deltaT = deltaT;
params.TR = TR;
params.NTrial = NTrial;
params.Ttotal = Ttotal;
T = round(Ttotal/deltaT); N = 3;
u = zeros(T,N); 
u(T/4:T/4*3,1) = 1;

S_list = -[0.1 0.2 0.5 0.8];
RandInput_list = [0.01 0.1 1 10 20 30 40 50 60 70 80 90 100];
V_pre_list = {}; 
for Input_idx = 1:length(RandInput_list)
	tmp = {};
	for S_idx = 1:length(S_list)
		V_pre_multi = zeros(Ttotal/deltaT,3,NTrial);
		for Trial = 1:NTrial
			S = S_list(S_idx);
			RandInput = RandInput_list(Input_idx);
			G_confounder = -1*eye(3);
			G_confounder(2,1) = S; G_confounder(3,1) = S;
			G_chain = -1*eye(3);
			G_chain(2,1) = S; G_chain(3,2) = S;
			G = G_confounder;
			% G = G_chain;
			% V_pre = Linear_simulation(G,deltaT,RandInput);
			V_pre = Linear_simulation(G,deltaT,RandInput,Ttotal,u);
			V_pre_multi(:,:,Trial) = V_pre;
		end
		tmp{S_idx} = V_pre_multi;
	end
	V_pre_list(Input_idx,:) = tmp;
	disp(['Progress: ' num2str(Input_idx) '/' num2str(length(RandInput_list))])
end
% save([savedir 'Confounder_Ts.mat'],'V_pre_list','S_list','RandInput_list','params','-v7.3')
% save([savedir 'Chain_Ts.mat'],'V_pre_list','S_list','RandInput_list','params','-v7.3')

save([savedir 'Confounder_HalfON_Ts.mat'],'V_pre_list','S_list','RandInput_list')

% Calculating DDC for V_pre_list S=0.5 and RandInput=0.5:
[Cov,precision,B,dCov] = estimators(V_pre_list{3,3},0,deltaT)
DDC=dCov*inv(B)