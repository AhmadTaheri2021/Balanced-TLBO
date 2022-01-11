%%-----------------------------------------------------------------------------%
% 
%  Copyright 2021 Ahmad Taheri. All Rights Reserved.
%   
%  Licensed under the Apache License, Version 2.0 (the "License");
%  you may not use this file except in compliance with the License.
%  You may obtain a copy of the License at
%      http://www.apache.org/licenses/LICENSE-2.0
%    
%  Unless required by applicable law or agreed to in writing, software
%  distributed under the License is distributed on an "AS IS" BASIS,
%  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%  See the License for the specific language governing permissions and
%  limitations under the License.
%
%%-----------------------------------------------------------------------------%

%------------------------------------
% Balanced Teaching Learning Based Optimization (BTLBO)
% The tutorial video for BTLBO can be found on youtube:
% https://youtu.be/OubHU3xe_Kc
%------------------------------------
clear;
clc;
%%            settings 

   nPop=50;         % Population Size
   D = 10;               % number of Decision variables 
   MaxFEs= 200000;       % Maximum number of function evaluations
   NumofExper = 1;   % Number of test
  % Benchmark = 3;  % 1:Classic 2:CEC2005 3:CEC2014 
   Func_id = 1; % CEC2014 1~30 / CEC2005 1~25 / Classic 1~22
  % FileName = [ FileName ' D30_NFE300K_30t ']; 

% =====================================================================================

 Opt = zeros(1,30); % 
 global initial_flag
 initial_flag = 0;
%% 
Function_name=['F' num2str(Func_id)];
%========== CEC2014 ==========
fhd=str2func('cec14_func');
CostFunction=Func_id;
LB=-100;%lb;
UB=100;%ub;
Opt = 100:100:3000;
%================
 
LB = LB.*ones(1,D);       %   Lower Bound
UB = UB.*ones(1,D);       %   Upper Bound

% Empty Solution Structure
empty_Solution.Position=[];
empty_Solution.Cost=[];

% Initialize Harmony Memory
Population=repmat(empty_Solution,nPop,1);


%% TLBO
SumBestCostBTLBO_=zeros(MaxFEs,1);
BestSolCostBTLBO= []; %zeros(MaxFEs,1);
SumBestCostTLBO_=zeros(MaxFEs,1);
BestSolCostTLBO= []; %zeros(MaxFEs,1);

%===================================================

for ii=1:NumofExper
    
 
  rand('state',sum(100*clock));
  initial_flag = 0; % should set the flag to 0 for each run, each function
  
% Create Initial Population
for i=1:nPop
   
    Population(i).Position= LB+rand(1,D).*(UB-LB); %LB+rand(1,D).*(UB-LB);
    
    Population(i).Cost = feval(fhd,Population(i).Position',CostFunction) -  (CostFunction*100); % CEC2014 F(X) - F(X*)
   %  Population(i).Cost = YourCostFunc(Population(i).Position);
end  
    
%%----------------------------------
tic;
%[BestCostRTLBO_,BestSolCostBTLBO(ii)]=TLBO_Main(D,MaxFEs,LB,UB,Population,nPop,CostFunction);  
[BestCostRTLBO_,BestSolCostBTLBO(ii)]=BTLBO_Algorithm(D,MaxFEs,LB,UB,Population,nPop,CostFunction);  
SumBestCostBTLBO_=SumBestCostBTLBO_+ BestCostRTLBO_(1:MaxFEs);
T2 = toc;
end


AveBestCostTLBO_=SumBestCostTLBO_ ./ NumofExper;
AveBestCostBTLBO_=SumBestCostBTLBO_ ./ NumofExper;
%% 
% Mean(1,1) = mean(BestSolCostTLBO);
Mean(1,2)=mean(BestSolCostBTLBO);

% SD(1,1)=std(BestSolCostTLBO);
SD(1,2)=std(BestSolCostBTLBO);

%filename=['BTLBO Result ' FileName ' _' Function_name '.mat']
%save(filename);
f1=figure;
semilogy(AveBestCostBTLBO_,'r-','LineWidth',2);
%hold on;
%semilogy(AveBestCostTLBO_,'b-','LineWidth',2);
grid on;
hold off;
xlabel('FEs');
str=['F(x) = ' Function_name];
ylabel(str);
%legend('BTLBO','TLBO');
legend('BTLBO');
