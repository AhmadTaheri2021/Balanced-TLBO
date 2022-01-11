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

function [BestCosts,BestSolCost]=BTLBO_Algorithm(D,MaxFEs,LB,UB,Population,nPop,CostFunction)

%------------------------------------
% Balanced Teaching Learning Based Optimization (BTLBO)
% The tutorial video for BTLBO can be found on youtube:
% https://youtu.be/OubHU3xe_Kc
%------------------------------------
%% Problem Definition
VarSize = [1 D]; %  Variables Matrix Size
%LB = LB.*ones(1,D);       %   Lower Bound
%UB = UB.*ones(1,D);       %   Upper Bound
%% Initialization 
% Initialize Population 
 pop = Population(1:nPop);
 [~,SortedIndx] = sort([pop.Cost]);
 pop = pop(SortedIndx(1:nPop));
% Initialize Best Solution
BestSol = pop(1);
% Select Teacher
Teacher = pop(1);
% Initialize Best Cost Record
BestCosts = ones(MaxFEs,1)+realmax;
BestCosts(1) = BestSol.Cost;
%% BTLBO Parameters
  FC = ones(1,nPop).* 0;
  Lambda = (nPop/100)^2;
  MaxFail =(MaxFEs/(20*nPop))*exp((1/MaxFEs)^Lambda);
  
  Rho = 0;
  FEs = 1;
%% BTLBO Main Loop
while FEs < MaxFEs 
    % Calculate the population Mean M // The mean value of all decision variables  
     Mean = zeros(1,D);
    for i=1:nPop      
        Mean = Mean + pop(i).Position;
    end
    Mean = Mean/nPop;   
 
    % Identify Teacher // Teacher is the best solution in population
      [~ , sortedIndx] = sort([pop.Cost]);
       Teacher = pop(sortedIndx(1));
 for i= 1 : nPop 

     alpha = FEs/MaxFEs - Rho;
     
   %  Choose Xj where j <> i    
        A = 1:nPop;
        A(i)=[];
        j = A(randi(nPop-1));

     Xnew = pop(i);
     
     %   Choose a Phase randomly from   1:Teaching Phase  2:Learning Phase  3:Tutoring Phase
     Ph = randi([1 3]);
 
%% Teaching Phase 1
 if Ph == 1
   % according to Eq. (8), Eq. (9), Eq. (10)
   % TF = randi([1 2]); 
   S = randi([0 1]);   
   TF =   randi([1 2]);
  Mw =  (Mean + (( S .* pop(j).Position + (1-S) .* pop(i).Position)))./2;   
  Xnew.Position = pop(i).Position ...
            +  (S).*rand.*1.*(TF).*(Teacher.Position - ( Mw ))...
            + (1-S).*rand.* (Teacher.Position - (TF).*( Mw ));
    
 end
 
%% Learning Phase 2
  if Ph == 2
   %  according to Eq. (4), Eq. (5)      
     Step = ( pop(i).Position - pop(j).Position);
        if pop(j).Cost < pop(i).Cost
            Step =  -  Step;
        end
        
        Xnew.Position = pop(i).Position +  rand(VarSize) .* Step;         
  end  % Phase 2      
  
%% Tutoring Phase 
  if  Ph == 3        
       % according to Eq. (11), Eq. (12), Eq. (13)          
       Indxs = randperm(D, ceil(D*min(1 * rand * exp(-(1-alpha))^2 ,1)));% 
       Xnew.Position(Indxs) = BestSol.Position(Indxs) + (1-(FEs/MaxFEs)).^2 .* unifrnd(-1,1,1,size(Indxs,2)) .*((pop(i).Position(Indxs) - pop(j).Position(Indxs)));      
  end  % Phase 3
                
    % Bound constraints control based on Eq.(17)
     LastPos = pop(i).Position;      
    Xnew = Clipping(LB,UB,Xnew,LastPos);
    % Evaluating Xi,new
    Xnew = Eval(Xnew,CostFunction);
       FEs = FEs + 1;   
       % Comparision
        if Xnew.Cost<pop(i).Cost
              pop(i) = Xnew;
              FC(i) = 0;
           
              if Xnew.Cost < Teacher.Cost
                Teacher = pop(i);
               end  
              
              if Xnew.Cost < BestSol.Cost
                BestSol = Xnew;
              end    
        else
           FC(i) = FC(i) + 1;  
           
        end    

    % 
    BestCosts(FEs) = BestSol.Cost;
 %   DiversityRTLBO_(FEs) =std([pop.Position]);% std([pop.Cost])/mean([pop.Cost]);
 
 end % nPop for 
  
 %% --------------- Restarting Phase ---------
  % 
 FT = (MaxFail * (0.005+((alpha)).^2)); 
       
  if size(find(FC >= (MaxFail/10) * (1+(FEs/MaxFEs).^(1-FEs/MaxFEs))),2) > nPop/2   
      % Restart 
           Rho = (FEs/MaxFEs)/2; 
    for i2 = 1 : nPop 
      pop(i2).Position = LB + rand(1,D).*(UB - LB);
      % Evaluation 
       pop(i2) = Eval(pop(i2),CostFunction);
       FEs = FEs + 1;
      % 
             FC(i2) = 0;  
            if pop(i2).Cost < BestSol.Cost
                 BestSol = pop(i2);
                 Teacher = pop(i2);
            end
    % 
    BestCosts(FEs) = BestSol.Cost;

    end
  else
     % Individual Restarting
   for i2 = 1 : nPop  
   if FC(i2) >= FT
      K = ceil(D*((rand * (exp(-(1-alpha)).^2))));
      Indxs = randperm(D,K);
      pop(i2).Position(Indxs) = LB(Indxs) + rand(1,size(Indxs,2)).*(UB(Indxs) - LB(Indxs));      
  
       % Evaluation 
    pop(i2) = Eval(pop(i2),CostFunction);
       FEs = FEs + 1;
      % 
             FC(i2) = 0;  
            if pop(i2).Cost < BestSol.Cost
                 BestSol = pop(i2);
                 Teacher = pop(i2);
            end
    % 
    BestCosts(FEs) = BestSol.Cost;
   end
   end
   
  end
            
 %% 
  formatSpec = '%10.7e';
  disp(['BTLBO  FEs ' num2str(FEs) ':  Best Cost = ' num2str(BestCosts(FEs),formatSpec) ' Teacher = ' num2str(Teacher.Cost,formatSpec)]);
 
end
% Return Xbest  
BestSolCost=BestSol.Cost;
end

function [newsol] = Clipping(LB,UB,newsol,LastPos)
% Bound constraints control based on Eq.(17)
      [~,underLB] = find(newsol.Position < LB);
      [~,uperUB] = find(newsol.Position > UB);
      if ~isempty(underLB)
      newsol.Position(underLB) =  LastPos(underLB) + unifrnd(0,1,1,size(underLB,2)).*  ( LB(underLB) - LastPos(underLB));
      end
      if ~isempty(uperUB)        
      newsol.Position(uperUB) =  LastPos(uperUB) + unifrnd(0,1,1,size(uperUB,2)).* ( UB(uperUB) - LastPos(uperUB)); 
      end
      
end

function [Xnew] = Eval(Xnew,CostFunction)
     % Evaluation
    %   switch FuncType
    %      case 1
     %   Xnew.Cost = CostFunction(Xnew.Position);
     %         case 2
    %    Xnew.Cost = benchmark_func(Xnew.Position,CostFunction) - Opt; % 2005
    %          case 3
      Xnew.Cost= feval('cec14_func',Xnew.Position',CostFunction)  - (CostFunction*100); % CEC2014 F(X) - F(X*)
    %    Xnew.Cost = YourCostFunc(Xnew.Position);
    %  end
            
end
