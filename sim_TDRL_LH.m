function R = sim_TDRL_LH_final(sim)

% This script simulates the experimental paradigms included in the paper by
%Sharpe et al. entitled:

% What’s past is prologue: prior experience shapes the neural
% circuits recruited for future learning.

    %
    % USAGE: results = sim_TDRL(sim)
    %
    % INPUTS:
    %   sim - string specifying simulation:
    %     '4A' -             Figure 4A: fear conditioning of single cue and context 
    %                        using a Mackintosh attentional rule
    %
    %     '4B - Low Alpha' - Figure 4B: conditioning as done in 4A. The context
    %                        alpha starts low due to prior appetitive training in
    %                        the conditioning box. Opto manipulation
    %                        results in modulation of the updates to the
    %                        eligibility traces
    %
    %     '4B - App Value' - Figure 4B: conditioning as done in 4A. The context
    %                        value starts out positively due to appetitive associations
    %                        from prior appetitive learning. Opto manipulation
    %                        results in modulation of the updates to the
    %                        eligibility traces 
    %                   
    %     '4C' -             Figure 4C: fear conditioning as in Figure 4A.
    %                        Contexts between appetitive and fear conditioning
    %                        are switched so neither the low alpha
    %                        or increase in appetitive value due to prior
    %                        conditioning applies. Opto manipulation
    %                        results in modulation of the updates to the
    %                        eligibility traces  
    %
    % OUTPUTS:
    %   R
    %
    %Matt Gardner April 2020
    


%This provides a transfer function for determining Pavlovian conditioned responses
%from value
CR = @(V,max,Voff,ks) max./(1 + exp(-ks*(-V - Voff))) - .1;
%parameters used, based on single cue fear conditioning data
Voff = .2;
ks = 5;
MaxFiber = 90;


switch sim
    case '4A'
        
        %Stimuli used in the fear conditioning Design. Each row of the arrays
        %represents a single time step. Each column represents whether an
        %element was present for the particular stimulus in a binary manner.
        %randTrials adds the ITI between trials. 
        %Column 1 : Auditory Cue
        %Column 2 : Context
        D.TB =  [0 1; 1 0; 0 1; 0 1];
        
        %This sets the initial alphas
        alpha = [0.3 0.3];
        
        R = struct();
            
        %stimuli are randomized within each execution of the TDRL model.
        %This is included for designs in which the order of stimulus
        %presentations affects the value.
        for i = 1:100
            
            %this randomizes trial order for each stage of conditioning
            [S{1},Log(1)] = RandTrials(3,D,'TB');
            [S{2},Log(2)] = RandTrials(6,D,'TB');
            
            %This gets a vector of the stage numbers. This is used incase
            %there are identical stimuli in each stage.
            stg = stagenum(S);
            
            M = cell2mat(S');
            
            %L keeps track of the logicals for the first state of each
            %stimulus
            [Lmerge] = merge_stage_logicals(Log);
            %logicals of the first state of each stimulus
            %These are the states that get used for plotting
            TBSTART = Lmerge.TB & stg == 1;
            L.tb = nstate(TBSTART,2) & stg == 1;
            L.b_1 = nstate(TBSTART,1) & stg == 1;
            TBSTART_e = Lmerge.TB & stg == 2;
            L.tb_e = nstate(TBSTART_e,2) & stg == 2;
            L.b_e = nstate(TBSTART_e,1) & stg == 2;
            FN = fieldnames(L);
            
            r = zeros(size(L.tb));
            %rewarded states
            r(nstate(L.tb,2)) = -1;
           
            %This implements tonic inhibition during just A in AX
            %trials
             opto = -0.7*(L.tb);
            
            %Initial values of each of the stimuli. 
            ui = [0; 0];
            
            results = linearTDRL_Lambda_LH_Inact_Mac(M,r,'Vinitial',ui,'alpha',alpha,'opto',opto,'lambda',0.7,'Mac_Trials',size(D.TB,1),'theta',0.2);
            for k = 1:numel(FN) 
                R(1).(FN{k})(i,:) = CR(results.V(L.(FN{k})),MaxFiber,Voff,ks); 
            end
            
        end
           
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        %The following code creates a panel with line plots of each stage of conditioning
        %as a different column, as well as a final column with a bar plot
        %of the behavior during the extinction test
        
        
        %This determines the subplot locations for the panel
        plotwidths = [0 .15 .1 .1 .1];
        plotheights = [.15 .15 .15];
        spcx = [0.1 0.01 0.01 0.1];
        spcy = [0 0.15 0.15];
        positions = cell(3,4);
        
        for sx = 1:4
            for sy = 1:3
                
                positions{4-sy,sx} = [sum(spcx(1:sx)) + sum(plotwidths(1:sx)), sum(spcy(1:sy)) + sum(plotheights(1:sy)), plotwidths(sx + 1),plotheights(sy)];
            end
        end
        
        V = zeros(numel(R),2);
        
        for k = 1:numel(R)
            V(k,:) = mean([R(k).tb_e(:,1) R(k).b_e(:,1)]);
        end
        
        figure
        hold on
        models = {''};
        m = 1;
            
            subplot('Position',positions{m,4})
            b = bar(V(m,:));
            b.FaceColor = 'flat';
            b.CData(2,:) = [1 0 0];
    
            set(gca,'FontSize',11,'XLim', [0 3],'YLim',[0 100],'XTickLabel', {'Tone' 'CXT'}, 'YTick', 0:20:100);
            title ('Test','FontSize',15)
            xlabel(models{m},'FontSize',15)
            ylabel('Freezing','FontSize',15);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,2})
            title(models{m})
            Resp = [];
            for i = 1:3
                range = i;
                Resp(i,:) = mean([mean(R(m).tb(:,range),2) mean(R(m).b_1(:,range),2)]);
            end
            plot(1:3,Resp(:,1),'-ob')
            hold on
            plot(1:3,Resp(:,2),'-or')
            set(gca,'FontSize',11, 'XTick', 1:3, 'XLim',[.5 3.5],'YLim',[0 100],'YTick',0:20:100);%0:10:50)
            title ('Cond','FontSize',15)
            xlabel('Trials','FontSize',15)
            ylabel('Freezing','FontSize',15);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,3})
            title(models{m})
            Resp = [];
            for i = 1:6
                range = i;%(i-1)*3+1:i*3;
                Resp(i,:) = mean([mean(R(m).tb_e(:,range),2) mean(R(m).b_e(:,range),2)]);
            end
            plot(1:6,Resp(:,1),'-ob')
            hold on
            plot (1:6,Resp(:,2),'-or')
            
            set(gca,'FontSize',11, 'XTick', 1:6, 'XLim',[0.5 6.5],'YLim',[0 100],'YTick',[]);
            %  ylabel('CR','FontSize',15);
           title ('Test','FontSize',15)
            xlabel('Trials','FontSize',15);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        
case '4B - Low Alpha'
        
        %Stimuli used in the fear conditioning Design. Each row of the arrays
        %represents a single time step. Each column represents whether an
        %element was present for the particular stimulus in a binary manner.
        %randTrials adds the ITI between trials. 
        %Column 1 : Auditory Cue
        %Column 2 : Context B
        D.TB =  [0 1; 1 0; 0 1; 0 1];
        
        %This sets the initial alphas
        alpha = [0.3 0.01];
        
        R = struct();
            
        %stimuli are randomized within each execution of the TDRL model.
        %This is included for designs in which the order of stimulus
        %presentations affects the value.
        for i = 1:100
            
            %this randomizes trial order for each stage of conditioning
            [S{1},Log(1)] = RandTrials(3,D,'TB');
            [S{2},Log(2)] = RandTrials(6,D,'TB');
            
            %This gets a vector of the stage numbers. This is used incase
            %there are identical stimuli in each stage.
            stg = stagenum(S);
            
            M = cell2mat(S');
            
            %L keeps track of the logicals for the first state of each
            %stimulus
            [Lmerge] = merge_stage_logicals(Log);
            %logicals of the first state of each stimulus
            %These are the states that get used for plotting
            TBSTART = Lmerge.TB & stg == 1;
            L.tb = nstate(TBSTART,2) & stg == 1;
            L.b_1 = nstate(TBSTART,1) & stg == 1;
            TBSTART_e = Lmerge.TB & stg == 2;
            L.tb_e = nstate(TBSTART_e,2) & stg == 2;
            L.b_e = nstate(TBSTART_e,1) & stg == 2;
            FN = fieldnames(L);
            
            r = zeros(size(L.tb));
            %rewarded states
            r(nstate(L.tb,2)) = -1;
           
            %This implements tonic inhibition during just A in AX
            %trials
             opto = -0.7*(L.tb);
            
            %Initial values of each of the stimuli. 
            ui = [0; 0];
            
            %%%CONTROL%%%%
            results = linearTDRL_Lambda_LH_Inact_Mac(M,r,'Vinitial',ui,'alpha',alpha,'opto',opto,'lambda',0.7,'Mac_Trials',size(D.TB,1),'theta',0.2);
            for k = 1:numel(FN) 
                R(1).(FN{k})(i,:) = CR(results.V(L.(FN{k})),MaxFiber,Voff,ks); 
            end

            %%%%EXPERIMENTAL%%%%%%%
            results = linearTDRL_Lambda_LH_Inact_Mac(M,r,'Vinitial',ui,'alpha',alpha,'opto',opto,'opto_update','update','lambda',0.7,'Mac_Trials',size(D.TB,1),'theta',0.2);
            for k = 1:numel(FN)
                    R(2).(FN{k})(i,:) = CR(results.V(L.(FN{k})),MaxFiber,Voff,ks); 
     
            end
            
        end
           
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        %The following code creates a panel with line plots of each stage of conditioning
        %as a different column, as well as a final column with a bar plot
        %of the behavior during the extinction test
        
        
        %This determines the subplot locations for the panel
        plotwidths = [0 .15 .1 .1 .1];
        plotheights = [.15 .15 .15];
        spcx = [0.1 0.01 0.01 0.1];
        spcy = [0 0.15 0.15];
        positions = cell(3,4);
        
        for sx = 1:4
            for sy = 1:3
                
                positions{4-sy,sx} = [sum(spcx(1:sx)) + sum(plotwidths(1:sx)), sum(spcy(1:sy)) + sum(plotheights(1:sy)), plotwidths(sx + 1),plotheights(sy)];
            end
        end
        
        V = zeros(numel(R),2);
        
        for k = 1:numel(R)
            V(k,:) = mean([R(k).tb_e(:,1) R(k).b_e(:,1)]);
        end
        
        figure
        hold on
        models = {'eYFP', 'NpHR'};
        for m = 1:2
            
            subplot('Position',positions{m,4})
            b = bar(V(m,:));
            b.FaceColor = 'flat';
            b.CData(2,:) = [1 0 0];
    
            set(gca,'FontSize',11,'XLim', [0 3],'YLim',[0 100],'XTickLabel', {'Tone' 'CXT'}, 'YTick', 0:20:100);
            title ('Test','FontSize',15)
            xlabel(models{m},'FontSize',15)
            ylabel('Freezing','FontSize',15);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,2})
            title(models{m})
            Resp = [];
            for i = 1:3
                range = i;
                Resp(i,:) = mean([mean(R(m).tb(:,range),2) mean(R(m).b_1(:,range),2)]);
            end
            plot(1:3,Resp(:,1),'-ob')
            hold on
            plot(1:3,Resp(:,2),'-or')
            set(gca,'FontSize',11, 'XTick', 1:3, 'XLim',[.5 3.5],'YLim',[0 100],'YTick',0:20:100);%0:10:50)
            title ('Cond','FontSize',15)
            xlabel('Trials','FontSize',15)
            ylabel('Freezing','FontSize',15);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,3})
            title(models{m})
            Resp = [];
            for i = 1:6
                range = i;%(i-1)*3+1:i*3;
                Resp(i,:) = mean([mean(R(m).tb_e(:,range),2) mean(R(m).b_e(:,range),2)]);
            end
            plot(1:6,Resp(:,1),'-ob')
            hold on
            plot (1:6,Resp(:,2),'-or')
            
            set(gca,'FontSize',11, 'XTick', 1:6, 'XLim',[0.5 6.5],'YLim',[0 100],'YTick',[]);
            %  ylabel('CR','FontSize',15);
           title ('Test','FontSize',15)
            xlabel('Trials','FontSize',15);
        end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

case '4B - App Value'
        
        %Stimuli used in the fear conditioning Design. Each row of the arrays
        %represents a single time step. Each column represents whether an
        %element was present for the particular stimulus in a binary manner.
        %randTrials adds the ITI between trials. 
        %Column 1 : Auditory Cue
        %Column 2 : Context B
        D.TB =  [0 1; 1 0; 0 1; 0 1];
        
        %This sets the initial alphas
        alpha = [0.3 0.3];
        
        R = struct();
            
        %stimuli are randomized within each execution of the TDRL model.
        %This is included for designs in which the order of stimulus
        %presentations affects the value.
        for i = 1:100
            
            %this randomizes trial order for each stage of conditioning
            [S{1},Log(1)] = RandTrials(3,D,'TB');
            [S{2},Log(2)] = RandTrials(6,D,'TB');
            
            %This gets a vector of the stage numbers. This is used incase
            %there are identical stimuli in each stage.
            stg = stagenum(S);
            
            M = cell2mat(S');
            
            %L keeps track of the logicals for the first state of each
            %stimulus
            [Lmerge] = merge_stage_logicals(Log);
            %logicals of the first state of each stimulus
            %These are the states that get used for plotting
            TBSTART = Lmerge.TB & stg == 1;
            L.tb = nstate(TBSTART,2) & stg == 1;
            L.b_1 = nstate(TBSTART,1) & stg == 1;
            TBSTART_e = Lmerge.TB & stg == 2;
            L.tb_e = nstate(TBSTART_e,2) & stg == 2;
            L.b_e = nstate(TBSTART_e,1) & stg == 2;
            FN = fieldnames(L);
            
            r = zeros(size(L.tb));
            %rewarded states
            r(nstate(L.tb,2)) = -1;
           
            %This implements tonic inhibition during just A in AX
            %trials
             opto = -0.7*(L.tb);
            
            %Initial values of each of the stimuli. Here the initial value
            %of the context is increased due to possible learning during
            %appetitive conditioning
            ui = [0; 0.3];
            
            %%%CONTROL%%%%
            results = linearTDRL_Lambda_LH_Inact_Mac(M,r,'Vinitial',ui,'alpha',alpha,'opto',opto,'lambda',0.7,'Mac_Trials',size(D.TB,1),'theta',0.2);
            for k = 1:numel(FN) 
                R(1).(FN{k})(i,:) = CR(results.V(L.(FN{k})),MaxFiber,Voff,ks); 
            end

            %%%%EXPERIMENTAL%%%%%%%
            results = linearTDRL_Lambda_LH_Inact_Mac(M,r,'Vinitial',ui,'alpha',alpha,'opto',opto,'opto_update','update','lambda',0.7,'Mac_Trials',size(D.TB,1),'theta',0.2);
            for k = 1:numel(FN)
                    R(2).(FN{k})(i,:) = CR(results.V(L.(FN{k})),MaxFiber,Voff,ks); 
     
            end
            
        end
           
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        %The following code creates a panel with line plots of each stage of conditioning
        %as a different column, as well as a final column with a bar plot
        %of the behavior during the extinction test
        
        
        %This determines the subplot locations for the panel
        plotwidths = [0 .15 .1 .1 .1];
        plotheights = [.15 .15 .15];
        spcx = [0.1 0.01 0.01 0.1];
        spcy = [0 0.15 0.15];
        positions = cell(3,4);
        
        for sx = 1:4
            for sy = 1:3
                
                positions{4-sy,sx} = [sum(spcx(1:sx)) + sum(plotwidths(1:sx)), sum(spcy(1:sy)) + sum(plotheights(1:sy)), plotwidths(sx + 1),plotheights(sy)];
            end
        end
        
        V = zeros(numel(R),2);
        
        for k = 1:numel(R)
            V(k,:) = mean([R(k).tb_e(:,1) R(k).b_e(:,1)]);
        end
        
        figure
        hold on
        models = {'eYFP', 'NpHR'};
        for m = 1:2
            
            subplot('Position',positions{m,4})
            b = bar(V(m,:));
            b.FaceColor = 'flat';
            b.CData(2,:) = [1 0 0];
    
            set(gca,'FontSize',11,'XLim', [0 3],'YLim',[0 100],'XTickLabel', {'Tone' 'CXT'}, 'YTick', 0:20:100);
            title ('Test','FontSize',15)
            xlabel(models{m},'FontSize',15)
            ylabel('Freezing','FontSize',15);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,2})
            title(models{m})
            Resp = [];
            for i = 1:3
                range = i;
                Resp(i,:) = mean([mean(R(m).tb(:,range),2) mean(R(m).b_1(:,range),2)]);
            end
            plot(1:3,Resp(:,1),'-ob')
            hold on
            plot(1:3,Resp(:,2),'-or')
            set(gca,'FontSize',11, 'XTick', 1:3, 'XLim',[.5 3.5],'YLim',[0 100],'YTick',0:20:100);%0:10:50)
            title ('Cond','FontSize',15)
            xlabel('Trials','FontSize',15)
            ylabel('Freezing','FontSize',15);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,3})
            title(models{m})
            Resp = [];
            for i = 1:6
                range = i;%(i-1)*3+1:i*3;
                Resp(i,:) = mean([mean(R(m).tb_e(:,range),2) mean(R(m).b_e(:,range),2)]);
            end
            plot(1:6,Resp(:,1),'-ob')
            hold on
            plot (1:6,Resp(:,2),'-or')
            
            set(gca,'FontSize',11, 'XTick', 1:6, 'XLim',[0.5 6.5],'YLim',[0 100],'YTick',[]);
            %  ylabel('CR','FontSize',15);
           title ('Test','FontSize',15)
            xlabel('Trials','FontSize',15);
        end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
 

case '4C'
        
        %Stimuli used in the fear conditioning Design. Each row of the arrays
        %represents a single time step. Each column represents whether an
        %element was present for the particular stimulus in a binary manner.
        %randTrials adds the ITI between trials. 
        %Column 1 : Auditory Cue
        %Column 2 : Context B
        D.TB =  [0 1; 1 0; 0 1; 0 1];
          
        %This sets the initial alphas
        alpha = [0.3 0.3];
        
        R = struct();
            
        %stimuli are randomized within each execution of the TDRL model.
        %This is included for designs in which the order of stimulus
        %presentations affects the value.
        for i = 1:100
            
            %this randomizes trial order for each stage of conditioning
            [S{1},Log(1)] = RandTrials(3,D,'TB');
            [S{2},Log(2)] = RandTrials(6,D,'TB');
            
            %This gets a vector of the stage numbers. This is used incase
            %there are identical stimuli in each stage.
            stg = stagenum(S);
            
            M = cell2mat(S');
            
            %L keeps track of the logicals for the first state of each
            %stimulus
            [Lmerge] = merge_stage_logicals(Log);
            %logicals of the first state of each stimulus
            %These are the states that get used for plotting
            TBSTART = Lmerge.TB & stg == 1;
            L.tb = nstate(TBSTART,2) & stg == 1;
            L.b_1 = nstate(TBSTART,1) & stg == 1;
            TBSTART_e = Lmerge.TB & stg == 2;
            L.tb_e = nstate(TBSTART_e,2) & stg == 2;
            L.b_e = nstate(TBSTART_e,1) & stg == 2;
            FN = fieldnames(L);
            
            r = zeros(size(L.tb));
            %rewarded states
            r(nstate(L.tb,2)) = -1;
           
            %This implements tonic inhibition during just A in AX
            %trials
             opto = -0.7*(L.tb);
            
            %Initial values of each of the stimuli. 
            ui = [0; 0];
            
            %%%CONTROL%%%%
            results = linearTDRL_Lambda_LH_Inact_Mac(M,r,'Vinitial',ui,'alpha',alpha,'opto',opto,'lambda',0.7,'Mac_Trials',size(D.TB,1),'theta',0.2);
            for k = 1:numel(FN) 
                R(1).(FN{k})(i,:) = CR(results.V(L.(FN{k})),MaxFiber,Voff,ks); 
            end

            %%%%EXPERIMENTAL%%%%%%%
            results = linearTDRL_Lambda_LH_Inact_Mac(M,r,'Vinitial',ui,'alpha',alpha,'opto',opto,'opto_update','update','lambda',0.7,'Mac_Trials',size(D.TB,1),'theta',0.2);
            for k = 1:numel(FN)
                    R(2).(FN{k})(i,:) = CR(results.V(L.(FN{k})),MaxFiber,Voff,ks); 
     
            end
            
        end
           
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        %The following code creates a panel with line plots of each stage of conditioning
        %as a different column, as well as a final column with a bar plot
        %of the behavior during the extinction test
        
        
        %This determines the subplot locations for the panel
        plotwidths = [0 .15 .1 .1 .1];
        plotheights = [.15 .15 .15];
        spcx = [0.1 0.01 0.01 0.1];
        spcy = [0 0.15 0.15];
        positions = cell(3,4);
        
        for sx = 1:4
            for sy = 1:3
                
                positions{4-sy,sx} = [sum(spcx(1:sx)) + sum(plotwidths(1:sx)), sum(spcy(1:sy)) + sum(plotheights(1:sy)), plotwidths(sx + 1),plotheights(sy)];
            end
        end
        
        V = zeros(numel(R),2);
        
        for k = 1:numel(R)
            V(k,:) = mean([R(k).tb_e(:,1) R(k).b_e(:,1)]);
        end
        
        figure
        hold on
        models = {'eYFP', 'NpHR'};
        for m = 1:2
            
            subplot('Position',positions{m,4})
            b = bar(V(m,:));
            b.FaceColor = 'flat';
            b.CData(2,:) = [1 0 0];
    
            set(gca,'FontSize',11,'XLim', [0 3],'YLim',[0 100],'XTickLabel', {'Tone' 'CXT'}, 'YTick', 0:20:100);
            title ('Test','FontSize',15)
            xlabel(models{m},'FontSize',15)
            ylabel('Freezing','FontSize',15);
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,2})
            title(models{m})
            Resp = [];
            for i = 1:3
                range = i;
                Resp(i,:) = mean([mean(R(m).tb(:,range),2) mean(R(m).b_1(:,range),2)]);
            end
            plot(1:3,Resp(:,1),'-ob')
            hold on
            plot(1:3,Resp(:,2),'-or')
            set(gca,'FontSize',11, 'XTick', 1:3, 'XLim',[.5 3.5],'YLim',[0 100],'YTick',0:20:100);%0:10:50)
            title ('Cond','FontSize',15)
            xlabel('Trials','FontSize',15)
            ylabel('Freezing','FontSize',15);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot('Position',positions{m,3})
            title(models{m})
            Resp = [];
            for i = 1:6
                range = i;%(i-1)*3+1:i*3;
                Resp(i,:) = mean([mean(R(m).tb_e(:,range),2) mean(R(m).b_e(:,range),2)]); %#ok<*AGROW>
            end
            plot(1:6,Resp(:,1),'-ob')
            hold on
            plot (1:6,Resp(:,2),'-or')
            
            set(gca,'FontSize',11, 'XTick', 1:6, 'XLim',[0.5 6.5],'YLim',[0 100],'YTick',[]);
            %  ylabel('CR','FontSize',15);
           title ('Test','FontSize',15)
            xlabel('Trials','FontSize',15);
        end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
 

end
end

function S1 = fstate(i,S) %#ok<*DEFNU>
%Provides a logical vector (N x 1) of the first state of 
%stimulus i by template matching i within the state matrix J. This script
%will check up to three rows of stimulus i within the state matrix J.
%Several rows of each stimulus are checked incase the first row or two of
%stimulus i are identical to another stimulus.
    S1 = false(size(S,1),1);
    
    if size(i,1) == 1
        S1 = all(S == i(1,:),2);
    elseif size(i,1) == 2
        S1 = all(S == i(1,:),2) & (all([S(2:end,:); zeros(1,size(i,2))] == i(2,:),2));
    elseif size(i,1) > 2  
        S1 = all(S == i(1,:),2) & (all([S(2:end,:); zeros(1,size(i,2))] == i(2,:),2)) & (all([S(3:end,:); zeros(2,size(i,2))] == i(3,:),2));
    end    
end

function Sn = nstate(X,n)
%Provides a logical vector of state n within a given stimulus by shifting
%the logicals of stimulus X by n values. For example, n = 2 provides an N x 1
%logical vector of the second row of the stimulus in X, where X is an N x 1
%logical vector. If n = -1, this returns the logicals for the prior state.

    if n >= 0
        Sn = [false(n-1,1); X(1:end-n+1)];
    else
        Sn = [X(-1*n+1:end); false(-1*n,1)];
    end
end

function Stage = stagenum(S)
%provides an N X 1 vector of the stage number of conditioning for
%the state matrix J. J is a cell array whose size is 1 X the number of
%stages.
    Stage = zeros(size(cell2mat(S'),1),1);
    
    ind_start = 1;
    for i = 1:numel(S)
        ind_end = ind_start + size(S{i},1) - 1;
        
        Stage(ind_start:ind_end) = i;
        ind_start = ind_end + 1;
        
    end    

end

function [out, Logicals] = RandTrials(n,D,varargin)
% This function randomizes presentations of stimuli for the linearTDSR function

%INPUTS:
%n is the number of presentations for each stimulus sequence
%D is the design matrix, each field represents a different stimulus

%the varargin are the unique stimulus sequences to be included 
    %each stimulus is an [s X D] matrix. (s = # of states within a single
    %stimulus sequence, D = # of features 

%OUTPUT: [n*s*(nargin - 1) X D] of randomized stimulus presentations
    
%Note that zeros(1 X D) are added between each presentation

out = [];
ncue = nargin - 2;
curr = 0;
Stimuli = fieldnames(D);

%Check whether the varargin strings match a fieldname of the structure
for v = 1:ncue
    if ~any(strcmp(varargin{v}, Stimuli))
        disp('one of the stimuli in the RandTrials function does not match the possible stimuli available in the design structure');
        return
    end
end    

%Get the size of the stimuli
stim_length = size(D.(varargin{1}),1);

%The Logicals keeps track of the first step of each stimulus
for s = 1:numel(Stimuli)
    Logicals.(Stimuli{s}) = false(n*ncue*stim_length,1);
end


for i = 1:n
   
    A = randperm(ncue,ncue);
    
    for j = 1:ncue
        
        Curr_Stim = D.(varargin{A(j)});
        
        str = size(Curr_Stim,1);
        
        %This sets the start of stimulus to true
        Logicals.(varargin{A(j)})(curr+1) = true;
        
        for k = 1:str
            curr = curr + 1;
            out(curr, :) = Curr_Stim(k,:);
        end
%         curr = curr + 1;
%         out(curr, :) = zeros(1, size(varargin{A(j)},2));
    end
end    
end    

function [L_out] = merge_stage_logicals(L)

fn = fieldnames(L);
for i = 1:numel(fn)
    L_out.(fn{i}) = [];
    for s = 1:numel(L)
        L_out.(fn{i}) = [L_out.(fn{i}); L(s).(fn{i})];
    end
end
end
