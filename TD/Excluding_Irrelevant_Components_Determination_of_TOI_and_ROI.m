% This code is to apply 2-way ANOVA on the decomposed subjects' signatures after
% applying tensor-decomposition.
% The significant bases will be saved to ANOVA/significant_entity,
% significant_action, significant_interaction
% in this version, ANOVA function was used for unbalanced conditions.

% This verion has seperate plots for each condition (i.e., human, function; human, non-fuctional;
% claw, functional; claw, nonfunctional) of only the significant/selected
% components. 

% This version also combines multiple temporal component together (by
% temporal times 'core') and doesnt times the spatial (only adds the spatial components together)

% Also, normalizes the data for spatial components (after multiplying 8by10 core
% so that it is more interpertable. Normalization is just for
% visualization not used prior to plotting

% Also excludes irrlevant components when determining the TOI and ROI


%% Data directory and hyperparameters
clc; clear all;
td_dir='./Results/lmlra_decomp'; %the directory of the decompostion after running Tensor_Construction_Tensor_Decomposition.m
load(strcat(td_dir,'/LMLRA_all')); %loading variables
noSamples_time=size(U_left{1},1); %number of time points
temp_r_l=size(U_left{1},2); %number of temporal components for left hemisphere
temp_r_r=size(U_right{1},2); %number of temporal components for right hemisphere
spac_r_l=size(U_left{2},2); %number of spatial components for left hemisphere
spac_r_r=size(U_right{2},2); %number of spatial components for right hemisphere
subj_size=size(U_left{3},2); % number of subjects 
fs=50; %sampling frequency 
conditions_title={'Human hand, function', 'Human hand, non-function', 'Mechanical hand, function', ...
    'Mechanical hand, non-function'};


%parameters for figures 
xbase=1;ybase=1; activations=zeros(4+xbase,8+ybase);
Infimg_left=imread('left pic.gif',1);
Infimg_left=imresize(Infimg_left,0.2);%(Infimg,0.09);
Infimg_right=imread('right pic.gif',1);
Infimg_right=imresize(Infimg_right,0.5);%(Infimg,0.09);
    

%Calculating temporal and spatial components for all subjects 
%core times subjects for left hemisphere  matrix time second matrix
subjectWeightsL(1:temp_r_l,1:spac_r_l,1:subj_size,1:subj_size)=0;
for t = 1:temp_r_l %number of temporal components
   for c = 1:spac_r_l %number of spatial components
      for s = 1:subj_size %number of subjects component
          for i = 1:subj_size %subject number
            subjectWeightsL(t,c,s,i) = Sol_left(t,c,s)* U_left{1,3}(i,s);
          end
      end
   end
end
%core times subjects for right hemisphere
subjectWeightsR(1:temp_r_r,1:spac_r_r,1:subj_size,1:subj_size)=0;
for t = 1:temp_r_r %number of temporal components
   for c = 1:spac_r_r %number of spatial components
      for s = 1:subj_size %number of subjects component
          for i = 1:subj_size %subject number
            subjectWeightsR(t,c,s,i) = Sol_right(t,c,s)* U_right{1,3}(i,s);
          end
      end
   end
end


%summing the components on subject number for left hemisphere
SumSubjectWeightsL(1:temp_r_l,1:spac_r_l,1:subj_size)=0;
for t = 1:temp_r_l %number of temporal components
   for c = 1:spac_r_l %number of spatial components
      for i = 1:subj_size %subjects number
          SumSubjectWeightsL(t,c,i) = sum(subjectWeightsL(t,c,:,i));
      end
   end
end

%summing the components on subject number for right hemisphere
SumSubjectWeightsR(1:temp_r_r,1:spac_r_r,1:subj_size)=0;
for t = 1:temp_r_r %number of temporal components
   for c = 1:spac_r_r %number of spatial components
      for i = 1:subj_size %subjects number
          SumSubjectWeightsR(t,c,i) = sum(subjectWeightsR(t,c,:,i));
      end
   end
end


%Finding ANOVA for the signatures 
r_components_l=size_core_left(1)*size_core_left(2);  %r_components_l=size(U_left{1},2); 
r_components_r=size_core_right(1)*size_core_right(2); %r_components_r=size(U_right{1},2);
F_entity=[]; F_action=[]; F_interaction=[];
T_anova_left = table('Size',[r_components_l 6],'VariableNames',{'F_entity', ...
    'P_entity','F_action','P_action','F_interaction','P_interaction'}, ...
    'VariableTypes',{'double', 'double', 'double', 'double', 'double', 'double'});
T_anova_right = table('Size',[r_components_r 6],'VariableNames',{'F_entity', ...
    'P_entity','F_action','P_action','F_interaction','P_interaction'}, ...
    'VariableTypes',{'double', 'double', 'double', 'double', 'double', 'double'});

%make two vectors for each factor of anova
entity_factor=condition_ind;
entity_factor(or(condition_ind==1,condition_ind==2))=1;
entity_factor(or(condition_ind==3,condition_ind==4))=2;
action_factor=condition_ind;
action_factor(or(condition_ind==1,condition_ind==3))=1;
action_factor(or(condition_ind==2,condition_ind==4))=2;

%anova for significant temporal and spatial commponent pairing  
 %Left hemisphere:
 index =1;
 index2cntL = zeros(size_core_left(1)*size_core_left(2),2);
for t=1:size_core_left(1) %number of temporal components %10 %spatial components 
    for c = 1:size_core_left(2) %number of spatial components
        [p,F_table]=anovan(squeeze(SumSubjectWeightsL(t,c,:)),{entity_factor,action_factor},...
            'model','interaction','varnames',{'entity','action'},'display','off');
        T_anova_left{index,1}=F_table{2,6}; T_anova_left{index,2}=p(1);
        T_anova_left{index,3}=F_table{3,6}; T_anova_left{index,4}=p(2);
        T_anova_left{index,5}=F_table{4,6}; T_anova_left{index,6}=p(3);
        index2cntL(index,:) = [t,c];  
        index = index+1;
    end 
end
%Right hemisphere:
index =1;
index2cntR = zeros(size_core_right(1)*size_core_right(2),2);
for t=1:size_core_right(1)%number of temporal components
    for c = 1:size_core_right(2) %number of spatial components
        [p,F_table]=anovan(squeeze(SumSubjectWeightsR(t,c,:)),{entity_factor,action_factor},...
            'model','interaction','varnames',{'entity','action'},'display','off');
        T_anova_right{index,1}=F_table{2,6}; T_anova_right{index,2}=p(1);
        T_anova_right{index,3}=F_table{3,6}; T_anova_right{index,4}=p(2);
        T_anova_right{index,5}=F_table{4,6}; T_anova_right{index,6}=p(3);
        index2cntR(index,:) = [t,c];
        index = index+1;
    end
end

% Copying the significant patterns or basis to ANOVA/significant_entity,
% significant_action, significant_interaction
entity_dir=strcat(td_dir,'/anova','/significant_entity/');
action_dir=strcat(td_dir,'/anova','/significant_action/');
inter_dir=strcat(td_dir,'/anova','/significant_interaction/');
if ~exist(strcat(td_dir,'/anova'))
    mkdir(strcat(td_dir,'/anova'))
    mkdir(entity_dir)
    mkdir(action_dir)
    mkdir(inter_dir)
end
save(strcat(td_dir,'/anova/anova_results'),'T_anova_left','T_anova_right');

%plotting significant components
for r=1:r_components_l 
    %Left
    if T_anova_left{r,2}<0.05 %find significant entity
        FigH = figure('Position', get(0, 'Screensize'),'visible',false);
        subplot(1,3,1); %plotting temporal component
            plot(U_left{1}(:,index2cntL(r,1)),'LineWidth',2)
            title(strcat('Temporal Component ', num2str(index2cntL(r,1))))
            xlabel('Time (sec)','Fontsize',18,'FontName','Times New Roman')
            ylabel('Magnitude','Fontsize',18,'FontName','Times New Roman')
            xticks(0:fs*2:noSamples_time)
            xticklabels(0:2:noSamples_time/fs)

            set(gca,'box','off', 'FontSize', 16);
       
        subplot(1,3,2); % plotting spatial component
            data=U_left{2}(:,index2cntL(r,2));
    
            activations(xbase+3,ybase+[2,4,6])=data(1:3); %the channels will be flipped when using imshow
            activations(xbase+2,ybase+[1,3,5,7])=data(4:7);
            activations(xbase+1,ybase+[2,4,6])=data(8:10);
            
            %ploting the emitters and detectors placed on left hemisphere
            plot(xbase+[1,5,3,7],ybase+[1,1,3,3],'ks','MarkerFaceColor','k','MarkerSize', 10)
            hold on
            plot(xbase+[3,7,1,5],ybase+[1,1,3,3],'ro','MarkerFaceColor','r','MarkerSize', 10)

            x = [1:8+xbase];
            y = [1:4+ybase];  
            %selecting the contour levels between the min and max
            %activations 
            levels=linspace(min(min(activations)),max(max(activations)),10); 
            contour(x,y,activations,levels,'ShowText','off','LineWidth',1.5) %narrow lines means negative levels. 
            contourcmap('hsv', levels,'Colorbar', 'on', ...
           'Location', 'horizontal', ...
           'TitleString', strcat('The 10 levels of magnitude are between ',num2str(min(min(data)))...
                ,' to ',num2str(max(max(data)))));
            axis equal                                % set axis units to be the same size
            axis square off
            imgh=imshow(repmat(Infimg_left,1,1,3), 'XData', [-11 14], 'YData', [-9 15]);
            imgh.AlphaData = .2;
            hold off
            xlim([0.77 9.716])
            ylim([0.8543 7.795])
            title(strcat('Spatial Component ', num2str(index2cntL(r,2))));
        subplot(1,3,3); % plotting subject 
             hold on;
            starti=1;
             subjectVal(1:subj_size)=0;
            for conditionI=1:length(conditions)
                endi=starti+length(find(condition_ind==conditionI))-1;
                subjectVal(starti:endi,1)=SumSubjectWeightsL(index2cntL(r,1),index2cntL(r,2),starti:endi);
                bar(conditionI, mean(subjectVal(starti:endi)));
                starti=endi+1;
            end
            xlabel('condition #','Fontsize',18,'FontName','Times New Roman')
            ylabel('Signature','Fontsize',18,'FontName','Times New Roman')
            set(gca,'box','off', 'FontSize', 16);
            legend(conditions_title,'Location','best');
        F    = getframe(FigH);
        imwrite(F.cdata, fullfile(entity_dir,...
                strcat('left mean',num2str(r),'.png')), 'png')
        close    
        
    end
    
    if T_anova_left{r,4}<0.05 %significant action sequence
        FigH = figure('Position', get(0, 'Screensize'),'visible',false);
        subplot(1,3,1); %plotting temporal 
            plot(U_left{1}(:,index2cntL(r,1)),'LineWidth',2)
            title(strcat('Temporal Component ', num2str(index2cntL(r,1))))
            xlabel('Time (sec)','Fontsize',18,'FontName','Times New Roman')
            ylabel('Magnitude','Fontsize',18,'FontName','Times New Roman')
            xticks(0:fs*2:noSamples_time)
            xticklabels(0:2:noSamples_time/fs)
            set(gca,'box','off', 'FontSize', 16);
       
        subplot(1,3,2); %plotting spatial 
            data=U_left{2}(:,index2cntL(r,2));
    
            activations(xbase+3,ybase+[2,4,6])=data(1:3); %the channels will be flipped when using imshow
            activations(xbase+2,ybase+[1,3,5,7])=data(4:7);
            activations(xbase+1,ybase+[2,4,6])=data(8:10);
            
            %ploting the emitters and detectors placed on left hemisphere
            plot(xbase+[1,5,3,7],ybase+[1,1,3,3],'ks','MarkerFaceColor','k','MarkerSize', 10)
            hold on
            plot(xbase+[3,7,1,5],ybase+[1,1,3,3],'ro','MarkerFaceColor','r','MarkerSize', 10)

            x = [1:8+xbase];
            y = [1:4+ybase];  
            %selecting the contour levels between the min and max
            %activations 
            levels=linspace(min(min(activations)),max(max(activations)),10);
            contour(x,y,activations,levels,'ShowText','off','LineWidth',1.5) %narrow lines means negative levels. 
            contourcmap('hsv', levels,'Colorbar', 'on', ...
           'Location', 'horizontal', ...
           'TitleString', strcat('The 10 levels of magnitude are between ',num2str(min(min(data)))...
                ,' to ',num2str(max(max(data)))));
            axis equal                                % set axis units to be the same size
            axis square off
            imgh=imshow(repmat(Infimg_left,1,1,3), 'XData', [-11 14], 'YData', [-9 15]);
            imgh.AlphaData = .2;
            hold off
            xlim([0.77 9.716])
            ylim([0.8543 7.795])
            title(strcat('Spatial Component ', num2str(index2cntL(r,2))));

        subplot(1,3,3);%plotting subject 
             hold on;
            starti=1;
             subjectVal(1:subj_size)=0;
            for conditionI=1:length(conditions)
                endi=starti+length(find(condition_ind==conditionI))-1;
                subjectVal(starti:endi,1)=SumSubjectWeightsL(index2cntL(r,1),index2cntL(r,2),starti:endi);
                bar(conditionI, mean(subjectVal(starti:endi)));
                starti=endi+1;
            end
            xlabel('condition #','Fontsize',18,'FontName','Times New Roman')
            ylabel('Signature','Fontsize',18,'FontName','Times New Roman')
            set(gca,'box','off', 'FontSize', 16);
            legend(conditions_title,'Location','best');

        F    = getframe(FigH);

        imwrite(F.cdata, fullfile(action_dir,...
                strcat('left mean',num2str(r),'.png')), 'png')
        close 

    end
    
    if T_anova_left{r,6}<0.05 % significant interaction 
         FigH = figure('Position', get(0, 'Screensize'),'visible',false);
        subplot(1,3,1); % plotting time 
            plot(U_left{1}(:,index2cntL(r,1)),'LineWidth',2)
            title(strcat('Temporal Component ', num2str(index2cntL(r,1))))
            xlabel('Time (sec)','Fontsize',18,'FontName','Times New Roman')
            ylabel('Magnitude','Fontsize',18,'FontName','Times New Roman')
            xticks(0:fs*2:noSamples_time)
            xticklabels(0:2:noSamples_time/fs)
            set(gca,'box','off', 'FontSize', 16);
       
        subplot(1,3,2); %plotting spatial 
     
            data=U_left{2}(:,index2cntL(r,2));
    
            activations(xbase+3,ybase+[2,4,6])=data(1:3); %the channels will be flipped when using imshow
            activations(xbase+2,ybase+[1,3,5,7])=data(4:7);
            activations(xbase+1,ybase+[2,4,6])=data(8:10);
            
            %ploting the emitters and detectors placed on left hemisphere
            plot(xbase+[1,5,3,7],ybase+[1,1,3,3],'ks','MarkerFaceColor','k','MarkerSize', 10)
            hold on
            plot(xbase+[3,7,1,5],ybase+[1,1,3,3],'ro','MarkerFaceColor','r','MarkerSize', 10)

            x = [1:8+xbase];
            y = [1:4+ybase];  
            %selecting the contour levels between the min and max
            %activations  
            levels=linspace(min(min(activations)),max(max(activations)),10); 
            contour(x,y,activations,levels,'ShowText','off','LineWidth',1.5) %narrow lines means negative levels. 
            contourcmap('hsv', levels,'Colorbar', 'on', ...
           'Location', 'horizontal', ...
           'TitleString', strcat('The 10 levels of magnitude are between ',num2str(min(min(data)))...
                ,' to ',num2str(max(max(data)))));
            axis equal                                % set axis units to be the same size
            axis square off
            imgh=imshow(repmat(Infimg_left,1,1,3), 'XData', [-11 14], 'YData', [-9 15]);
            imgh.AlphaData = .2;
            hold off
            xlim([0.77 9.716])
            ylim([0.8543 7.795])
            title(strcat('Spatial Component ', num2str(index2cntL(r,2))));
        subplot(1,3,3);% plotting subject
             hold on;
            starti=1;
             subjectVal(1:subj_size)=0;
            for conditionI=1:length(conditions)
                endi=starti+length(find(condition_ind==conditionI))-1;
                subjectVal(starti:endi,1)=SumSubjectWeightsL(index2cntL(r,1),index2cntL(r,2),starti:endi);
                bar(conditionI, mean(subjectVal(starti:endi)));
                starti=endi+1;
            end
           xlabel('condition #','Fontsize',18,'FontName','Times New Roman')
            ylabel('Signature','Fontsize',18,'FontName','Times New Roman')
            set(gca,'box','off', 'FontSize', 16);
            legend(conditions_title,'Location','best');
        F    = getframe(FigH);
        imwrite(F.cdata, fullfile(inter_dir,...
                strcat('left mean',num2str(r),'.png')), 'png')
        close 

    end
end
% plotting the significant entity, action, and interaction components on the right hemisphere 
for r=1:r_components_r  
    %Right
    if T_anova_right{r,2}<0.05 %plot significant entity 
        FigH = figure('Position', get(0, 'Screensize'),'visible',false);
        subplot(1,3,1); %plotting temporal components
            plot(U_right{1}(:,index2cntR(r,1)),'LineWidth',2)
            title(strcat('Temporal Component ', num2str(index2cntR(r,1))))
            xlabel('Time (sec)','Fontsize',18,'FontName','Times New Roman')
            ylabel('Magnitude','Fontsize',18,'FontName','Times New Roman')
            xticks(0:fs*2:noSamples_time)
            xticklabels(0:2:noSamples_time/fs)
            set(gca,'box','off', 'FontSize', 16);
        subplot(1,3,2); %plotting right spatial components 
            data=U_right{2}(:,index2cntR(r,2));

            activations(xbase+1,ybase+[6,4,2])=data(8:10);
            activations(xbase+2,ybase+[7,5,3,1])=data(4:7);
            activations(xbase+3,ybase+[6,4,2])=data(1:3); %the channels will be plotted as is when using imshow

            %ploting the emitters and detectors placed on left hemisphere
            plot(xbase+[1,5,3,7],ybase+[1,1,3,3],'ro','MarkerFaceColor','k','MarkerSize', 10)
            hold on
            plot(xbase+[3,7,1,5],ybase+[1,1,3,3],'ks','MarkerFaceColor','r','MarkerSize', 10)

            x = [1:8+xbase];
            y = [1:4+ybase];  
            %selecting the contour levels between the min and max
            %activations  
            levels=linspace(min(min(activations)),max(max(activations)),10); 
            contour(x,y,activations,levels,'ShowText','off','LineWidth',1.5) %narrow lines means negative levels. 
            contourcmap('hsv', levels,'Colorbar', 'on', ...
           'Location', 'horizontal', ...
           'TitleString', strcat('The 10 levels of magnitude are between ',num2str(min(min(activations)))...
                ,' to ',num2str(max(max(activations)))));
            axis equal                                % set axis units to be the same size
            axis square off
            imgh=imshow(repmat(Infimg_right,1,1,3), 'XData', [-4 21], 'YData', [-9 15]);
            imgh.AlphaData = .2;
            hold off
            xlim([0.77 9.716])
            ylim([0.8543 7.795]) 
            title(strcat('Spatial Component ', num2str(index2cntL(r,2))));  
         subplot(1,3,3);% plotting subject
             hold on;
            starti=1;
             subjectVal(1:subj_size)=0;
            for conditionI=1:length(conditions)
                endi=starti+length(find(condition_ind==conditionI))-1;
                subjectVal(starti:endi,1)=SumSubjectWeightsR(index2cntL(r,1),index2cntL(r,2),starti:endi);
                bar(conditionI, mean(subjectVal(starti:endi)));
                starti=endi+1;
            end
            xlabel('condition #','Fontsize',18,'FontName','Times New Roman')
            ylabel('Signature','Fontsize',18,'FontName','Times New Roman')
            set(gca,'box','off', 'FontSize', 16);
            legend(conditions_title,'Location','best');
         
        F    = getframe(FigH);
        imwrite(F.cdata, fullfile(entity_dir,...
                strcat('right mean',num2str(r),'.png')), 'png')
        close  
                
    end
    
    if T_anova_right{r,4}<0.05 %plot significant action for right hemisphere
        FigH = figure('Position', get(0, 'Screensize'),'visible',false);
        subplot(1,3,1); %plotting temporal components
            plot(U_right{1}(:,index2cntR(r,1)),'LineWidth',2)
            title(strcat('Temporal Component ', num2str(index2cntR(r,1))))
            xlabel('Time (sec)','Fontsize',18,'FontName','Times New Roman')
            ylabel('Magnitude','Fontsize',18,'FontName','Times New Roman')
            xticks(0:fs*2:noSamples_time)
            xticklabels(0:2:noSamples_time/fs)
            set(gca,'box','off', 'FontSize', 16);
        subplot(1,3,2); %plotting right spatial components 
            data=U_right{2}(:,index2cntR(r,2));

            activations(xbase+1,ybase+[6,4,2])=data(8:10);
            activations(xbase+2,ybase+[7,5,3,1])=data(4:7);
            activations(xbase+3,ybase+[6,4,2])=data(1:3); %the channels will be plotted as is when using imshow

            %ploting the emitters and detectors placed on left hemisphere
            plot(xbase+[1,5,3,7],ybase+[1,1,3,3],'ro','MarkerFaceColor','k','MarkerSize', 10)
            hold on
            plot(xbase+[3,7,1,5],ybase+[1,1,3,3],'ks','MarkerFaceColor','r','MarkerSize', 10)

            x = [1:8+xbase];
            y = [1:4+ybase];  
            %selecting the contour levels between the min and max
            %activations  
            levels=linspace(min(min(activations)),max(max(activations)),10);  
            contour(x,y,activations,levels,'ShowText','off','LineWidth',1.5) %narrow lines means negative levels. 
            contourcmap('hsv', levels,'Colorbar', 'on', ...
           'Location', 'horizontal', ...
           'TitleString', strcat('The 10 levels of magnitude are between ',num2str(min(min(activations)))...
                ,' to ',num2str(max(max(activations)))));
            axis equal                                % set axis units to be the same size 
            axis square off
            imgh=imshow(repmat(Infimg_right,1,1,3), 'XData', [-4 21], 'YData', [-9 15]);
            imgh.AlphaData = .2;
            hold off 
            xlim([0.77 9.716])
            ylim([0.8543 7.795]) 
            title(strcat('Spatial Component ', num2str(index2cntL(r,2))));  
         subplot(1,3,3);%plotting subject
             hold on;
            starti=1;
             subjectVal(1:subj_size)=0;
            for conditionI=1:length(conditions)
                endi=starti+length(find(condition_ind==conditionI))-1;
                subjectVal(starti:endi,1)=SumSubjectWeightsR(index2cntL(r,1),index2cntL(r,2),starti:endi);
                bar(conditionI, mean(subjectVal(starti:endi)));
                starti=endi+1;
            end
             xlabel('condition #','Fontsize',18,'FontName','Times New Roman')
            ylabel('Signature','Fontsize',18,'FontName','Times New Roman')
            set(gca,'box','off', 'FontSize', 16);
            legend(conditions_title,'Location','best');
                    
        F    = getframe(FigH);
        imwrite(F.cdata, fullfile(action_dir,...
                strcat('right mean',num2str(r),'.png')), 'png')
        close 
    end
    
    if T_anova_right{r,6}<0.05 %plot significant interactions for the right hemisphere 
       FigH = figure('Position', get(0, 'Screensize'),'visible',false);
        subplot(1,3,1); %plotting temporal components
            plot(U_right{1}(:,index2cntR(r,1)),'LineWidth',2)
            title(strcat('Temporal Component ', num2str(index2cntR(r,1))))
            xlabel('Time (sec)','Fontsize',18,'FontName','Times New Roman')
            ylabel('Magnitude','Fontsize',18,'FontName','Times New Roman')
            xticks(0:fs*2:noSamples_time)
            xticklabels(0:2:noSamples_time/fs)
            set(gca,'box','off', 'FontSize', 16);

        subplot(1,3,2); %plotting right spatial components 
            data=U_right{2}(:,index2cntR(r,2));

            activations(xbase+1,ybase+[6,4,2])=data(8:10);
            activations(xbase+2,ybase+[7,5,3,1])=data(4:7);
            activations(xbase+3,ybase+[6,4,2])=data(1:3); %the channels will be plotted as is when using imshow

            %ploting the emitters and detectors placed on left hemisphere
            plot(xbase+[1,5,3,7],ybase+[1,1,3,3],'ro','MarkerFaceColor','k','MarkerSize', 10)
            hold on
            plot(xbase+[3,7,1,5],ybase+[1,1,3,3],'ks','MarkerFaceColor','r','MarkerSize', 10)

            x = [1:8+xbase];
            y = [1:4+ybase];  
            %selecting the contour levels between the min and max
            %activations  
            levels=linspace(min(min(activations)),max(max(activations)),10); 
            contour(x,y,activations,levels,'ShowText','off','LineWidth',1.5) %narrow lines means negative levels. 
            contourcmap('hsv', levels,'Colorbar', 'on', ...
           'Location', 'horizontal', ...
           'TitleString', strcat('The 10 levels of magnitude are between ',num2str(min(min(activations)))...
                ,' to ',num2str(max(max(activations)))));
            axis equal                                % set axis units to be the same size
            axis square off
            imgh=imshow(repmat(Infimg_right,1,1,3), 'XData', [-4 21], 'YData', [-9 15]);
            imgh.AlphaData = .2;
            hold off
            xlim([0.77 9.716])
            ylim([0.8543 7.795]) 
            title(strcat('Spatial Component ', num2str(index2cntL(r,2)))); 

         subplot(1,3,3); %plotting subject 
                hold on;
                starti=1;
                subjectVal(1:subj_size)=0;
                for conditionI=1:length(conditions)
                    endi=starti+length(find(condition_ind==conditionI))-1;
                    subjectVal(starti:endi,1)=SumSubjectWeightsR(index2cntL(r,1),index2cntL(r,2),starti:endi);
                    bar(conditionI, mean(subjectVal(starti:endi)));
                    starti=endi+1;
                end
                xlabel('condition #','Fontsize',18,'FontName','Times New Roman')
                ylabel('Signature','Fontsize',18,'FontName','Times New Roman')
                set(gca,'box','off', 'FontSize', 16);
                legend(conditions_title,'Location','best');

        F    = getframe(FigH);
        imwrite(F.cdata, fullfile(inter_dir,...
        strcat('right mean',num2str(r),'.png')), 'png')
        close 
    end
end


%putting zeros when the temporal and spatial component isnt significant in the
%temporal component by spatial component by subjects
ZeroL =SumSubjectWeightsL(:,:,:); %left hemi 
ZeroR =SumSubjectWeightsR(:,:,:); %right hemi
ZeroL(:,:,:) =0; %left hemi 
ZeroR(:,:,:) =0; %right hemi

sig_component_entity_l=find(T_anova_left{:,2}<0.05); %Significant components based on the entity factor
sig_component_action_l=find(T_anova_left{:,4}<0.05); %Significant components based on the action factor
sig_component_int_l=find(T_anova_left{:,6}<0.05); %Significant components based on the interaction

sig_component_entity_r=find(T_anova_right{:,2}<0.05); %Significant components based on the entity factor
sig_component_action_r=find(T_anova_right{:,4}<0.05); %Significant components based on the action factor
sig_component_int_r=find(T_anova_right{:,6}<0.05); %Significant components based on the interaction


% Excluding Irrelevant Components
time_exclusion = 15; % time point (in seconds) that you want to used for cutoff 
time_exclu_adjust = (time_exclusion+2)*fs; % convert to time points 
epoch_exclude = .01; %excludes temporal components if the mean magnitude during stimulus presentation is below this values  
temp_freq_exclude =.1; % mean frequency (in Hz) should not be above this values  


%zeroing out core for only significant entity, action and interaction for
%left and right hemisphers
 %count of components that were stat sig and passed exclusion parameters 
c_entity_array_l=zeros(length(sig_component_entity_l),3); 
c_action_array_l=zeros(length(sig_component_action_l),3); 
c_inter_array_l=zeros(length(sig_component_int_l),3); 
c_entity_array_r=zeros(length(sig_component_entity_r),3); 
c_action_array_r=zeros(length(sig_component_action_r),3); 
c_inter_array_r=zeros(length(sig_component_int_r),3); 

%finding only significant components that passed the exclusion criteria 
%left hemisphere 
only_sig_component_entity_l=ZeroL; %main effect of entity 
    for i = 1:length(sig_component_entity_l)
        r=sig_component_entity_l(i); 
        if (mean(abs(U_left{1}((2*fs):time_exclu_adjust,index2cntL(r,1))))>epoch_exclude)==1 ...
            && (meanfreq(U_left{1}(:,index2cntL(r,1)),fs)<temp_freq_exclude)==1
                only_sig_component_entity_l(index2cntL(r,1),index2cntL(r,2),:)=Sol_left(index2cntL(r,1),index2cntL(r,2),:);
                c_entity_array_l(i,1:2) =  index2cntL(r,:); %temporal and spatial component number that passed selection criterion 
                c_entity_array_l(i,3)=r; %the component number after combining temporal and spatial componesn
        end
    end
only_sig_component_action_l=ZeroL; % main effect of action sequence 
      for i = 1:length(sig_component_action_l)
        r=sig_component_action_l(i); %it could be i = 1:2 (left7 and left9)
         if (mean(abs(U_left{1}((2*fs):time_exclu_adjust,index2cntL(r,1))))>epoch_exclude)==1 ...
            && (meanfreq(U_left{1}(:,index2cntL(r,1)),fs)<temp_freq_exclude)==1
                only_sig_component_action_l(index2cntL(r,1),index2cntL(r,2),:)=Sol_left(index2cntL(r,1),index2cntL(r,2),:);
                c_action_array_l(i,1:2) =  index2cntL(r,:);
                c_action_array_l(i,3)=r; 
         end
      end
 only_sig_component_int_l=ZeroL; %interaction effect 
    for i = 1:length(sig_component_int_l)
        r=sig_component_int_l(i); %%% it coule be i=1 (left13)
         if (mean(abs(U_left{1}((2*fs):time_exclu_adjust,index2cntL(r,1))))>epoch_exclude)==1 ...
            && (meanfreq(U_left{1}(:,index2cntL(r,1)),fs)<temp_freq_exclude)==1
                only_sig_component_int_l(index2cntL(r,1),index2cntL(r,2),:)=Sol_left(index2cntL(r,1),index2cntL(r,2),:);
                c_inter_array_l(i,1:2) =  index2cntL(r,:);
                 c_inter_array_l(i,3)=r; 
         end
    end 
%right hemisphere
only_sig_component_entity_r=ZeroR;
    for i = 1:length(sig_component_entity_r) 
        r=sig_component_entity_r(i);  
         if (mean(abs(U_right{1}((2*fs):time_exclu_adjust,index2cntR(r,1))))>epoch_exclude)==1 ...
            && (meanfreq(U_right{1}(:,index2cntR(r,1)),fs)<temp_freq_exclude)==1
             only_sig_component_entity_r(index2cntR(r,1),index2cntR(r,2),:)=Sol_right(index2cntR(r,1),index2cntR(r,2),:);
             c_entity_array_r(i,1:2) =  index2cntR(r,:);
             c_entity_array_r(i,3)=r; 
         end
    end
% none of the significant components were included for right hemisphere action
% all the code for right hemisphere action are commented out. You will have to uncomment it to make it work 
% only_sig_component_action_r=ZeroR;
    for i = 1:length(sig_component_action_r)
        r=sig_component_action_r(i); 
        if (mean(abs(U_right{1}((2*fs):time_exclu_adjust,index2cntR(r,1))))>epoch_exclude)==1 ...
            && (meanfreq(U_right{1}(:,index2cntR(r,1)),fs)<temp_freq_exclude)==1
                only_sig_component_action_r(index2cntR(r,1),index2cntR(r,2),:)=Sol_right(index2cntL(r,1),index2cntR(r,2),:);
                c_action_array_r(i,1:2)=  index2cntR(r,:);
                c_action_array_r(i,3)=r; 
        end
    end
 only_sig_component_int_r=ZeroR;
    for i = 1:length(sig_component_int_r)
        r=sig_component_int_r(i); 
        if (mean(abs(U_right{1}((2*fs):time_exclu_adjust,index2cntR(r,1))))>epoch_exclude)==1 ...
            && (meanfreq(U_right{1}(:,index2cntR(r,1)),fs)<temp_freq_exclude)==1
            only_sig_component_int_r(index2cntR(r,1),index2cntR(r,2),:)=Sol_right(index2cntR(r,1),index2cntR(r,2),:);
            c_inter_array_r(i,1:2) =  index2cntR(r,:);
            c_inter_array_r(i,3)=r; 
        end 
    end
   
% subject factors times zeroed out core 
%left hemisphere
  Subj_only_sig_component_entity_l(1:temp_r_l,1:spac_r_l,1:subj_size,1:subj_size)=0;
for t = 1:temp_r_l %number of temporal components
   for c = 1:spac_r_l %number of spatial components
      for s = 1:subj_size %number of subjects component
          for i = 1:subj_size %subject number
             Subj_only_sig_component_entity_l(t,c,s,i) = only_sig_component_entity_l(t,c,s)* U_left{1,3}(i,s);
          end
      end
   end
end
 Subj_only_sig_component_action_l(1:temp_r_r,1:spac_r_l,1:subj_size,1:subj_size)=0;
for t = 1:temp_r_l %number of temporal components
   for c = 1:spac_r_l%number of spatial components
      for s = 1:subj_size %number of subjects component
          for i = 1:subj_size %subject number
             Subj_only_sig_component_action_l(t,c,s,i) = only_sig_component_action_l(t,c,s)* U_left{1,3}(i,s);
          end
      end
   end
end
 Subj_only_sig_component_int_l(1:temp_r_l,1:spac_r_l,1:subj_size,1:subj_size)=0;
for t = 1:temp_r_l %number of temporal components
   for c = 1:spac_r_l %number of spatial components
      for s = 1:subj_size %number of subjects component
          for i = 1:subj_size %subject number
             Subj_only_sig_component_int_l(t,c,s,i) = only_sig_component_int_l(t,c,s)* U_left{1,3}(i,s);
          end
      end
   end
end
%right hemisphere
 Subj_only_sig_component_entity_r(1:temp_r_r,1:spac_r_r,1:subj_size,1:subj_size)=0;
for t = 1:temp_r_r %number of temporal components
   for c = 1:spac_r_r %number of spatial components
      for s = 1:subj_size %number of subjects component
          for i = 1:subj_size %subject number
             Subj_only_sig_component_entity_r(t,c,s,i) = only_sig_component_entity_r(t,c,s)* U_right{1,3}(i,s);
          end
      end
   end
end
%  Subj_only_sig_component_action_r(1:temp_r_r,1:spac_r_r,1:subj_size,1:subj_size)=0;
% for t = 1:temp_r_r %number of temporal components
%    for c = 1:spac_r_r %number of spatial components
%       for s = 1:subj_size %number of subjects component
%           for i = 1:subj_size %subject number
%              Subj_only_sig_component_action_r(t,c,s,i) = only_sig_component_action_r(t,c,s)* U_right{1,3}(i,s);
%           end
%       end
%    end
% end
 Subj_only_sig_component_int_r(1:temp_r_r,1:spac_r_r,1:subj_size,1:subj_size)=0;
for t = 1:temp_r_r %number of temporal components
   for c = 1:spac_r_r %number of spatial components
      for s = 1:subj_size %number of subjects component
          for i = 1:subj_size %subject number
             Subj_only_sig_component_int_r(t,c,s,i) = only_sig_component_int_r(t,c,s)* U_right{1,3}(i,s);
          end
      end
   end
end

%avg/summing the components on subject number so this makes a summed zeroed 'core' which is temporal weight by spatial weight by subject number tensor 
%left hemisphere
SumSubjectWeightsL_entity(1:temp_r_l,1:spac_r_l,1:subj_size)=0;
for t = 1:temp_r_l %number of temporal components
   for c = 1:spac_r_l %number of spatial components
      for i = 1:subj_size %subjects number
          SumSubjectWeightsL_entity(t,c,i) = sum(Subj_only_sig_component_entity_l(t,c,:,i));
      end
   end
end
SumSubjectWeightsL_action(1:temp_r_l,1:spac_r_l,1:subj_size)=0;
for t = 1:temp_r_l %number of temporal components
   for c = 1:spac_r_l %number of spatial components
      for i = 1:subj_size %subjects number
          SumSubjectWeightsL_action(t,c,i) = sum(Subj_only_sig_component_action_l(t,c,:,i));
      end
   end
end
SumSubjectWeightsL_Int(1:temp_r_l,1:spac_r_l,1:subj_size)=0;
for t = 1:temp_r_l %number of temporal components
   for c = 1:spac_r_l %number of spatial components
      for i = 1:subj_size %subjects number
          SumSubjectWeightsL_Int(t,c,i) = sum(Subj_only_sig_component_int_l(t,c,:,i));
      end
   end
end
%avg/summing the components on subject number right hemisphere
SumSubjectWeightsR_entity(1:temp_r_r,1:spac_r_r,1:subj_size)=0;
for t = 1:temp_r_r %number of temporal components
   for c = 1:spac_r_r %number of spatial components
      for i = 1:subj_size %subjects number
          SumSubjectWeightsR_entity(t,c,i) = sum(Subj_only_sig_component_entity_r(t,c,:,i));
      end
   end
end
% SumSubjectWeightsR_action(1:temp_r_r,1:spac_r_r,1:subj_size)=0;
% for t = 1:temp_r_r %number of temporal components
%    for c = 1:spac_r_r %number of spatial components
%       for i = 1:subj_size %subjects number
%           SumSubjectWeightsR_action(t,c,i) = sum(Subj_only_sig_component_action_r(t,c,:,i));
%       end
%    end
% end
SumSubjectWeightsR_Int(1:temp_r_r,1:spac_r_r,1:subj_size)=0;
for t = 1:temp_r_r %number of temporal components
   for c = 1:spac_r_r %number of spatial components
      for i = 1:subj_size %subjects number
          SumSubjectWeightsR_Int(t,c,i) = sum(Subj_only_sig_component_int_r(t,c,:,i));
      end
   end
end

%making a 8 by 10 core by summing up all the subjects for left hemisphere
%makes individual plots for each condition
%left entity effect with mean for each condition 
temp_by_spatial_core_entity_l_hum_fun=ZeroL(:,:,1);
    for i = 1:length(sig_component_entity_l)
        r=sig_component_entity_l(i); 
        temp_by_spatial_core_entity_l_hum_fun(index2cntL(r,1),index2cntL(r,2))=...
            mean(SumSubjectWeightsL_entity(index2cntL(r,1),index2cntL(r,2),(condition_ind==1)));
    end
temp_by_spatial_core_entity_l_hum_nonfun=ZeroL(:,:,1);
    for i = 1:length(sig_component_entity_l)
        r=sig_component_entity_l(i); 
        temp_by_spatial_core_entity_l_hum_nonfun(index2cntL(r,1),index2cntL(r,2))=...
            mean(SumSubjectWeightsL_entity(index2cntL(r,1),index2cntL(r,2),(condition_ind==2)));
    end
temp_by_spatial_core_entity_l_claw_fun=ZeroL(:,:,1);
    for i = 1:length(sig_component_entity_l)
        r=sig_component_entity_l(i); 
        temp_by_spatial_core_entity_l_claw_fun(index2cntL(r,1),index2cntL(r,2))=...
            mean(SumSubjectWeightsL_entity(index2cntL(r,1),index2cntL(r,2),(condition_ind==3)));
    end
temp_by_spatial_core_entity_l_claw_nonfun=ZeroL(:,:,1);
    for i = 1:length(sig_component_entity_l)
        r=sig_component_entity_l(i); 
        temp_by_spatial_core_entity_l_claw_nonfun(index2cntL(r,1),index2cntL(r,2))=...
            mean(SumSubjectWeightsL_entity(index2cntL(r,1),index2cntL(r,2),(condition_ind==4)));
    end    
% left action effect with mean for each condition   
temp_by_spatial_core_action_l_hum_fun=ZeroL(:,:,1);
    for i = 1:length(sig_component_action_l)
        r=sig_component_action_l(i); 
        temp_by_spatial_core_action_l_hum_fun(index2cntL(r,1),index2cntL(r,2))=...
            mean(SumSubjectWeightsL_action(index2cntL(r,1),index2cntL(r,2),(condition_ind==1)));
    end
temp_by_spatial_core_action_l_hum_nonfun=ZeroL(:,:,1);
    for i = 1:length(sig_component_action_l)
        r=sig_component_action_l(i); 
        temp_by_spatial_core_action_l_hum_nonfun(index2cntL(r,1),index2cntL(r,2))=...
            mean(SumSubjectWeightsL_action(index2cntL(r,1),index2cntL(r,2),(condition_ind==2)));
    end
temp_by_spatial_core_action_l_claw_fun=ZeroL(:,:,1);
    for i = 1:length(sig_component_action_l)
        r=sig_component_action_l(i); 
        temp_by_spatial_core_action_l_claw_fun(index2cntL(r,1),index2cntL(r,2))=...
            mean(SumSubjectWeightsL_action(index2cntL(r,1),index2cntL(r,2),(condition_ind==3)));
    end
temp_by_spatial_core_action_l_claw_nonfun=ZeroL(:,:,1);
    for i = 1:length(sig_component_action_l)
        r=sig_component_action_l(i); 
        temp_by_spatial_core_action_l_claw_nonfun(index2cntL(r,1),index2cntL(r,2))=...
            mean(SumSubjectWeightsL_action(index2cntL(r,1),index2cntL(r,2),(condition_ind==4)));
    end
% left interaction effect with meant for each condition  
temp_by_spatial_core_int_l_hum_fun=ZeroL(:,:,1);
    for i = 1:length(sig_component_int_l)
        r=sig_component_int_l(i); 
       temp_by_spatial_core_int_l_hum_fun(index2cntL(r,1),index2cntL(r,2))=...
           mean(SumSubjectWeightsL_Int(index2cntL(r,1),index2cntL(r,2),(condition_ind==1)));
    end
temp_by_spatial_core_int_l_hum_nonfun=ZeroL(:,:,1);
    for i = 1:length(sig_component_int_l)
        r=sig_component_int_l(i); 
       temp_by_spatial_core_int_l_hum_nonfun(index2cntL(r,1),index2cntL(r,2))=...
           mean(SumSubjectWeightsL_Int(index2cntL(r,1),index2cntL(r,2),(condition_ind==2)));
    end    
 temp_by_spatial_core_int_l_claw_fun=ZeroL(:,:,1);
    for i = 1:length(sig_component_int_l)
        r=sig_component_int_l(i); 
       temp_by_spatial_core_int_l_claw_fun(index2cntL(r,1),index2cntL(r,2))=...
           mean(SumSubjectWeightsL_Int(index2cntL(r,1),index2cntL(r,2),(condition_ind==3)));
    end
  temp_by_spatial_core_int_l_claw_nonfun=ZeroL(:,:,1);
    for i = 1:length(sig_component_int_l)
        r=sig_component_int_l(i); 
       temp_by_spatial_core_int_l_claw_nonfun(index2cntL(r,1),index2cntL(r,2))=...
           mean(SumSubjectWeightsL_Int(index2cntL(r,1),index2cntL(r,2),(condition_ind==4)));
    end   
%summing up all the subjects for right hemisphere
%right entity  
temp_by_spatial_core_entity_r_hum_fun=ZeroR(:,:,1);
    for i = 1:length(sig_component_entity_r)
        r=sig_component_entity_r(i); 
        temp_by_spatial_core_entity_r_hum_fun(index2cntR(r,1),index2cntR(r,2))=...
            mean(SumSubjectWeightsR_entity(index2cntR(r,1),index2cntR(r,2),(condition_ind==1)));
    end
temp_by_spatial_core_entity_r_hum_nonfun=ZeroR(:,:,1);
    for i = 1:length(sig_component_entity_r)
        r=sig_component_entity_r(i); 
        temp_by_spatial_core_entity_r_hum_nonfun(index2cntR(r,1),index2cntR(r,2))=...
            mean(SumSubjectWeightsR_entity(index2cntR(r,1),index2cntR(r,2),(condition_ind==2)));
    end
temp_by_spatial_core_entity_r_claw_fun=ZeroR(:,:,1);
    for i = 1:length(sig_component_entity_r)
        r=sig_component_entity_r(i); 
        temp_by_spatial_core_entity_r_claw_fun(index2cntR(r,1),index2cntR(r,2))=...
            mean(SumSubjectWeightsR_entity(index2cntR(r,1),index2cntR(r,2),(condition_ind==3)));
    end    
temp_by_spatial_core_entity_r_claw_nonfun=ZeroR(:,:,1);
    for i = 1:length(sig_component_entity_r)
        r=sig_component_entity_r(i); 
        temp_by_spatial_core_entity_r_claw_nonfun(index2cntR(r,1),index2cntR(r,2))=...
            mean(SumSubjectWeightsR_entity(index2cntR(r,1),index2cntR(r,2),(condition_ind==4)));
    end   
% right action    
% temp_by_spatial_core_action_r_hum_fun=ZeroR(:,:,1);
%     for i = 1:length(sig_component_action_r)
%         r=sig_component_action_r(i); 
%         temp_by_spatial_core_action_r_hum_fun(index2cntR(r,1),index2cntR(r,2))=...
%             mean(SumSubjectWeightsR_action(index2cntR(r,1),index2cntR(r,2),(condition_ind==1)));
%     end
% temp_by_spatial_core_action_r_hum_nonfun=ZeroR(:,:,1);
%     for i = 1:length(sig_component_action_r)
%         r=sig_component_action_r(i); 
%         temp_by_spatial_core_action_r_hum_nonfun(index2cntR(r,1),index2cntR(r,2))=...
%             mean(SumSubjectWeightsR_action(index2cntR(r,1),index2cntR(r,2),(condition_ind==2)));
%     end
% temp_by_spatial_core_action_r_claw_fun=ZeroR(:,:,1);
%     for i = 1:length(sig_component_action_r)
%         r=sig_component_action_r(i); 
%         temp_by_spatial_core_action_r_claw_fun(index2cntR(r,1),index2cntR(r,2))=...
%             mean(SumSubjectWeightsR_action(index2cntR(r,1),index2cntR(r,2),(condition_ind==3)));
%     end
% temp_by_spatial_core_action_r_claw_nonfun=ZeroR(:,:,1);
%     for i = 1:length(sig_component_action_r)
%         r=sig_component_action_r(i); 
%         temp_by_spatial_core_action_r_claw_nonfun(index2cntR(r,1),index2cntR(r,2))=...
%             mean(SumSubjectWeightsR_action(index2cntR(r,1),index2cntR(r,2),(condition_ind==4)));
%     end    
%     
% right interaction     
temp_by_spatial_core_int_r_hum_fun=ZeroR(:,:,1);
    for i = 1:length(sig_component_int_r)
        r=sig_component_int_r(i); 
        temp_by_spatial_core_int_r_hum_fun(index2cntR(r,1),index2cntR(r,2))=...
            mean(SumSubjectWeightsR_Int(index2cntR(r,1),index2cntR(r,2),(condition_ind==1)));
    end
 temp_by_spatial_core_int_r_hum_nonfun=ZeroR(:,:,1);
    for i = 1:length(sig_component_int_r)
        r=sig_component_int_r(i); 
        temp_by_spatial_core_int_r_hum_nonfun(index2cntR(r,1),index2cntR(r,2))=...
            mean(SumSubjectWeightsR_Int(index2cntR(r,1),index2cntR(r,2),(condition_ind==2)));
    end   
temp_by_spatial_core_int_r_claw_fun=ZeroR(:,:,1);
    for i = 1:length(sig_component_int_r)
        r=sig_component_int_r(i); 
        temp_by_spatial_core_int_r_claw_fun(index2cntR(r,1),index2cntR(r,2))=...
            mean(SumSubjectWeightsR_Int(index2cntR(r,1),index2cntR(r,2),(condition_ind==3)));
    end    
temp_by_spatial_core_int_r_claw_nonfun=ZeroR(:,:,1);
    for i = 1:length(sig_component_int_r)
        r=sig_component_int_r(i); 
        temp_by_spatial_core_int_r_claw_nonfun(index2cntR(r,1),index2cntR(r,2))=...
            mean(SumSubjectWeightsR_Int(index2cntR(r,1),index2cntR(r,2),(condition_ind==4)));
    end    



% finding temporal profile 
%left entity
    timept_only_sig_component_entity_l_hum_fun= sum(U_left{1,1}* temp_by_spatial_core_entity_l_hum_fun,2);
    timept_only_sig_component_entity_l_hum_nonfun= sum(U_left{1,1}* temp_by_spatial_core_entity_l_hum_nonfun,2);
    timept_only_sig_component_entity_l_claw_fun= sum(U_left{1,1}* temp_by_spatial_core_entity_l_claw_fun,2);
    timept_only_sig_component_entity_l_claw_nonfun= sum(U_left{1,1}* temp_by_spatial_core_entity_l_claw_nonfun,2);
%left action 
    timept_only_sig_component_action_l_hum_fun=sum(U_left{1,1}* temp_by_spatial_core_action_l_hum_fun,2);
    timept_only_sig_component_action_l_hum_nonfun=sum(U_left{1,1}* temp_by_spatial_core_action_l_hum_nonfun,2);
    timept_only_sig_component_action_l_claw_fun=sum(U_left{1,1}* temp_by_spatial_core_action_l_claw_fun,2);
    timept_only_sig_component_action_l_claw_nonfun=sum(U_left{1,1}* temp_by_spatial_core_action_l_claw_nonfun,2);
%left interaction 
    timept_only_sig_component_int_l_hum_fun = sum(U_left{1,1}* temp_by_spatial_core_int_l_hum_fun,2);
    timept_only_sig_component_int_l_hum_nonfun = sum(U_left{1,1}* temp_by_spatial_core_int_l_hum_nonfun,2);
    timept_only_sig_component_int_l_claw_fun = sum(U_left{1,1}* temp_by_spatial_core_int_l_claw_fun,2); 
    timept_only_sig_component_int_l_claw_nonfun = sum(U_left{1,1}* temp_by_spatial_core_int_l_claw_nonfun,2);  
%right entity
     timept_only_sig_component_entity_r_hum_fun = sum(U_right{1,1}*temp_by_spatial_core_entity_r_hum_fun,2);
     timept_only_sig_component_entity_r_hum_nonfun = sum(U_right{1,1}*temp_by_spatial_core_entity_r_hum_nonfun,2);
     timept_only_sig_component_entity_r_claw_fun = sum(U_right{1,1}*temp_by_spatial_core_entity_r_claw_fun,2);
     timept_only_sig_component_entity_r_claw_nonfun = sum(U_right{1,1}*temp_by_spatial_core_entity_r_claw_nonfun,2);
%right action     
%      timept_only_sig_component_action_r_hum_fun  = sum(U_right{1,1}*temp_by_spatial_core_action_r_hum_fun,2);
%      timept_only_sig_component_action_r_hum_nonfun  = sum(U_right{1,1}*temp_by_spatial_core_action_r_hum_nonfun,2);
%      timept_only_sig_component_action_r_claw_fun  = sum(U_right{1,1}*temp_by_spatial_core_action_r_claw_fun,2);
%      timept_only_sig_component_action_r_claw_nonfun  = sum(U_right{1,1}*temp_by_spatial_core_action_r_claw_nonfun,2);
 %right interaction 
     timept_only_sig_component_int_r_hum_fun = sum(U_right{1,1}* temp_by_spatial_core_int_r_hum_fun,2);
     timept_only_sig_component_int_r_hum_nonfun = sum(U_right{1,1}* temp_by_spatial_core_int_r_hum_nonfun,2);
     timept_only_sig_component_int_r_claw_fun = sum(U_right{1,1}* temp_by_spatial_core_int_r_claw_fun,2);
     timept_only_sig_component_int_r_claw_nonfun = sum(U_right{1,1}* temp_by_spatial_core_int_r_claw_nonfun,2);
     

%getting ROI
%left entity
    chan_only_sig_component_entity_l_hum_fun= sum(temp_by_spatial_core_entity_l_hum_fun*U_left{1,2}');
    chan_only_sig_component_entity_l_hum_nonfun= sum(temp_by_spatial_core_entity_l_hum_nonfun*U_left{1,2}'); 
    chan_only_sig_component_entity_l_claw_fun= sum(temp_by_spatial_core_entity_l_claw_fun*U_left{1,2}');
    chan_only_sig_component_entity_l_claw_nonfun= sum(temp_by_spatial_core_entity_l_claw_nonfun*U_left{1,2}');
%left action
    chan_only_sig_component_action_l_hum_fun= sum(temp_by_spatial_core_action_l_hum_fun*U_left{1,2}');
    chan_only_sig_component_action_l_hum_nonfun= sum(temp_by_spatial_core_action_l_hum_nonfun*U_left{1,2}');
    chan_only_sig_component_action_l_claw_fun= sum(temp_by_spatial_core_action_l_claw_fun*U_left{1,2}');
    chan_only_sig_component_action_l_claw_nonfun= sum(temp_by_spatial_core_action_l_claw_nonfun*U_left{1,2}');
%left interaction 
    chan_only_sig_component_int_l_hum_fun= sum(temp_by_spatial_core_int_l_hum_fun*U_left{1,2}');
    chan_only_sig_component_int_l_hum_nonfun= sum(temp_by_spatial_core_int_l_hum_nonfun*U_left{1,2}');
    chan_only_sig_component_int_l_claw_fun= sum(temp_by_spatial_core_int_l_claw_fun*U_left{1,2}');
    chan_only_sig_component_int_l_claw_nonfun= sum(temp_by_spatial_core_int_l_claw_nonfun*U_left{1,2}');
%right entity
    chan_only_sig_component_entity_r_hum_fun= sum(temp_by_spatial_core_entity_r_hum_fun*U_right{1,2}');
    chan_only_sig_component_entity_r_hum_nonfun= sum(temp_by_spatial_core_entity_r_hum_nonfun*U_right{1,2}');
    chan_only_sig_component_entity_r_claw_fun= sum(temp_by_spatial_core_entity_r_claw_fun*U_right{1,2}');
    chan_only_sig_component_entity_r_claw_nonfun= sum(temp_by_spatial_core_entity_r_claw_nonfun*U_right{1,2}');
%right action
%     chan_only_sig_component_action_r_hum_fun= sum(temp_by_spatial_core_action_r_hum_fun*U_right{1,2}');
%     chan_only_sig_component_action_r_hum_nonfun= sum(temp_by_spatial_core_action_r_hum_nonfun*U_right{1,2}');
%     chan_only_sig_component_action_r_claw_fun= sum(temp_by_spatial_core_action_r_claw_fun*U_right{1,2}');
%     chan_only_sig_component_action_r_claw_nonfun= sum(temp_by_spatial_core_action_r_claw_nonfun*U_right{1,2}');
%right interaction 
    chan_only_sig_component_int_r_hum_fun= sum(temp_by_spatial_core_int_r_hum_fun*U_right{1,2}');
    chan_only_sig_component_int_r_hum_nonfun= sum(temp_by_spatial_core_int_r_hum_nonfun*U_right{1,2}');
    chan_only_sig_component_int_r_claw_fun= sum(temp_by_spatial_core_int_r_claw_fun*U_right{1,2}');
    chan_only_sig_component_int_r_claw_nonfun= sum(temp_by_spatial_core_int_r_claw_nonfun*U_right{1,2}');
    
%normalizing values for plotting 
for i = 1:spac_r_l
    %left entity 
        chan_only_sig_component_entity_l_hum_fun_norm(i)= ...
            (chan_only_sig_component_entity_l_hum_fun(i)...
            -mean(chan_only_sig_component_entity_l_hum_fun))/...
            std(chan_only_sig_component_entity_l_hum_fun);
        chan_only_sig_component_entity_l_hum_nonfun_norm(i)= ...
            (chan_only_sig_component_entity_l_hum_nonfun(i)...
            -mean(chan_only_sig_component_entity_l_hum_nonfun))/...
            std(chan_only_sig_component_entity_l_hum_nonfun);
        chan_only_sig_component_entity_l_claw_fun_norm(i)= ...
            (chan_only_sig_component_entity_l_claw_fun(i)...
            -mean(chan_only_sig_component_entity_l_claw_fun))/...
            std(chan_only_sig_component_entity_l_claw_fun);
        chan_only_sig_component_entity_l_claw_nonfun_norm(i)= ...
            (chan_only_sig_component_entity_l_claw_nonfun(i)...
            -mean(chan_only_sig_component_entity_l_claw_nonfun))/...
            std(chan_only_sig_component_entity_l_claw_nonfun);
    %left action
        chan_only_sig_component_action_l_hum_fun_norm(i)= ...
            (chan_only_sig_component_action_l_hum_fun(i)...
            -mean(chan_only_sig_component_action_l_hum_fun))/...
            std(chan_only_sig_component_action_l_hum_fun);
        chan_only_sig_component_action_l_hum_nonfun_norm(i)= ...
            (chan_only_sig_component_action_l_hum_nonfun(i)...
            -mean(chan_only_sig_component_action_l_hum_nonfun))/...
            std(chan_only_sig_component_action_l_hum_nonfun);
        chan_only_sig_component_action_l_claw_fun_norm(i)= ...
            (chan_only_sig_component_action_l_claw_fun(i)...
            -mean(chan_only_sig_component_action_l_claw_fun))/...
            std(chan_only_sig_component_action_l_claw_fun);
        chan_only_sig_component_action_l_claw_nonfun_norm(i)= ...
            (chan_only_sig_component_action_l_claw_nonfun(i)...
            -mean(chan_only_sig_component_action_l_claw_nonfun))/...
            std(chan_only_sig_component_action_l_claw_nonfun);
    %left interaction 
        chan_only_sig_component_int_l_hum_fun_norm(i)= ...
            (chan_only_sig_component_int_l_hum_fun(i)...
            -mean(chan_only_sig_component_int_l_hum_fun))/...
            std(chan_only_sig_component_int_l_hum_fun);
        chan_only_sig_component_int_l_hum_nonfun_norm(i)= ...
            (chan_only_sig_component_int_l_hum_nonfun(i)...
            -mean(chan_only_sig_component_int_l_hum_nonfun))/...
            std(chan_only_sig_component_int_l_hum_nonfun);
        chan_only_sig_component_int_l_claw_fun_norm(i)= ...
            (chan_only_sig_component_int_l_claw_fun(i)...
            -mean(chan_only_sig_component_int_l_claw_fun))/...
            std(chan_only_sig_component_int_l_claw_fun);
        chan_only_sig_component_int_l_claw_nonfun_norm(i)= ...
            (chan_only_sig_component_int_l_claw_nonfun(i)...
            -mean(chan_only_sig_component_int_l_claw_nonfun))/...
            std(chan_only_sig_component_int_l_claw_nonfun);
end
for i = 1:spac_r_r
    %right entity
        chan_only_sig_component_entity_r_hum_fun_norm(i)= ...
            (chan_only_sig_component_entity_r_hum_fun(i)...
            -mean(chan_only_sig_component_entity_r_hum_fun))/...
            std(chan_only_sig_component_entity_r_hum_fun);
        chan_only_sig_component_entity_r_hum_nonfun_norm(i)= ...
            (chan_only_sig_component_entity_r_hum_nonfun(i)...
            -mean(chan_only_sig_component_entity_r_hum_nonfun))/...
            std(chan_only_sig_component_entity_r_hum_nonfun);
        chan_only_sig_component_entity_r_claw_fun_norm(i)= ...
            (chan_only_sig_component_entity_r_claw_fun(i)...
            -mean(chan_only_sig_component_entity_r_claw_fun))/...
            std(chan_only_sig_component_entity_r_claw_fun);
        chan_only_sig_component_entity_r_claw_nonfun_norm(i)= ...
            (chan_only_sig_component_entity_r_claw_nonfun(i)...
            -mean(chan_only_sig_component_entity_r_claw_nonfun))/...
            std(chan_only_sig_component_entity_r_claw_nonfun);
    %right action
%         chan_only_sig_component_action_r_hum_fun_norm(i)= ...
%             (chan_only_sig_component_action_r_hum_fun(i)...
%             -mean(chan_only_sig_component_action_r_hum_fun))/...
%             std(chan_only_sig_component_action_r_hum_fun);
%         chan_only_sig_component_action_r_hum_nonfun_norm(i)= ...
%             (chan_only_sig_component_action_r_hum_nonfun(i)...
%             -mean(chan_only_sig_component_action_r_hum_nonfun))/...
%             std(chan_only_sig_component_action_r_hum_nonfun);
%         chan_only_sig_component_action_r_claw_fun_norm(i)= ...
%             (chan_only_sig_component_action_r_claw_fun(i)...
%             -mean(chan_only_sig_component_action_r_claw_fun))/...
%             std(chan_only_sig_component_action_r_claw_fun);
%         chan_only_sig_component_action_r_claw_nonfun_norm(i)= ...
%             (chan_only_sig_component_action_r_claw_nonfun(i)...
%             -mean(chan_only_sig_component_action_r_claw_nonfun))/...
%             std(chan_only_sig_component_action_r_claw_nonfun);
    %right interaction
         chan_only_sig_component_int_r_hum_fun_norm(i)= ...
            (chan_only_sig_component_int_r_hum_fun(i)...
            -mean(chan_only_sig_component_int_r_hum_fun))/...
            std(chan_only_sig_component_int_r_hum_fun);
        chan_only_sig_component_int_r_hum_nonfun_norm(i)= ...
            (chan_only_sig_component_int_r_hum_nonfun(i)...
            -mean(chan_only_sig_component_int_r_hum_nonfun))/...
            std(chan_only_sig_component_int_r_hum_nonfun);
        chan_only_sig_component_int_r_claw_fun_norm(i)= ...
            (chan_only_sig_component_int_r_claw_fun(i)...
            -mean(chan_only_sig_component_int_r_claw_fun))/...
            std(chan_only_sig_component_int_r_claw_fun);
        chan_only_sig_component_int_r_claw_nonfun_norm(i)= ...
            (chan_only_sig_component_int_r_claw_nonfun(i)...
            -mean(chan_only_sig_component_int_r_claw_nonfun))/...
            std(chan_only_sig_component_int_r_claw_nonfun);
end 
    
    
%plotting significant temporal profile for each condition 
%left entity
TD_anova_plot(timept_only_sig_component_entity_l_hum_fun,...
    chan_only_sig_component_entity_l_hum_fun_norm,true,td_dir,entity_dir,...
    '  Human Hand, Function','Left_entity_Human_Hand_Function')  
TD_anova_plot(timept_only_sig_component_entity_l_hum_nonfun,...
    chan_only_sig_component_entity_l_hum_nonfun_norm,true,td_dir,entity_dir,...
    ' Human Hand, Non-function','Left_entity_Human Hand_Non-function') 
TD_anova_plot(timept_only_sig_component_entity_l_claw_fun,...
    chan_only_sig_component_entity_l_claw_fun_norm,true,td_dir,entity_dir,...
    ' Mechanical Hand, Function','Left_entity_Mechanical_Hand_Function')  
TD_anova_plot(timept_only_sig_component_entity_l_claw_nonfun,...
    chan_only_sig_component_entity_l_claw_nonfun_norm,true,td_dir,entity_dir,...
    ' Mechanical Hand, Non-function', 'Left_entity_Mechanical_Hand_Non-function')  
%left action
TD_anova_plot(timept_only_sig_component_action_l_hum_fun,...
    chan_only_sig_component_action_l_hum_fun_norm,true,td_dir,action_dir,...
    ' Human Hand, Function', 'Left_action_Human_Hand_Function') 
TD_anova_plot(timept_only_sig_component_action_l_hum_nonfun,...
    chan_only_sig_component_action_l_hum_nonfun_norm,true,td_dir,action_dir,...
    'Human Hand, Non-function','Left_action_Human_Hand_Non-function')  
TD_anova_plot(timept_only_sig_component_action_l_claw_fun,...
    chan_only_sig_component_action_l_claw_fun_norm,true,td_dir,action_dir,...
    ' Mechanical Hand, Function', 'Left_action_Mechanical_Hand_Function')   
TD_anova_plot(timept_only_sig_component_action_l_claw_nonfun,...
    chan_only_sig_component_action_l_claw_nonfun_norm,true,td_dir,action_dir,...
    ' Mechanical Hand, Non-function','Mechanical_Hand_Non_function')   
%left interaction
TD_anova_plot(timept_only_sig_component_int_l_hum_fun,...
    chan_only_sig_component_int_l_hum_fun_norm,true,td_dir,inter_dir,...
    ' Human Hand, Function', 'Left_Interaction_Human_Hand_Function')  
TD_anova_plot(timept_only_sig_component_int_l_hum_nonfun,...
    chan_only_sig_component_int_l_hum_nonfun_norm,true,td_dir,inter_dir,...
    'Human Hand, Non-function','Left_Interaction_Human_Hand_Non-function')  
TD_anova_plot(timept_only_sig_component_int_l_claw_fun,...
    chan_only_sig_component_int_l_claw_fun_norm,true,td_dir,inter_dir,...
    'Mechanical Hand, Function', 'Left_Interaction_Mechanical_Hand_Function')  
TD_anova_plot(timept_only_sig_component_int_l_claw_nonfun,...
    chan_only_sig_component_int_l_claw_nonfun_norm,true,td_dir,inter_dir,...
    'Mechanical Hand, Non-function', 'Left_Interaction_Mechanical_Hand_Non-function')   
%right entity
TD_anova_plot(timept_only_sig_component_entity_r_hum_fun,...
    chan_only_sig_component_entity_r_hum_fun_norm,false,td_dir,entity_dir,...
    'Human Hand, Function', 'Right_entity_Human_Hand_Function') 
TD_anova_plot(timept_only_sig_component_entity_r_hum_nonfun,...
    chan_only_sig_component_entity_r_hum_nonfun_norm,false,td_dir,entity_dir,...
    'Human Hand, Non-function','Right_entity_Human_Hand_Non-function') 
TD_anova_plot(timept_only_sig_component_entity_r_claw_fun,...
    chan_only_sig_component_entity_r_claw_fun_norm,false,td_dir,entity_dir,...
    'Mechanical Hand, Function', 'Right_entity_Mechanical_Hand_Function') 
TD_anova_plot(timept_only_sig_component_entity_r_claw_nonfun,...
    chan_only_sig_component_entity_r_claw_nonfun_norm,false,td_dir,entity_dir,...
    'Mechanical Hand, Non-function', 'Right_entity_Mechanical_Hand_Non-function') 
%right action
% TD_anova_plot(timept_only_sig_component_action_r_hum_fun,...
%     chan_only_sig_component_action_r_hum_fun_norm,false,td_dir,action_dir,...
%     ' Right action: Human Hand, Function', 'Right_action_Human_Hand_Function')
% TD_anova_plot(timept_only_sig_component_action_r_hum_nonfun,...
%     chan_only_sig_component_action_r_hum_nonfun_norm,false,td_dir,action_dir,...
%     ' Right action: Human Hand, Non-function', 'Right_action_Human_Hand_Non-function')
% TD_anova_plot(timept_only_sig_component_action_r_claw_fun,...
%     chan_only_sig_component_action_r_claw_fun_norm,false,td_dir,action_dir,...
%     ' Right action: Mechanical Hand, Function', 'Right_action_Mechanical_Hand_Function')
% TD_anova_plot(timept_only_sig_component_action_r_claw_nonfun,...
%     chan_only_sig_component_action_r_claw_nonfun_norm,false,td_dir,action_dir,...
%     ' Right action: Mechanical Hand, Non-function', 'Right_action_Mechanical_Hand_Non-function')
%right interaction
TD_anova_plot(timept_only_sig_component_int_r_hum_fun,...
    chan_only_sig_component_int_r_hum_fun,false,td_dir,inter_dir,...
    'Human Hand, Function', 'Right_Interaction_Human_Hand_Function')
TD_anova_plot(timept_only_sig_component_int_r_hum_nonfun,...
    chan_only_sig_component_int_r_hum_nonfun,false,td_dir,inter_dir,...
    'Human Hand, Non-function','Right_Interaction_Human_Hand_Non-function')
TD_anova_plot(timept_only_sig_component_int_r_claw_fun,...
    chan_only_sig_component_int_r_claw_fun,false,td_dir,inter_dir,...
    'Mechanical Hand, Function',  'Right_Interaction_Mechanical Hand_Function')
TD_anova_plot(timept_only_sig_component_int_r_claw_nonfun,...
    chan_only_sig_component_int_r_claw_nonfun,false,td_dir,inter_dir,...
    'Mechanical Hand, Non-function','Right_Interaction_Mechanical_Hand_Non-function')


save(fullfile(strcat(td_dir,'/anova'),'all variables'));

        







