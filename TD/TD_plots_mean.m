function TD_plots_mean(U,r_components,left,group_ind,groups,saveDir,fs)
    %Generate a figure for each CPD component contains. Plotting channels activations on the left and
    %right hemisphere. Taking group indeces to use in bar plots
    % Plotting the mean of the signatures of each group of subjects
    noSamples_time=length(U{1});
    noSamples_freq=length(U{2});
    xbase=1;ybase=1; activations=zeros(4+xbase,8+ybase);
    Infimg_left=imread('.\left pic.gif',1);
    Infimg_left=imresize(Infimg_left,0.2);%(Infimg,0.09);
    Infimg_right=imread('.\right pic.gif',1);
    Infimg_right=imresize(Infimg_right,0.5);%(Infimg,0.09);
    if exist(fullfile(saveDir,'figures'))~=7
        mkdir(fullfile(saveDir,'figures'));
    end
    
%plotting temporal component   
     for i=1:r_components(1)
        FigT = figure('Position', get(0, 'Screensize'),'visible',false);
        plot(U{1}(1:fs:end,i),'LineWidth',2)
        xlabel('Time (sec)','Fontsize',18,'FontName','Times New Roman')
        ylabel('Magnitude','Fontsize',18,'FontName','Times New Roman')
        xticks(0:noSamples_time/fs)
        xticklabels(0:noSamples_time/fs)
        set(gca,'box','off', 'FontSize', 16);
        
        F    = getframe(FigT);
        if left
            imwrite(F.cdata, fullfile(saveDir,'figures',...
                strcat('left_Temporal_mean',num2str(i),'.png')), 'png')
            close 
        else
            imwrite(F.cdata, fullfile(saveDir,'figures',...
                strcat('right_Temporal_mean',num2str(i),'.png')), 'png')
            close 
        end
    end
        
 %plotting spatial components 
  for i=1:r_components(2)   
      if i ==1 || rem(i-1,4) == 0
            count = 1;
            FigC = figure('Position', get(0, 'Screensize'),'visible',false);
      else
            count = count +1;
      end
        if left
            subplot(1,4,count);
            data=U{2}(:,i);
            
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
            
        else
            subplot(1,5,count);
            data=U{2}(:,i);
            
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
        end
        if rem(i,4) == 0|| i == r_components(2)
            F    = getframe(FigC);
            if left
                imwrite(F.cdata, fullfile(saveDir,'figures',...
                    strcat('left_Spatial_mean',num2str(i),'.png')), 'png')
                close 
            else
                imwrite(F.cdata, fullfile(saveDir,'figures',...
                    strcat('right_Spatial_mean',num2str(i),'.png')), 'png')
                close 
            end
        end
  end
  
  %plotting subject components
  for i=1:r_components(3)
      if i ==1 || rem(i-1,5) == 0
            count = 1;
            FigS = figure('Position', get(0, 'Screensize'),'visible',false);
      else
            count = count +1;
      end
        subplot(1,5,count);
        hold on;
        starti=1;
        for groupI=1:length(groups)
            endi=starti+length(find(group_ind==groupI))-1;
            bar(groupI,mean(U{3}(starti:endi,i)));
            starti=endi+1;
        end
        xlabel('Subject #','Fontsize',18,'FontName','Times New Roman')
        ylabel('Signature','Fontsize',18,'FontName','Times New Roman')
        set(gca,'box','off', 'FontSize', 16);
        legend(groups,'Location','best');
        if rem(i,5) == 0|| i == r_components(3)
            F    = getframe(FigS);
            if left
                imwrite(F.cdata, fullfile(saveDir,'figures',...
                    strcat('left_Subject_mean',num2str(i),'.png')), 'png')
                close 
            else
                imwrite(F.cdata, fullfile(saveDir,'figures',...
                    strcat('right_subject_mean',num2str(i),'.png')), 'png')
                close 
            end
        end
        
  end 

% plotting individual subject components 
  for i=1:r_components(3)
      
        FigS = figure('Position', get(0, 'Screensize'),'visible',false);
      
        subplot(1,1,1);
        hold on;
        starti=1;
        for groupI=1:length(groups)
            endi=starti+length(find(group_ind==groupI))-1;
            bar(groupI,mean(U{3}(starti:endi,i)));
            starti=endi+1;
        end
        xlabel('Subject #','Fontsize',18,'FontName','Times New Roman')
        ylabel('Signature','Fontsize',18,'FontName','Times New Roman')
        set(gca,'box','off', 'FontSize', 16);
        legend(groups,'Location','best');
        
        F    = getframe(FigS);
        
        if left
            imwrite(F.cdata, fullfile(saveDir,'figures',...
                strcat('left_subject_indiv_mean',num2str(i),'.png')), 'png')
            close 
        else
            imwrite(F.cdata, fullfile(saveDir,'figures',...
                strcat('right_subject_indiv_mean',num2str(i),'.png')), 'png')
            close 
        end          
  end
  
end