function CPD_plots_mean(U,r_components,left,condition_ind,conditions,saveDir,tick_time,tick_freq)
    %Generate a figure for each CPD component contains subplots of temporal (time),
    %spectral (frequency), spatial (channels), and subjects. Plotting channels 
    %activations on the left and right hemisphere. Taking condition indeces
    %to use in bar plots.
    %Plotting the mean of each condition of subjects

    noSamples_time=length(U{1});
    noSamples_freq=length(U{2});
    xbase=1;ybase=1; activations=zeros(4+xbase,8+ybase);
    Infimg_left=imread('.\left pic.gif',1);
    Infimg_left=imresize(Infimg_left,0.2); 
    Infimg_right=imread('.\right pic.gif',1);
    Infimg_right=imresize(Infimg_right,0.5); 
    if exist(fullfile(saveDir,'figures'))~=7
        mkdir(fullfile(saveDir,'figures'));
    end
    for i=1:r_components
        FigH = figure('Position', get(0, 'Screensize'),'visible',false);

        % plotting temporal
        subplot(1,4,1);
        plot(U{1}(:,i),'LineWidth',2)
        xlabel('Time (sec)','Fontsize',18,'FontName','Times New Roman')
        ylabel('Magnitude','Fontsize',18,'FontName','Times New Roman')
        xticks(0:1/tick_time:noSamples_time)
        xticklabels(0:noSamples_time*tick_time)
        set(gca,'box','off', 'FontSize', 16);
        
        % plotting spectral
        subplot(1,4,2);
        plot(U{2}(:,i),'LineWidth',2)
        xlabel('Frequency (Hz)','Fontsize',18,'FontName','Times New Roman')
        ylabel('Magnitude','Fontsize',18,'FontName','Times New Roman')
        xticks(0:tick_freq:noSamples_freq)
        xticklabels(0:noSamples_freq/tick_freq)
        set(gca,'box','off', 'FontSize', 16);
        
        % plotting spatial 
        if left
            subplot(1,4,3);
            data=U{3}(:,i);
            
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
            %activations in each CPD component. Not between all components
            %because they are not related.
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
            subplot(1,4,3);
            data=U{3}(:,i);
            
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
            %activations in each CPD component. Not between all components
            %because they are not related.
            levels=linspace(min(min(activations)),max(max(activations)),10); 
            contour(x,y,activations,levels,'ShowText','off','LineWidth',1.5) %narrow lines means negative levels. 
            contourcmap('hsv', levels,'Colorbar', 'on', ...
           'Location', 'horizontal', ...
           'TitleString', strcat('The 10 levels of magnitude are between ',num2str(min(min(activations)))...
                ,' to ',num2str(max(max(activations)))));
            axis equal                                % set axis units to be the same size
            %box on                                    % display bounding box
            axis square off
            imgh=imshow(repmat(Infimg_right,1,1,3), 'XData', [-4 21], 'YData', [-9 15]);
            imgh.AlphaData = .2;
            hold off
            %zoom(2)
            xlim([0.77 9.716])
            ylim([0.8543 7.795]) 

        end

        % plotting mean subject within a condition 
        subplot(1,4,4);
        hold on;
        starti=1;
        for conditionI=1:length(conditions)
            endi=starti+length(find(condition_ind==conditionI))-1;
            bar(conditionI,mean(U{4}(starti:endi,i)));
            starti=endi+1;
        end
        xlabel('Subject #','Fontsize',18,'FontName','Times New Roman')
        ylabel('Signature','Fontsize',18,'FontName','Times New Roman')
        set(gca,'box','off', 'FontSize', 16);
        legend(conditions,'Location','best');
        
        F    = getframe(FigH);
        if left
            imwrite(F.cdata, fullfile(saveDir,'figures',...
                strcat('left_mean',num2str(i),'.png')), 'png')
            close 
        else
            imwrite(F.cdata, fullfile(saveDir,'figures',...
                strcat('right_mean',num2str(i),'.png')), 'png')
            close 
        end
    end
end