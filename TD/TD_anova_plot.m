function TD_anova_plot(timept_only_sig_component,...
    chan_only_sig_component,left,cpd_dir,dir,name,namePNG)

load(strcat(cpd_dir,'/LMLRA_all'));%, 'U_left','U_right','condition_ind',...
    %'tensor_right', 'tensor_left','Sol_left','Sol_right','size_core_left',...
    %'size_core_right','conditions_title','conditionI');
noSamples_time=size(U_left{1},1);
%no_freq_samp=size(U_left{2},1);
fs=50;
window_size=3*fs; overlap=window_size*.9; %the length of each recording about 6.7 seconds
signal_length=1351; nfft=2^12;
time_sp=0:(window_size-overlap)/fs:(signal_length-window_size)/fs;
freq_sp=0:fs/nfft:1;%only the frequencies between 0 and 1

tick_time=(window_size-overlap)/fs; tick_freq=nfft/fs; 
conditions_title={'Human hand, function', 'Human hand, non-function', 'Mechanical hand, function', ...
    'Mechanical hand, non-function'};


xbase=1;ybase=1; activations=zeros(4+xbase,8+ybase);
Infimg_left=imread('left pic.gif',1);
Infimg_left=imresize(Infimg_left,0.2);%(Infimg,0.09);
Infimg_right=imread('right pic.gif',1);
Infimg_right=imresize(Infimg_right,0.5);%(Infimg,0.09);

% plotting temporal profile 
 FigH = figure('Position', get(0, 'Screensize'),'visible',false);
            plot(timept_only_sig_component,'LineWidth',8,'color','black')
            title(strcat(name), 'FontSize', 24)
            xlabel('Time (sec)','Fontsize',18,'FontName','Times New Roman')
            ylabel('Magnitude','Fontsize',18,'FontName','Times New Roman')
            xticks(0:fs*2:noSamples_time)
            xticklabels(-2:2:24)

            set(gca,'box','off', 'FontSize', 20);
             F    = getframe(FigH);
       if left
                imwrite(F.cdata, fullfile(dir,strcat('left_time_selectedComponent_norm_',namePNG,'.png')), 'png')
                saveas(FigH,fullfile(dir,strcat('left_time_selectedComponent_norm_',namePNG,'.fig')), 'fig')
                 close 
       else
                imwrite(F.cdata, fullfile(dir,strcat('right_time_selectedComponent_norm_',namePNG,'.png')), 'png')
                saveas(FigH,fullfile(dir,strcat('right_time_selectedComponent_norm_',namePNG,'.fig')), 'fig')
                 close 
       end 


%plotting ROI
if left
  FigH = figure('Position', get(0, 'Screensize'),'visible',false);
           
            data=chan_only_sig_component; 
    
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
            contour(x,y,activations,levels,'ShowText','off','LineWidth',7) %narrow lines means negative levels. 
            contourcmap('hsv', levels,'Colorbar', 'on', ...
           'Location', 'horizontal', ...
           'TitleString', strcat('The 10 levels of magnitude are between ',num2str(min(min(data)))...
                ,' to ',num2str(max(max(data)))));
            axis equal                                % set axis units to be the same size
            axis square off
            imgh=imshow(repmat(Infimg_left,1,1,3), 'XData', [-11 14], 'YData', [-9 15]);
            imgh.AlphaData = .4;
            hold off
            xlim([0.77 9.716])
            ylim([0.8543 7.795])
            title(strcat('Spatial'));

else
        FigH = figure('Position', get(0, 'Screensize'),'visible',false);
            data=chan_only_sig_component;

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
            contour(x,y,activations,levels,'ShowText','off','LineWidth',7) %narrow lines means negative levels. 
            contourcmap('hsv', levels,'Colorbar', 'on', ...
           'Location', 'horizontal', ...
           'TitleString', strcat('The 10 levels of magnitude are between ',num2str(min(min(activations)))...
                ,' to ',num2str(max(max(activations)))));
            axis equal                                % set axis units to be the same size
            axis square off
            imgh=imshow(repmat(Infimg_right,1,1,3), 'XData', [-4 21], 'YData', [-9 15]);
            imgh.AlphaData = .4;
            hold off
            xlim([0.77 9.716])
            ylim([0.8543 7.795]) 
            title(strcat('Spatial Component '));  
end
F    = getframe(FigH);
if left
        imwrite(F.cdata, fullfile(dir,strcat('left_channels_selectedComponent_norm_',namePNG,'.png')), 'png')
         saveas(FigH,fullfile(dir,strcat('left_channels_selectedComponent_norm_',namePNG,'.fig')), 'fig')
         close 
else
        imwrite(F.cdata, fullfile(dir,strcat('right_channels_selectedComponent_norm_',namePNG,'.png')), 'png')
        saveas(FigH,fullfile(dir,strcat('right_channels_selectedComponent_norm_',namePNG,'.fig')), 'fig')
         close 
end 
end



    
