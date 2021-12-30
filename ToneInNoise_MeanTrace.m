function ToneInNoise_MeanTrace(DFF,ID,Type)

 if ~exist('Type','var')
     Type = 'ToneNoisePassive'
 end 




    mu = nanmean(DFF(:,:,:),2);
    std = nanstd(DFF(:,:,:),[],2)/sqrt(size(DFF,2));
    hold on   
    shadedErrorBar([],mu,std,'k')
   
    aa = axis;
     plot([aa(1), aa(2)], [0 0 ] ,'--k')


    if exist('ID', 'var')
        if isnumeric(ID)
            title(['Neuron' num2str(ID)])
        else
            title(ID)
        end
    end
    
    % plot lines according to type
%     
%     switch Type
%         case 'ToneNoisePassive'
%             
%             plot([30 30], [aa(3),aa(4)],'--')
%             plot([60 60], [aa(3),aa(4)],'--g')
%             plot([90 90], [aa(3),aa(4)],'--g')
%             plot([120 120], [aa(3),aa(4)],'--r')
%             xlim([0 150])
%         case 'ToneNoiseActive' 
%             plot([30 30], [aa(3),aa(4)],'--g')
%             plot([60 60], [aa(3),aa(4)],'--g')
%             xlim([0 120])
%         
%             
%     end 
    
  end 