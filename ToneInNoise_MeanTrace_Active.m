function ToneInNoise_MeanTrace_Active(DFF,ID,Type)

    mu = nanmean(DFF(:,:,:),2);
    std = nanstd(DFF(:,:,:),[],2)/sqrt(size(DFF,2));
    hold on   
    shadedErrorBar([],mu,std,'k')
    %plot(mean(Vec_DFF(:,:,nn),2),'k','LineWidth',2)
    aa = axis;
    plot([aa(1), aa(2)], [0 0 ] ,'--k')
    plot([30 30], [aa(3),aa(4)],'--g')
    plot([60 60], [aa(3),aa(4)],'--g')
    %plot([90 90], [aa(3),aa(4)],'--g')
    %plot([120 120], [aa(3),aa(4)],'--r')
    xlim([0 120])
    if exist('ID', 'var')
        if isnumeric(ID)
            title(['Neuron' num2str(ID)])
        else
            title(ID)
        end
    end
  end 