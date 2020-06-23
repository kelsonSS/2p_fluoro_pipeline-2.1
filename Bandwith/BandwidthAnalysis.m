function [BD,P,Levels] = BandwidthAnalysis(DF,Type,Classes,Lvl)

%  DF is a df_by_level (LevelX Freq X neuron) object and returns the bandwidth for each
%  neuron as a level X neuron object.

% types 
% BRFS - 'binary receptive field sum'
%
% RFS -'receptive field sum'
%
%

if ~exist('Lvl', 'var'); Lvl = .75; end

if ~ exist('Type','var'); Type = 'BRFS'; end 

if ~ exist('Classes','var'); Classes = 0; end 


Levels = [];

if ~isstruct(DF)
    errordlg('This function currently only accepts structs')
end 

    Freqs = unique(DF.FreqLevelOrder{:,1});
    Levels = sort(unique(DF.FreqLevelOrder{:,2}),'descend');
     
    active = DF.active{:,2};
   
    if iscell(DF.df_by_level)
        
        DF_shapes = cellfun(@size,DF.df_by_level,'UniformOutput',0);
        DF_shapes = cell2mat(cellfun(@(x) x(1:2)',DF_shapes,'UniformOutput',0));
        DF_shape = DF_shapes(:);
        DF_shape(DF_shapes == 0) = [];
        DF_shape = reshape(DF_shape,2,[])';
        rc = min(DF_shape);  % row column number for extraction
    end
    
    try
        df_by_level = DF.df_by_level;
        if  iscell(df_by_level)
            df_by_level = [];
            for ii = 1:length(DF.df_by_level);
                try
                df_by_level = cat(3,df_by_level,DF.df_by_level{ii}(1:rc(1),1:rc(2),:));
                catch
                end ;
            end
        end
    df_by_level = abs(df_by_level);    
        
        
        
        
    catch
        
        warning('Multiple FRA SHAPES DETECTED Aborting! \n')
        BD = []
        return
    end
  
    
    if Classes
      for Class = 1:length(DF.Classes)
        Class_idx = DF.Class_idx == Class;
        final_idx = (active>0) & Class_idx;
        BD{Class,1} = FindBandwidth(df_by_level(:,:,final_idx),Freqs,Lvl,Type,DF.Classes{Class});
        BD{Class,3} = DF.Classes{Class};
     end
     
    else
      DF = df_by_level(:,:,active>0);
      BD{1,1} = FindBandwidth(DF,Freqs,Lvl,Type);
      BD{1,3} = 'All'
    end 
  
   
    if  strmatch(Type, 'interp')
        for class_idx = 1:size(BD,1)
            BD2{class_idx,1} = cellfun(@(X) max(X(:,2)-X(:,1)) , BD{class_idx,1});
            BD3{class_idx,1} = cellfun(@(X) sum(X(:,2)-X(:,1)) , BD{class_idx,1});
            
            
        end
        BD2(:,3) = BD(:,3);
        BD3(:,3) = BD(:,3);
        
        
        
        PlotBandwidth(BD2);title('Max')
     BD = AnalyzeBandwidth(BD2);
        PlotBandwidth(BD3);title('Sum')
    else
    BD  = AnalyzeBandwidth(BD);
          PlotBandwidth(BD)
    end


end


function BD = AnalyzeBandwidth(BD)
     for ii =  1:size(BD,1)
        BD{ii,2} = findSignficance(BD{ii,1});
     end 


end 






function BD = FindBandwidth(DF,Freqs,Lvl,Type,ClassName)


if ~exist('ClassName','var')
    ClassName = 'All'
end 

m =  max(max(abs(DF)));
DF = DF./m;

lvl_idx = DF > Lvl;

switch Type
    case 'interp'
        Freqs = Freqs(1:size(DF,2));
        Freqs_lg2 = log2(Freqs);
        Freqs_lg2 = Freqs_lg2 - min(Freqs_lg2);
        interp_factor = .1; % interpolate by .X of an octave Ex .1 = 10th octave
        Freqs_lg2_interp = [min(Freqs_lg2):interp_factor:max(Freqs_lg2)];
        
        for neuron = 1:size(DF,3)
            V = DF(:,:,neuron);
            %interpolation
            DF_interp(:,:,neuron) = interp1(Freqs_lg2',V',Freqs_lg2_interp')';
            
            for lvl_idx = 1:size(DF_interp,1)
                % find_bandwidth
                inds= find(  DF_interp(lvl_idx,:,neuron) > Lvl );
                
                if isempty(inds)
                    bands(1,:) = [0 0];
                else
                % multipeaked neuron test
                % if there are multiple peaks it will show up as a beak in the
                % indicies
                % for example: inds = 1 2 3 10 11 12 20 21 22
                %              diff = 1 1 7 1  1  8  1  1
                %            brea ks = [3 6]
                %    bands{1} = Freqs_lg2_interp(inds([1,3])
                %    bands{2} = Freqs_lg2_interp(inds([4,6])
                %    bands{3} = Freqs_lg2_interp(inds([7,end])
                
                
                breaks = find(diff(inds) > 1)  ;
                
                if ~isempty(breaks)
                    
                    bands(1,:) = Freqs_lg2_interp(...
                        inds([1,breaks(1)]));
                    for ii = 1:length(breaks)
                        
                        if ii == length(breaks)
                            bands(ii+1,:) = Freqs_lg2_interp(...
                                inds(...
                                [breaks(ii)+1,end]) );
                        else
                            bands(ii+1,:) = Freqs_lg2_interp(...
                                inds(...
                                [breaks(ii)+1,breaks(ii+1)]) );
                        end
                    end
                    
                else
                    
                    bands(1,:) = Freqs_lg2_interp([inds(1),inds(end)] );
                end
            end 
                BD{lvl_idx,neuron} = bands;
                
                clear bands
                clear breaks
                clear inds
            end
            
        end
        
        
        
        
        
    case 'RFS'
       BD =  squeeze(sum( DF .* lvl_idx,2)); 
    case 'BRFS'
       BD = squeeze(sum(lvl_idx,2)); 
       figure
       imagesc(sum(lvl_idx,3))
       title( sprintf('%s Average FRA',ClassName),'Interpreter','none')
    otherwise
        errordlg('Illegal Type Keyword')
        return
end 
end



function PlotBandwidth(BD)
    figure
    hold on 
    if ~iscell(BD)
        errorbar([],nanmean(BD,2),nanstd(BD,[],2) / sqrt(size(BD,2)) * 1.96)
    end
    
    if iscell(BD{1,1})
        return
    end 
        
        
    for ii =  1:size(BD,1)
      
       errorbar([],nanmean(BD{ii,1},2),nanstd(BD{ii,1},[],2) / sqrt(size(BD{ii,1},2)))
       
    end 
     
     legend(BD(:,3),'Interpreter','none')
     xticks(0:1:4)
     xlabel('level')
     ylabel('Bandwidth (half-octaves) ')
    
    
    
end 




