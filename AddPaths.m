 function AddPaths(dirOrg)       
    
    %*****************************
    % Set current data folder
    %*****************************
    global dirData    
    
    dirData = ['..' filesep '2DChebClassData'];
    %***********************************
    % Initialization of main data folder
    %***********************************
    
    global dirDataOrg    %Main data folder 
    dirDataOrg  = dirData; 
    if(nargin >= 1)
        ChangeDirData([dirDataOrg filesep dirOrg],'ORG');
    end
    
    %******************************
    % Add subfolders to matlab path
    %******************************
    addpath(genpath(pwd));                       
    
    %*********************************************
    % Initialization of recompute and load options
    %*********************************************
    
    global recomputeAll  %if true, all data is recomputed, even if 
                         %   precomputed data is found in dirData folder
    global loadAll       %if true, all data is loaded if available in 
                         %   dirData folder. This overrides recomputeAll
    %(see also DataStorage.m for use of recomputeAll / loadAll)
    
    if(isempty(recomputeAll))
       recomputeAll = false;
    end
    
%     if(isempty(loadAll))
%         no = fprintf('Do you want to be asked if data should be recomputed? (press any key for -yes-)\n');            
%         if(getkeywait(2) == -1)            
%             loadAll = true;
%         else
%             cprintf('*g','Thanks. You will be asked if data should be recomputed.\n');
%             loadAll = false;
%         end
%     end    

    loadAll = true;
    
    %******************
    % Display Info
    %******************
    if(recomputeAll)
        cprintf('*m','No precomputed data will be used.\n');
    elseif(loadAll)
        cprintf('*m','Precomputed data will be used if available.\n');
    end
    
end