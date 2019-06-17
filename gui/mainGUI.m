function mainGUI( obj )
%MAINGUI Summary of this function goes here
%   Detailed explanation goes here
    % ============================ GUI ====================================
    % --- Create the starting window
    % =====================================================================  
        
    mainF = figure('Name','Plot Control',...
        'NumberTitle','off','units',...
        'normalized','outerposition',[0 0.2 0.2 0.3]);
    
    % - Show scan label for viewer -
    uicontrol('Parent',mainF,'Units','normalized','Style','text',...
        'Position',[0.05 0.9 0.8 0.07],'String',['Sample: ' obj.metaData.sampleName ' | Scan #' num2str(obj.metaData.scanNumber)]); 
    
    % Button elements + callback ==========================================
    
    % - Show single frame -
    hSingle = uicontrol('Parent',mainF,'Style','edit','Units','normalized',...
        'Position',[0.55 0.75 0.2 0.1],'String','1');
    
    uicontrol('Parent',mainF,'Units','normalized',...
        'Position',[0.05 0.75 0.4 0.1],'String','Single Frame',...
        'Callback',@run_showScanSingle);                
    
    function run_showScanSingle(hObj,callbackdata)
        index = str2num(get(hSingle,'String'));
        obj.showScanSingle(index);
    end

    % =====================================================================
  
    % - Show average pattern - 
    uicontrol('Parent',mainF,'Units','normalized',...
        'Position',[0.05 0.65 0.4 0.1],'String','Average Pattern',...
        'Callback',@run_showScanAverage);                
    
    function run_showScanAverage(hObj,callbackdata)
        obj.showScanAverage;
    end
    
    % =====================================================================
    
    % - Show integlal plot of the scan - 
    uicontrol('Parent',mainF,'Units','normalized',...
        'Position',[0.05 0.55 0.4 0.1],'String','Integral Scan',...
        'Callback',@run_showScanIntegral);
    
    function run_showScanIntegral(hObj,callbackdata)
        obj.showScanIntegral;
    end
    
    % =====================================================================
    
    % - Scroll through data - 
    hScroll = uicontrol('Parent',mainF,'Style','edit','Units','normalized',...
        'Position',[0.55 0.45 0.2 0.1],'String','1');
            
    uicontrol('Parent',mainF,'Units','normalized',...
        'Position',[0.05 0.45 0.4 0.1],'String','Scroll Data',...
        'Callback',@run_showScanScroll);
    
    function run_showScanScroll(hObj,callbackdata)
        maxVal = str2double(get(hScroll,'String'));
        obj.showScanScroll(maxVal);
    end
    
    % =====================================================================
    if strcmp(obj.metaData.scanType1,'mesh')
        % - Show STXM map - 
        uicontrol('Parent',mainF,'Units','normalized',...
            'Position',[0.05 0.35 0.4 0.1],'String','Scan STXM',...
            'Callback',@run_showScanSTXM);

        % =====================================================================

         % - Show live STXM map - 
        uicontrol('Parent',mainF,'Units','normalized',...
            'Position',[0.05 0.25 0.4 0.1],'String','Live STXM',...
            'Callback',@run_showScanSTXMLive);    
    end
    
    % =====================================================================
    if strcmp(obj.metaData.scanType,'d2scan')
        % - Show reflectivity curve - 
        uicontrol('Parent',mainF,'Units','normalized',...
            'Position',[0.05 0.25 0.4 0.1],'String','Reflectivity curve',...
            'Callback',@run_showReflectivity);  
    end
    
    function run_showScanSTXM(hObj,callbackdata)
        obj.showScanSTXM;
    end
    
    function run_showReflectivity(hObj,callbackdata)
        obj.showReflectivity;
    end

    function run_showScanSTXMLive(hObj,callbackdata)
        obj.showScanSTXMLive;
    end

    % =====================================================================
    % - Save gif - 
    uicontrol('Parent',mainF,'Units','normalized',...
        'Position',[0.05 0.15 0.4 0.1],'String','Save GIF',...
        'Callback',@run_saveGif);
    
    % - Save data - 
    uicontrol('Parent',mainF,'Units','normalized',...
        'Position',[0.05 0.05 0.4 0.1],'String','Save Data',...
        'Callback',@run_saveScan);
    
    function run_saveGif(hObj,callbackdata)
        obj.saveGif;
    end

    function run_saveScan(hObj,callbackdata)
        obj.saveScan;
    end

end

