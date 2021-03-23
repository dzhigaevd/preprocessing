classdef mainGUI < handle

    
    %class properties - access is private so nothing else can access these
    %variables. Useful in different sitionations
    properties (Access = private)                
        gui_h;
    end
% ( obj )

    methods
         %function - class constructor - creates and init's the gui
        function this = mainGUI
            
            %make the gui handle and store it locally
            this.gui_h = guihandles(standard_fig);
            
            %set the callback functions to the two edit text box's
            set(this.gui_h.density_box, 'callback', @(src, event) Edit_density(this, src, event));
            set(this.gui_h.volume_box , 'callback', @(src, event) Edit_volume(this, src, event));
            
            %set the callback functions to the two buttons (calculate &
            %reset)
            set(this.gui_h.calculate_btn, 'callback', @(src, event) Calculate_callback(this, src, event));
            set(this.gui_h.reset_btn,     'callback', @(src, event) Reset_callback   (this, src, event));
            
            %Set the selection change Fcn for the radio button box. This
            %function will be called when the selection changes within the
            %box
            set(this.gui_h.unitgroup, 'selectionchangefcn', @(src, event) Ui_callback      (this, src, event));
            
            %sets the figure close function. This lets the class know that
            %the figure wants to close and thus the class should cleanup in
            %memory as well
            set(this.gui_h.figure1,  'closerequestfcn', @(src,event) Close_fcn(this, src, event));
            
            %reset the gui (not needed, but this is used to duplicate the
            %functionality of the matlab example)
            this = Reset(this);
            
        end
    end
    mainF = figure('Name','Plot Control',...
        'NumberTitle','off','units',...
        'normalized','outerposition',[0 0.2 0.2 0.3]);
    
    % Button elements =====================================================
    % - Show single frame -
    h = uicontrol('Parent',mainF,'Style','edit','Units','normalized',...
        'Position',[0.55 0.75 0.2 0.1],'String','1');
    
    uicontrol('Parent',mainF,'Units','normalized',...
        'Position',[0.05 0.75 0.4 0.1],'String','Single Frame',...
        'Callback',@runShowSingle);                
    
    % - Show average pattern - 
    uicontrol('Parent',mainF,'Units','normalized',...
        'Position',[0.05 0.65 0.4 0.1],'String','Average Pattern',...
        'Callback',@runShowScanAverage);                
    
    % - Show integlal plot of the scan - 
    uicontrol('Parent',mainF,'Units','normalized',...
        'Position',[0.05 0.55 0.4 0.1],'String','Integral Plot',...
        'Callback',@runShowScanIntegral);
    
    % - Scroll through data - 
    h = uicontrol('Parent',mainF,'Style','edit','Units','normalized',...
        'Position',[0.55 0.45 0.2 0.1],'String','1');
    
    uicontrol('Parent',mainF,'Units','normalized',...
        'Position',[0.05 0.45 0.4 0.1],'String','Integral Plot',...
        'Callback',@runShowScrollData);
        
    % - Show STXM map - 
    uicontrol('Parent',mainF,'Units','normalized',...
        'Position',[0.05 0.35 0.4 0.1],'String','Integral Plot',...
        'Callback',@runShowScanSTXM);
    
     % - Show live STXM map - 
    uicontrol('Parent',mainF,'Units','normalized',...
        'Position',[0.05 0.25 0.4 0.1],'String','Integral Plot',...
        'Callback',@runLiveSTXM);
    
%     function runShowSingle(hObj,callbackdata)
%         index = str2num(get(h,'String'));
%         obj.showSingle(index);
%     end
% 
%     function runShowScanAverage(hObj,callbackdata)
%         obj.showScanAverage;
%     end
% 
%     function runShowScanIntegral(hObj,callbackdata)
%         obj.showScanIntegral;
%     end
% 
%     function runShowScrollData(hObj,callbackdata)
%         maxVal = str2num(get(h,'String'));
%         obj.showScrollData(maxVal);
%     end
%     
%     function runShowScanSTXM(hObj,callbackdata)
%         obj.showScanSTXM;
%     end
% 
%     function runLiveSTXM(hObj,callbackdata)
%         obj.liveSTXM;
%     end
% 

end

