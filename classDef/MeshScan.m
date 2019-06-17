classdef MeshScan < Scan
    %DIFFDATA Summary of this class goes here
    %   Detailed explanation goes here
        
    properties
        STXMmap; % intensity sum of each diffraction pattern as a function of ptychographic scan position
    end
        
    % Constructor of the object
    methods
        function obj = MeshScan(dataArray,fioParam)
            inputParser
            obj@Scan(dataArray,fioParam);
        end
    end
    
    % For calculations
    methods 
        function handles = showScanSTXMLive(obj,mode)     
            global KEY_IS_PRESSED
            KEY_IS_PRESSED = 0;
            
            if nargin == 1
                mode = 'rect';
            end

            handles.figHandle = figure;
            subplot(1,2,1); imagesc(log10(averageScan(obj))); axis image; title('Integrated intensity');

            map = zeros(obj.metaData.nV,obj.metaData.nH);

            if strcmp(mode,'rect')
                h = imrect;
            else
                h = impoly;
            end
            
            while 1   
                if strcmp(mode,'rect')
                    pos = round(getPosition(h)); %[xmin ymin width height]
                else
                    mask = createMask(h);
                end

                kk = 1;
                for ii = 1:obj.metaData.slowMotorPoints
                    for jj = 1:obj.metaData.fastMotorPoints
                        if strcmp(mode,'rect')
                            map(ii,jj) = sum(sum(obj.data(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3),kk)));
                        else
                            map(ii,jj) = sum(sum(obj.data(:,:,kk).*mask));
                        end
                        kk = kk+1;                        
                    end
                end
                subplot(1,2,2); imagesc(obj.metaData.xVector,obj.metaData.yVector,map); 
                xlabel(obj.metaData.fastMotorName);
                ylabel(obj.metaData.slowMotorName);
                title('PRESS CTRL+C TO STOP!');
                drawnow;
            end           
        end
        
%         function handles = liveIntegrate(obj,mode)
%             if nargin == 1
%                 mode = 'rect';
%             end
%             
%             handles.imHandle = figure;
%             subplot(1,2,1); imagesc(log10(sumPos(obj))); axis image; title('Integrated intensity');
%             
%             integral = zeros(1,obj.dataSize(3));
%             
%             if strcmp(mode,'rect')
%                 h = imrect;
%             else
%                 h = impoly;
%             end
%             
%             while 1   
%                 if strcmp(mode,'rect')
%                     pos = round(getPosition(h)); %[xmin ymin width height]
%                 else
%                     mask = createMask(h);
%                 end
% 
%                 for kk = 1:obj.dataSize(3)                    
%                     if strcmp(mode,'rect')
%                         integral(kk) = sum(sum(obj.data(pos(2):pos(2)+pos(4),pos(1):pos(1)+pos(3),kk)));
%                     else
%                         integral(kk) = sum(sum(obj.data(:,:,kk).*mask));
%                     end                                        
%                 end
%                 subplot(1,2,2); plot(integral);title('Integrated intensity within ROI');
%                 drawnow;
%             end
%         end
        
        function STXMmap = STXM(obj)
            % Full integrated map
            kk = 1;
            STXMmap = zeros(obj.metaData.nV,obj.metaData.nH);
            
            for ii = 1:obj.metaData.slowMotorPoints
                for jj = 1:obj.metaData.fastMotorPoints
                    STXMmap(ii,jj) = sum(sum(obj.data(:,:,kk)));
                    kk = kk+1;                        
                end
            end           
        end                
    end
        
    % For display
    methods                                                                                     
        function handles = showScanSTXM(obj) 
            handles.figHandle = figure;
            handles.imHandle = imagesc(obj.metaData.xVector,obj.metaData.yVector,STXM(obj));           
            xlabel(obj.metaData.fastMotorName);
            ylabel(obj.metaData.slowMotorName);            
            colorbar
        end
    end
    
end

