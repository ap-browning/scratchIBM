function varargout = Plot_ScatterImage(varargin)
%SCATTERIMAGE
%I = Plot_ScatterImage(Data) creates an image with a white background, and red
%discs centred at locations specified by Data(:,1) and Data(:,2)
%I = Plot_ScatterImage(Experiment,Time) creates an image that
%superimposes the experimental image (cropped) with the cell locations (in yellow)

warning('off','imageio:tiffmexutils:libtiffWarning');

Width   = 1440;
Height  = 1900;

Radius  = 12;

IBMColor = [1,0,0];
DatColor = [1,1,0];


    switch nargin
        
        % Data is an input
        case 1
            
            Data = varargin{1};
            
            Dots = GetBinary(Data(:,1),Data(:,2));
            
            R = zeros(size(Dots)); G = R; B = R;
            R(Dots)     = IBMColor(1);
            G(Dots)     = IBMColor(2);
            B(Dots)     = IBMColor(3);
            R(~Dots)    = 1;
            G(~Dots)    = 1;
            B(~Dots)    = 1;
            
            I   = cat(3,R,G,B); 
            I   = imgaussfilt(I,1);
            
            if nargin == 2
                imwrite(I,varargin{2});
            end
            
        % Filenames are input
        case 2
            
            DensityRep  = varargin{1};
            Time        = varargin{2};
            
            % Check inputs
            if ~strcmp(Time(end),'h')
                Time = [num2str(Time),'h'];
            end
            
            Data    = csvread(['../Data/',DensityRep,'/PC3_',DensityRep,'_',Time,'.csv']);
            ExpI    = imread(['../Data/',DensityRep,'/PC3_',DensityRep,'_',Time,'.jpg']);            
            
            CurSz   = size(ExpI);
            
            ExpI    = imresize(ExpI,1440/CurSz(2));
            ExpI    = ExpI(1:1900,:,:);
            ExpI    = permute(ExpI,[2,1,3]);
            
            Dots    = GetBinary(Data(:,1),Data(:,2));
            
            R       = ExpI(:,:,1);
            G       = ExpI(:,:,2);
            B       = ExpI(:,:,3);
            
            R(Dots) = DatColor(1) * 255;
            G(Dots) = DatColor(2) * 255;
            B(Dots) = DatColor(3) * 255;
            
            I       = cat(3,R,G,B);
            I       = imgaussfilt(I,1);
                       
    end
    
    
    function Domain = GetBinary(X,Y)
        
        Domain = false(Width,Height);
        Mask   = false(2*Radius+1,2*Radius+1);
        
        for i = 1:2*Radius+1
           for j = 1:2*Radius+1

               if (i-Radius-1)^2 + (j-Radius-1)^2 <= Radius^2
                  Mask(i,j) = true; 
               end

           end
        end
        
        N = length(X);

        % Place circles
        for i = 1:N

            xp      = round(X(i));
            yp      = round(Y(i));

            if xp >= Radius+1 && xp <= (Width - Radius-1) 
                xlook = xp-Radius:xp+Radius;  
            elseif xp < Radius+1  
                xlook = [abs(xp-Radius+ Width):Width,1:xp+Radius];  
            elseif xp > Width-Radius-1  
                xlook = [xp-Radius:Width,1:mod(xp+Radius,Width)];
            end

            if yp >= Radius+1 && yp <= (Height - Radius-1) 
                ylook = yp-Radius:yp+Radius;  
            elseif yp < Radius+1  
                ylook = [abs(yp-Radius+Height):Height,1:yp+Radius];  
            elseif yp > Height-Radius-1  
                ylook = [yp-Radius:Height,1:mod(yp+Radius,Height)];
            end

            Domain(xlook,ylook) = Domain(xlook,ylook) + Mask;

        end
        
    end

    varargout{1} = I;
    if nargout == 2
        varargout{2} = length(Data);
    end
    
end