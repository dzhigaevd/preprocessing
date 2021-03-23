%% Binning

binning_size = 4;
data_unbinned = data;

for ii = 1:size(data,3)
    for jj = 1:size(data,4)
        for kk = 1:size(data,5)
            convoluted = conv2(data(:,:,ii,jj,kk), ones(binning_size));
            convoluted_size = size(convoluted);
            data_binned(:,:,ii,jj,kk) = convoluted(binning_size:binning_size:convoluted_size(1), binning_size:binning_size:convoluted_size(2));
        end
    end
end  

data = data_binned;
clear data_binned;

%% Q-space calculation

data = double(data);

% Experiment parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NanoMax convention:
% gamma - horizontal detector
% delta - vertical detector
% gonphi - rotation about vertical axis
% gontheta - rotation about horizontal axis
nanomax.photon_energy   = 15000;
nanomax.gonphi      = [5 5.02 5.04 5.06 5.08 5.10 5.12 5.14 5.16 5.18 5.20 5.22 5.24 5.26 5.28...
                   5.30 5.32 5.34 5.36 5.38 5.40 5.42 5.44 5.46 5.48 5.50 5.52 5.54 5.56...
                   5.58 5.60]; % [deg] - can be a range of angles
nanomax.gontheta        = 0; % [deg] - can be a range of angles

nanomax.radius          = 1; % [m]
nanomax.delta           = 2.1; % [deg] these angles are corrected with the sign respecting the rotation rules
nanomax.gamma           = 13; % [deg] 
if exist('binning_size') ~= 0  
    nanomax.detector_pitch  = 55e-6*binning_size; % [m]
else
    nanomax.detector_pitch  = 55e-6;
end

nanomax.direct_beam     = [200,131];

nanomax.dPhi            = 0.02; % degrees

% Constants. They are needed for correct labeling of axes
h                       = 4.1357e-15;                                  % Plank's constant
c                       = 2.99792458e8;                                % Speed of light in vacuum

wavelength = h*c/nanomax.photon_energy;

k = 2*pi/wavelength; % wave vector

[hd,vd] = meshgrid(-size(data,2)*nanomax.detector_pitch/2:nanomax.detector_pitch:(size(data,2)/2-1)*nanomax.detector_pitch,...
    -size(data,1)*nanomax.detector_pitch/2:nanomax.detector_pitch:(size(data,1)/2-1)*nanomax.detector_pitch);

hd = hd+(size(data,2)/2-nanomax.direct_beam(2)).*nanomax.detector_pitch;
vd = vd+(size(data,1)/2-nanomax.direct_beam(1)).*nanomax.detector_pitch;
zd = ones(size(vd)).*nanomax.radius;

d = [hd(:),vd(:),zd(:)]';

r = squeeze(sqrt(sum(d.^2,1)));

hq = k*d(1,:)./r;
vq = k*d(2,:)./r;
zq = k*(1-d(3,:)./r);

q = [hq;vq;zq];

% Check the coordinate mapping - real space m, direct scattering
test = data(:,:,1,1,1);
figure;
scatter3(q(1,:),q(2,:),q(3,:),10,test(:),'fill','s');
xlabel('qh');
ylabel('qv');
zlabel('qz');
view(7.7549,53.9161);
% PROVED!
%

% Sample orientation matrix. Bounds the sample crystal with the laboratory frame
% Angles alpha beta gamma were manually adjusted so that known peaks 
% are exactly in their places

% X is horizontal, perp to the beam, Y is vertical

Rh = [1         0                    0; % detector rotation around horizintal axis 
      0         cosd(nanomax.delta) -sind(nanomax.delta);
      0         sind(nanomax.delta)  cosd(nanomax.delta)]; 

Rv = [cosd(nanomax.gamma)  0  sind(nanomax.gamma); % detector rotation around vertical axis 
      0                    1  0;
      -sind(nanomax.gamma) 0  cosd(nanomax.gamma)];

Rz = [cosd(0) -sind(0) 0; % detector rotation around beam axis 
      sind(0)  cosd(0)  0;
      0        0         1];

U = Rh*Rv*Rz;
    
qR = (U*q); % correct so far in real space

% Check the coordinate mapping - real space m, Bragg condition
figure;
scatter3(qR(1,:),qR(2,:),qR(3,:),10,test(:),'fill','s');
xlabel('qh');
ylabel('qv');
zlabel('qz');
view(7.7549,53.9161);
% PROVED!

%
% Initial coordinate of ki
ki = [0,0,k]';

kf = U*ki;

Q = kf-ki;

% Lab coordinate system: accosiated with the ki
QLab(1,:) = (qR(1,:)+Q(1));
QLab(2,:) = (qR(2,:)+Q(2));
QLab(3,:) = (qR(3,:)+Q(3));

scan.data_meta.QLab = QLab;

% Check the lab space mapping
% figure;
% scatter3(QLab(1,:),QLab(2,:),QLab(3,:),10,scan.data_average(:),'fill','s');
% xlabel('Qx');
% ylabel('Qy');
% zlabel('Qz');
% view(7.7549,53.9161);
% PROVED!

%
% Small corrections to misalignment of the sample
% Here the rocking curve should be introduced
% alpha
% beta

% Gonphi correction
alpha = 0;
beta = 0;
for ii = 1:length(nanomax.gonphi) %or gontheta
    % Rotations to bring the q vector into sample coordinate system
    Rsh = [1         0                       0; % detector rotation around horizintal axis 
           0         cosd(nanomax.gontheta) -sind(nanomax.gontheta);
           0         sind(nanomax.gontheta)  cosd(nanomax.gontheta)]; 

    Rsv = [cosd(-nanomax.gonphi(ii)+alpha)  0  sind(-nanomax.gonphi(ii)+alpha); % detector rotation around vertical axis 
           0                      1  0;
          -sind(-nanomax.gonphi(ii)+alpha)  0  cosd(-nanomax.gonphi(ii)+alpha)];

    Rsz = [cosd(0+beta) -sind(0+beta) 0; % detector rotation around beam axis 
           sind(0+beta)  cosd(0+beta)  0;
           0        0         1];

    Rs = Rsh*Rsv*Rsz; 

    % Sample coordinate system: accosiated with the ki
    QSample(:,:,ii) = Rs*QLab;
%     Qs = Rs*Q;
% 
%     modQSample = squeeze(sqrt(sum(Qs.^2,1))); %change alpha, beta until =H111
end

% scan.data_meta.QSample = QSample;

% Check the sample space mapping
figure;
for ii = 1:length(nanomax.gonphi)
    Q1 = squeeze(QSample(1,:,ii));
    Q2 = squeeze(QSample(2,:,ii));
    Q3 = squeeze(QSample(3,:,ii));
    t  = squeeze(sum(data(:,:,:,:,ii),[3,4]));
    
    scatter3(Q1(:),Q2(:),Q3(:),10,t(:),'fill','s');
    xlabel('Qx');ylabel('Qy');zlabel('Qz');
%     view(7.7549,53.9161);
    hold on;
end
view(0,0);
    view(0,90);
% PROVED!


% interpolate Q
scaleCoefficient = 1;

qx = squeeze(QSample(1,:,:));
qy = squeeze(QSample(2,:,:));
qz = squeeze(QSample(3,:,:));

dqX = scaleCoefficient*(2*pi*nanomax.detector_pitch/(nanomax.radius*wavelength));
dqY = scaleCoefficient*(2*pi*nanomax.detector_pitch/(nanomax.radius*wavelength));
dqZ = scaleCoefficient*(2*pi*nanomax.detector_pitch/(nanomax.radius*wavelength));

[Qx,Qy,Qz] = meshgrid(min(qx(:)):dqX:max(qx(:)), min(qy(:)):dqY:max(qy(:)),min(qz(:)):dqZ:max(qz(:)));
mkdir('Documents/data/analysis/qspace_full');

%% ... and interpolation (which takes forever)
for ii = 1:size(data,4)
    tic;
    parfor (jj = 1:size(data,3),30)  
        t = squeeze(data(:,:,jj,ii,:));          
        F = TriScatteredInterp(qx(:),qy(:),qz(:),t(:));       
        lab_data = F(Qx,Qy,Qz);
        lab_data(isnan(lab_data))=0;
%                 denom = sum(lab_data(:));
%                 COM_qx(jj,ii) = sum(sum(sum(lab_data.*Qx)))/denom;
%                 COM_qz(jj,ii) = sum(sum(sum(lab_data.*Qz)))/denom;
%                 COM_qy(jj,ii) = sum(sum(sum(lab_data.*Qy)))/denom;

        q_data(:,:,jj,ii,:) = lab_data;
            %sprintf('Done [%d %d] /n',jj,ii);
            %com_data(:,jj,ii) = ndimCOM(lab_data,'auto');
    end    
    toc
end
   

save('Documents/data/analysis/qspace_full/q_data.mat','q_data','-v7.3');
% save('Documents/data/analysis/qspace_full/com_qx.mat','COM_qx','-v7.3');
% save('Documents/data/analysis/qspace_full/com_qz.mat','COM_qz','-v7.3');
% save('Documents/data/analysis/qspace_full/com_qy.mat','COM_qy','-v7.3');

%%  Visualisation (3D Bragg peak)
clc;

v(:,:,:) = q_data(:,:,1,1,:);

p = patch(isosurface(Qx,Qy,Qz,v));
hold on;

m = data(:,:,1,1,1);
%scatter3(QSample(1,:,15),QSample(2,:,15),QSample(3,:,15),10,m(:),'fill','s');

% Plot settings
p.FaceColor = 'blue';
p.EdgeColor = 'none';
p.FaceVertexAlphaData = 0.5;
p.FaceAlpha = 'flat'; 
xlabel('Qx')
ylabel('Qy')
zlabel('Qz')
daspect([1 1 1])
%view(0,90); % Along beam direction
view(3);
axis tight
camlight 
lighting gouraud

%% Strain and tilt

clc;

q_data_average = squeeze(mean(q_data,[3 4])); 
v = q_data_average;

COM_average = centerofmass(Qx,Qy,Qz,v);

figure();
p = patch(isosurface(Qx,Qy,Qz,v));
hold on;
plot3(out(1),out(2),out(3),'g*');

% Plot settings
p.FaceColor = 'blue';
p.EdgeColor = 'none';
p.FaceVertexAlphaData = 0.5;
p.FaceAlpha = 'flat'; 
xlabel('Qx')
ylabel('Qy')
zlabel('Qz')
daspect([1 1 1])
%view(0,90); % Along beam direction
view(3);
axis tight
camlight 
lighting gouraud

for i=1:size(q_data,3)
    for j=1:size(q_data,4)
        v(:,:,:) = q_data(:,:,i,j,:);
        COM = centerofmass(Qx,Qy,Qz,v);
        COMx = COM(1);
        COMy = COM(2);
        COMz = COM(3);
        
        q_mod_average = sqrt(COM_average(1)^2+COM_average(2)^2+COM_average(3)^2);
        q_mod = sqrt(COMx.^2+COMy.^2+COMz.^2);

        q_tilt_z(i,j) = asind(COMy./q_mod)-atand(COM_average(2)/COM_average(1));
        q_tilt_y(i,j) = asind((COMz)./q_mod)-atand(COM_average(3)/COM_average(1));
        q_strain(i,j) = q_mod_average./q_mod-1;
    end
end

for ii = 1:size(q_data,4)
    for jj = 1:size(q_data,3)
        map(ii,jj)= sum(sum(sum(squeeze(q_data(:,:,jj,ii,:)))));
    end
end
%%
map_threshold = 1.2e+06;
map_show = map>map_threshold;
%map_show = ones(10,13);

figure(); 
im0 = subplot(1,4,1); imagesc(map'.*map_show');
cmap = jet(256); cmap(1,:) = [0,0,0]; colormap(im0, cmap); colorbar; 
axis image; title('Bragg intensity');

im1 = subplot(1,4,2); imagesc(q_tilt_z.*map_show'); 
cmap = jet(256); cmap(179,:) = [0,0,0]; colormap(im1, cmap); colorbar; 
axis image; title('Tilt (y), [deg]');

im2 = subplot(1,4,3); imagesc(q_tilt_y.*map_show'); 
cmap = jet(256); cmap(1,:) = [0,0,0]; colormap(im2, cmap); colorbar; 
axis image; title('Tilt (z), [deg]');

im3 = subplot(1,4,4); imagesc(100*q_strain.*map_show'); 
cmap = jet(256); cmap(62,:) = [0,0,0]; colormap(im3, cmap); colorbar; 
axis image; title('Strain, [%]');