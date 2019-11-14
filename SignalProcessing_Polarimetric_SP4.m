% 
% Simple program to calculate power and plot for PAR data 
% *** it is assumed that you load the data before running *** 
% 
% Calculate power using the mean of absolute value squared 
% Position indicator (in km)... convert polar (az and range) to rectangular 

%For dual-pol radars, select which polarization to view here

load('iq_PX-1000_20130520_200702_e02.60.mat')

if ~exist('X')
        fprintf('Dual Pol system, changing X_h/X_v to X...\n');
        X = X_h;
        X = X_v;
end

%The PAR data is in int16 format to save space, this if statement
% will convert the int16 to double format. This is an important
% step as basic matlab commands cannot be performed on int16 data
if ~isfloat(X)
	fprintf('Converting data into floating point numbers ...\n');
	X = double(X);
end
%There has been a little change on the name of the radar OU-PRIME
if strcmp(radar,'OUPRIME'), radar = 'OU-PRIME'; end

el_rad = el/180*pi; 
[r,az_rad] = meshgrid(((0:num_gates-1)*delr+r_min)/1e3,az_set/180*pi); 
x = r*cos(el_rad).*sin(az_rad); 
y = r*cos(el_rad).*cos(az_rad); 
z = r*sin(el_rad); 

R0=mean(X(:,:,1:num_pulses).*conj(X(:,:,1:num_pulses)),3);
% R1=mean(X(:,:,2:num_pulses).*conj(X(:,:,1:num_pulses-1)),3);
% R2=mean(X(:,:,3:num_pulses).*conj(X(:,:,1:num_pulses-2)),3);


%%%%%%%%%% Calculate The Power %%%%%%%%%%%%%%%%%%
pow=10*log10(R0);

acf_array = zeros(360,num_gates);
psd = zeros(360,num_gates);
%velocity = zeros(360, num_gates);

for m = 1:360
    for n = 1:num_gates
        acf = xcorr(squeeze(X_h(m,n,1:5)), 'biased');
        %acf_array(m,n,:) = xcorr(X_h(m,n,1:num_pulses).*conj(X_h(i,j,1:num_pulses)));
        acf_array(m,n) = acf(6);
        %psd(m,n) = abs(fftshift(fft(acf(50))));
        velocity(m,n) = -1*(lambda/(4*pi*pri)).*angle(acf_array(m,n));
    end
end

for m = 1:360
    for n = 1:num_gates
        acf_v = xcorr(squeeze(X_v(m,n,:)), 'biased');
        acf_array_v(m,n) = acf_v(50);
       
    end
end

el_rad = el/180*pi;
[r,az_rad] = meshgrid(([0:num_gates-1]*delr+r_min)/1e3,az_set/180*pi);
r_cal_dB = 20*log10(1e3*r);
pow_dB = 10*log10(acf_array);
%dBZ = pow_dB;
dBZ = (pow_dB+r_cal_dB) -80;

pow_dB_v = 10*log10(acf_array_v);
dBZ_v = (pow_dB_v + r_cal_dB)-80;

%Compute Zdr
zdr = dBZ - dBZ_v;

%Compute rho_hv





%Now compute power

powh = 10*log10(acf_array);
powh2 = 10*log10(psd);


% %Rho_hv
% cc = zeros(360,num_gates);
% lag = zeros(360,num_gates);
% 
% 
% for i = 1:360
%     for j = 1:num_gates
%         cc = xcorr(squeeze(X_h(i,j,:)),squeeze(X_v(i,j,:)), 'biased');
%         lag(i,j) = cc(50);
%         
%     end
% end
% 
% rho_hv = abs(lag)./sqrt(acf_array.*acf_array_v);




%%%%%%%%%%%% Plot the Data %%%%%%%%%%%%%%%%%%%%%%
% figure; 
% set(gcf,'render','painters');
% if strcmp(radar,'OU-PRIME'), naz_max = 720; else naz_max=360; end
% if length(az_set)<naz_max
% 	pcolor(x,y,rho_hv);
% else
% 	pcolor([x;x(1,:)],[y;y(1,:)],[rho_hv;rho_hv(1,:)]);
% end;
% shading flat
% axis equal
% axis([-30 30 -30 30])
% %colormap(1)
% colormap = boonlib('carbmap',2)
% xlabel('\bf \fontsize{11} Zonal Distance (km)')
% ylabel('\bf \fontsize{11} Meridional Distance (km)')
% title(['\bf \fontsize{12}',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees (CC)']);
% caxis([0.7,1])
% colorbar;
% 
% %%%%%%%%%%%%  Produce map layover %%%%%%%%%%%%%%%%%%
% bmapover(gca,[],radar,{'OK'});
% 
% clear X;

%%Plotting for Q1
%%%%%%%%%%%% Plot the Data %%%%%%%%%%%%%%%%%%%%%%
figure; 
set(gcf,'render','painters');
if strcmp(radar,'OU-PRIME'), naz_max = 720; else naz_max=360; end
if length(az_set)<naz_max
	pcolor(x,y,dBZ);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[dBZ;dBZ(1,:)]);
end;
shading flat
axis equal
axis([-30 30 -30 30])
colormap(1)
%colormap = boonlib('rbmap',2)
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title(['\bf \fontsize{12}',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees (H-pol range corrected (5 pulses)']);
caxis([-10,60])
colorbar;

%%%%%%%%%%%%  Produce map layover %%%%%%%%%%%%%%%%%%
bmapover(gca,[],radar,{'OK'});

clear X;


