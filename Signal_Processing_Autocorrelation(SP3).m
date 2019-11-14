% 
load('iq_PX-1000_20130520_200702_e02.60.mat')

% Simple program to calculate power and plot for PAR data 
% *** it is assumed that you load the data before running *** 
% 
% Calculate power using the mean of absolute value squared 
% Position indicator (in km)... convert polar (az and range) to rectangular 

%For dual-pol radars, select which polarization to view here
if ~exist('X')
        fprintf('Dual Pol system, changing X_h/X_v to X...\n');
        X = X_h;
        %Y = X_v;
end


if ~exist('Y')
        fprintf('Dual Pol system, changing X_h/X_v to X...\n');
        %X = X_h;
        Y = X_v;
end

%The PAR data is in int16 format to save space, this if statement
% will convert the int16 to double format. This is an important
% step as basic matlab commands cannot be performed on int16 data
if ~isfloat(X)
	fprintf('Converting data into floating point numbers ...\n');
	X = double(X);
end


if ~isfloat(Y)
	fprintf('Converting data into floating point numbers ...\n');
	Y = double(Y);
end

%There has been a little change on the name of the radar OU-PRIME
if strcmp(radar,'OUPRIME'), radar = 'OU-PRIME'; end

el_rad = el/180*pi; 
[r,az_rad] = meshgrid(((0:num_gates-1)*delr+r_min)/1e3,az_set/180*pi); 
x = r*cos(el_rad).*sin(az_rad); 
y = r*cos(el_rad).*cos(az_rad); 
z = r*sin(el_rad); 

R0_x=mean(X(:,:,1:num_pulses).*conj(X(:,:,1:num_pulses)),3);
R0_y=mean(Y(:,:,1:num_pulses).*conj(Y(:,:,1:num_pulses)),3);
% R1=mean(X(:,:,2:num_pulses).*conj(X(:,:,1:num_pulses-1)),3);
% R2=mean(X(:,:,3:num_pulses).*conj(X(:,:,1:num_pulses-2)),3);


%%%%%%%%%% Calculate The Power %%%%%%%%%%%%%%%%%%
pow_x=10*log10(R0_x);
pow_y = 10*log10(R0_y);

powdiff = pow_x - pow_y;


% plot([1:num_pulses]*pri,real(squeeze(X_h(1,6,:))), ...
%     [1:num_pulses]*pri,imag(squeeze(X_h(1,6,:))));
% xlabel('Time (sec)');
% ylabel('Q(t) I(t)');
% title('Ground Clutter: Range = 0.18 km; Azimuth Angle = 1 degree');
% legend('Q(t)', 'I(t)');
% grid;

% %Part 3

% I =  real(squeeze(X_h(53,451,:)));
% Q =  imag(squeeze(X_h(53,451,:)));
% 
% plot(I,Q,'.');
% xlabel('Amplitude (I(t))');
% ylabel('Amplitude (Q(t))');
% title('In-Phase and Quadrature Sequence (Tornado)');
% axis([-1000 1000 -1000 1000]);
% axis square;



%%%Actual concatonation Part 3

% 
% center_r = 441;
% center_az = 38;
% itot = real(squeeze(X_h(center_az,center_r,:)));
% qtot = imag(squeeze(X_h(center_az,center_r,:)));
% 
% for i = [-1,1]
%     for j = [-1,1]
%         itot = [itot ; real(squeeze(X_h(center_az+i,center_r+j,:)))];
%         qtot = [qtot ; imag(squeeze(X_h(center_az+i,center_r+j,:)))];
%     end
% end
% 
% figure    
% plot(itot,qtot,'.')
% % axis([-30 30 -30 30])
% title('I(t) Q(t) Concatonation Tornado')
% xlabel('Amplitude (I(t))');
% ylabel('Amplitude (Q(t))');
% grid;






% 
%  %point = real(squeeze(X_h(90,500,:)))+ imag(squeeze(X_h(90,500,:))) ;
% point = squeeze(X_h(38,441,:));
% autocorr_weather = (1/50)*abs(xcorr(point)/xcorr(point,0));
% %autocorr_weather = (1/50)*abs(xcorr(point));
% plot(autocorr_weather);
% xlabel('Lag');
% ylabel('Autocorrelation');
% title('Autocorrelation of Tornado, Range = 13.2 km, Azimuth = 38 degrees');

% %subplot(2,1,1);
% point_I = real(squeeze(X_h(90,15,:)));
% point_Q = imag(squeeze(X_h(90,15,:)));
% h= histogram(point_I, 10);
% xlabel('Lag');
% ylabel('PDF');
% title('I(t) PDF (Ground Clutter)')


% 
% %%%%%Problem 5
% point_I = (real(squeeze(X_h(90,500,:))));
% point_Q = imag(squeeze(X_h(90,500,:)));
% amplitude = sqrt(point_I.^2 + point_Q.^2);
% power = point_I.^2 + point_Q.^2;
% phase = atan(Q./I);
% histogram(power,10);
% title('Clear Air (Power) PDF');
% ylabel('PDF');


%%%%%Problem 6

% point = real(squeeze(X_h(38,441,:)))+ imag(squeeze(X_h(38,441,:)));
% autocorrelation = xcorr(point);
% f = fftshift(autocorrelation);
% plot(f);
% xlabel('Frequency');
% ylabel('Amplitude');


%%%%%Correct problem 6

% extract time series data

%point = squeeze(X_h(38,441,:));
% time_series = squeeze(X_h(359,500,:));
% 
% %ACF... note that matlab does not have the (1/M) term
% 
% [acf,lags]=xcorr(time_series);
% 
% %setup velocity axis for plotting S(f)
% % 
% npts = 128; % Number of points for zero-padding of fft
% va = lambda/4/pri;
% vel=-[-npts/2:(npts/2)-1]*2*va/npts;
% % 
% % %Spectral estimate
% % 
 d = rectwin(99);
% Sf=abs(fftshift(fft(acf.*d,npts)));
% plot(vel,Sf);
% xlabel('Velocity (m/s)');
% ylabel('Amplitude');
% title('PSD In Storm');


%%%%%Problem 7%%%%%%%%%%%
% 
%summed = zeros(360,1010); %Sum pulses and loop through each gate and azimuth angle

% 
% for i = 1:360
%     for j = 1:1010
%         summed(i,j) =  [mean(X_h(i,j,1:num_pulses).*conj(X_h(i,j,1:num_pulses)),3)];
%         
%     end
% end 
% 




acf_array = zeros(360,num_gates);
psd = zeros(360,num_gates);

for m = 1:360
    for n = 1:num_gates
        acf = xcorr(squeeze(X_h(m,n,:)), 'biased');
        %acf_array(m,n,:) = xcorr(X_h(m,n,1:num_pulses).*conj(X_h(i,j,1:num_pulses)));
        acf_array(m,n) = acf(50);
        psd(m,n) = abs(fftshift(fft(acf(50))));
    end
end

Now compute power

powh = 10*log10(acf_array);
powh2 = 10*log10(psd);


 figure; 
set(gcf,'render','painters');
if strcmp(radar,'OU-PRIME'), naz_max = 720; else naz_max=360; end
if length(az_set)<naz_max
	pcolor(x,y,powh2);
else
	pcolor([x;x(1,:)],[y;y(1,:)],[powh2;powh2(1,:)]);
end;
shading flat
axis equal
axis([-30 30 -30 30])
colormap(jet);
xlabel('\bf \fontsize{11} Zonal Distance (km)')
ylabel('\bf \fontsize{11} Meridional Distance (km)')
title('PSD Power');
colorbar;








% %%%%%%%%%%%% Plot the Data %%%%%%%%%%%%%%%%%%%%%%
% figure; 
% set(gcf,'render','painters');
% if strcmp(radar,'OU-PRIME'), naz_max = 720; else naz_max=360; end
% if length(az_set)<naz_max
% 	pcolor(x,y,pow_x);
% else
% 	pcolor([x;x(1,:)],[y;y(1,:)],[pow_x;pow_x(1,:)]);
% end;
% shading flat
% axis equal
% axis([-30 30 -30 30])
% colormap(jet);
% xlabel('\bf \fontsize{11} Zonal Distance (km)')
% ylabel('\bf \fontsize{11} Meridional Distance (km)')
% title(['\bf \fontsize{12}',datestr(scan_time),' El=',num2str(el,'%5.2f'),' degrees H-pol']);
% colorbar;%%Part 3 concatenated


% I =  real(squeeze(X_h(10:20,500,:)));
% Q =  imag(squeeze(X_h(10:20,500,:)));
% 
% plot(I,Q,'.');
% xlabel('Time (sec)');
% ylabel('I(t) Q(t)');
% title('I(t) Q(t) (Storm); 5 deg Azimuth Radius (centered at 15 km, 345 degrees ');
% axis([-1000 1000 -1000 1000]);
% axis square;

% 
% %%%%%%%%%%%%  Produce map layover %%%%%%%%%%%%%%%%%%
% bmapover(gca,[],radar,{'OK'});
% 
% clear X;
% 
% 
% 
