clc
clear

nxn=35;  % number of data points
fs=10^4;  % sampling frequency
%% --------------------------original signal-----------------------------------------
x=randi([0 8],1,32)+1i*randi([0 8],1,32);% generate signal
figure(1)
stem(abs(x))
title('input signal to the RF front end x[n]')
%% -----------------------------IFFT-------------------------------------------
xifft=ifft(x,32);% ifft
%% -----------------------------add CP------------------------------------------
xn=[xifft(29:32) xifft];% add CP
%% ------------------------------DAC------------------------------------------
t=0:1/fs:nxn+1-2/fs;
xp=zeros(1,length(t));
for i=0:nxn
    xp(fs*i+1)=xn(i+1);
end
figure(2)
plot(t,abs(xp))
title('impulse train x_p(t)')
%% ----------------------------window a signal-----------------------------------
xc=zeros(1,length(t));
for i=1:length(t)
    xc(i)=xn(fix(i/fs)+1);
end                        % Convert to square wave
figure(3)
plot(t,abs(xc))
title('after DAC x_c(t)')
%% -------------------------seperate real part and imaginary part-------------
xr=real(xc);
xi=imag(xc);
figure(4)
subplot(2,1,1)
plot(t,xr)
title('real part x_r')
subplot(2,1,2)
plot(t,xi)
title('imaginary part x_i')
%% ------------------------modulation----------------------------------------
xcos=xr.*cos(150.*t); % wc=150
xsin=xi.*sin(150.*t);
figure(5)
plot(t,xcos)
hold on
plot(t,xsin)
title('after modulation x_{cos} x_{sin}')
legend('x_{cos}','x_{sin}')
%% ---------------------sum-----------------------------------------------
xtr=xcos+xsin;
figure(6)
plot(t,xtr)
title('after summation x_{tr}')
%% block 4
h=zeros(1,length(t));
h(1)=0.5;
h(fs*1.5+1)=0.4;
h(fs*2.5+1)=0.3;
h(fs*3+1)=0.3;
figure(7)
plot(t,h)
title('channel impulse response h(t)')
%% ---------------------------receive------------------------------------
XTR=fft(xtr,length(t));
H=fft(h,length(t));
y=ifft(XTR.*H,length(t));
figure(8)
plot(t,y)
title('received signal after channel')

yn=ifft(fft(y,length(t))./H);
figure(9)
plot(t,yn)
title('divided by H(jw) time domain')
%% -----------------------demodulation----------------------------
yr=2*yn.*cos(150.*t);
yi=2*yn.*sin(150.*t);
figure(10)
subplot(2,1,1)
plot(t,yr)
hold on
plot(t,xr)
title('real part y_r before LPF')
subplot(2,1,2)
plot(t,yi)
hold on
plot(t,xi)
title('imaginary part y_i before LPF')
%% --------------------low pass filter----------------------------
[b,a]=butter(1,100/(fs/2));
yrf=filter(b,a,yr);
yif=filter(b,a,yi);
figure(11)
subplot(2,1,1)
plot(t,yrf)
hold on
plot(t,xr)
title('real part y_r after LPF')
subplot(2,1,2)
plot(t,yif)
hold on
plot(t,xi)
title('imaginary part y_i after LPF')
%% -------------------------sampling---------------------------------
t2=0:nxn;
yrd=zeros(1,nxn+1); % yrd means singal y, real part, discrete
for i=0:nxn
    for ii=fs*i+1:fs*i+fs-1
         yrd(i+1)=yrd(i+1)+yrf(ii);
     end
     yrd(i+1)=yrd(i+1)/fs;
end
yid=zeros(1,nxn+1); % yid means signal y, imaginary part, discrete
for i=0:nxn
    for ii=fs*i+1:fs*i+fs-1
         yid(i+1)=yid(i+1)+yif(ii);
     end
     yid(i+1)=yid(i+1)/fs;
end
figure(12)
subplot(2,1,1)
stem(t2,yrd,'*')
hold on
stem(t2,real(xn))
title('real part y_r after sampling(blue)')
subplot(2,1,2)
stem(t2,yid,'*')
hold on
stem(t2,imag(xn))
title('imaginary part y_i after sampling(blue)')
%% ---------------------------remove CP + FFT----------------------
yrr=[yrd(5:36)];% yrr: signal y, real part, CP removed
yir=[yid(5:36)];% yrr: signal y, imaginary part, CP removed
yfft=fft(yrr+1i*yir,32);
figure(13)
subplot(3,1,1)
stem(abs(x))
title('|x|')
subplot(3,1,2)
stem(abs(yfft))
title('|y|')
subplot(3,1,3)
stem(abs(x-yfft))
title('|x-y|')
