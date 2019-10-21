close all;
clear all;


%length of the signal
N=512;

%% produce signal using frequencies k (sinusoids) 
number_of_freq=2;
k=1+unidrnd(N-1,1,number_of_freq); % generate random frequncies from 1 to 500
amp=1+rand(1,number_of_freq)*9; % generate random amplitudes from 1 to 10

n=0:N-1;
  

x=exp(-(n-250).^2/10000).*sum(diag(amp)*sin(2*pi*k'*n/N)); % generate complex sinusoids vector
%x=sum(diag(amp)*sin(2*pi*k'*n/N))+5*rand(1,N); % generate complex sinusoids vector + noise

%% Sparse signal in frequency domain. 
xf=fftshift(fft(fftshift((x))));

%creating dft matrix
B=dftmtx(N);
Binv=inv(B);

%Taking DFT of the signal
xf=B*x';

%Selecting random rows\col of the DFT matrix
q=randperm(N); 


figure;
subplot(3,4,[1 4]);
plot(x);
hold on;
grid on;
xlabel('Samples');
ylabel('Amplitude');
xlim([0 N]);
title(['Original Signal ' num2str(N) ' samples with two different frequency sinsuoids']);

subplot(3,4,[5 6]);
plot(abs(xf));
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Original Signal, Discrete Fourier Transform');
xlim([0 N]);


%K=0;
%xp=rand(length(xf),1);
%while sum(abs(xf)./abs(xp))>2*number_of_freq %increase K until converged
    
    %Number of random observations to take
for K=10:5:256;
 %K=K+5;
    
    %creating measurement matrix
    A=Binv(q(1:K),:);
    
    %taking random time measurements
    y=(A*xf);
    
    %Calculating Initial guess
    x0=A'*y;
   
   %x0=rand(N,1);
  
    % Solve:  min_x ||x||_L1  s.t.  Ax = b
    xp=l1eq(x0,A,[],y,1e-5);
    
    %recovered signal in time domain
    xprec=real(Binv*xp);
    
    %Show me the money!
    subplot(3,4,[1 4]);
    hold on;
    plot(q(1:K),x(q(1:K)),'r.','MarkerSize',14);hold on;
    %grid on;
    %xlabel('Samples');
    %ylabel('Amplitude');
    xlim([0 N]);
    %title('Original Signal,1024 samples with two different frequency sinsuoids');
    
    subplot(3,4,[7 8])
    plot(abs(xp),'r');
    grid on;
    xlabel('Samples');
    ylabel('Amplitude');
    xlim([0 N]);
    title(sprintf('Recovered Signal, Discrete Fourier Transform sampled with %d samples',K));
    
    subplot(3,4,[9 12]);
    plot(xprec,'r')
    grid on;
    xlabel('Samples');
    ylabel('Amplitude');
    title(sprintf('Recovered Signal in Time Domain'));
    xlim([0 N]);
    getframe(gcf);
end
 
%%
return
figure;
% stem(q(1:K),x(q(1:K)),'r.','MarkerSize',20);hold on;
% plot(x0); hold on; plot(rand(500,1)+3); hold on; plot(rand(500,1)-3); hold on; plot(rand(500,1)+6);
subplot(1,2,1);
plot(x);hold on 
plot(q(1:K),x(q(1:K)),'r.','MarkerSize',14);hold on;
 

subplot(1,2,2);
stem(k,amp)
 %grid on;
    %xlabel('Samples');
    %
ylim([-12 12])

 figure;
 stem(real(y)) ;
 ylim([-4 8])