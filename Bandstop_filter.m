M = 45; As = 60; n=0:1:M-1;
beta = 0.1102*(As-8.7)+0.3;
w_kai = (kaiser(M,beta))'; wc1 = pi/3; wc2 = 2*pi/3;
hd = ideal_lp(wc1,M) + ideal_lp(pi,M) - ideal_lp(wc2,M);
h = hd .* w_kai; [db,mag,pha,grd,w] = freqz_m(h,1);
subplot(2,2,1); stem(n,hd); title('Ideal Impulse Response')
axis([-1 M -0.2 0.8]); xlabel('n'); ylabel('hd(n)')
subplot(2,2,2); stem(n,w_kai);title('Kaiser Window')
axis([-1 M 0 1.1]); xlabel('n'); ylabel('w(n)')
subplot(2,2,3); stem(n,h);title('Actual Impulse Response')
axis([-1 M -0.2 0.8]); xlabel('n'); ylabel('h(n)')
subplot(2,2,4);plot(w/pi,db); axis([0 1 -80 10]);
title('Magnitude Response in dB');grid;
xlabel('frequency in pi units'); ylabel('Decibels')
function F = ideal_lp(wc,M)
t = (M-1)/2;
x = 0: (M-1);
m = x - t + eps;
F = sin(wc*m) ./ (pi*m);
end
function [db,mag,pha,grd,w] = freqz_m(b,a)
[H,w] = freqz(b,a,1000,'whole');
H = (H(1:1:501))'; w = (w(1:1:501))';
mag = abs(H);
db = 20*log10((mag+eps)/max(mag));
pha = angle(H);
grd = grpdelay(b,a,w);
end
