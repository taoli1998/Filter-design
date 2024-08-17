ws1 = 0.2*pi; wp1 = 0.35*pi; wp2 = 0.65*pi; ws2 = 0.8*pi;
tr_width = min((wp1-ws1),(ws2-wp2)); M = ceil(11*pi/tr_width) + 1;
n=0:1:M-1; wc1 = (ws1+wp1)/2; wc2 = (wp2+ws2)/2;
hd = ideal_lp(wc2,M) - ideal_lp(wc1,M);
w_bla = (blackman(M))'; h = hd .* w_bla;
[db,mag,pha,grd,w] = freqz_m(h,1); delta_w = 2*pi/1000;
Rp = -min(db(wp1/delta_w+1:1:wp2/delta_w)); % Actua; Passband Ripple
As = -round(max(db(ws2/delta_w+1:1:501))); % Min Stopband Attenuation
%Plots
subplot(2,2,1); stem(n,hd); title('Ideal Impulse Response')
axis([0 M-1 -0.4 0.5]); xlabel('n'); ylabel('hd(n)')
subplot(2,2,2); stem(n,w_bla);title('Blackman Window')
axis([0 M-1 0 1.1]); xlabel('n'); ylabel('w(n)')
subplot(2,2,3); stem(n,h);title('Actual Impulse Response')
axis([0 M-1 -0.4 0.5]); xlabel('n'); ylabel('h(n)')
subplot(2,2,4);plot(w/pi,db);axis([0 1 -150 10]);
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