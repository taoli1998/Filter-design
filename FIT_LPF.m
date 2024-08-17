wp = 0.2*pi; ws = 0.3*pi; tr_width = ws - wp;
M = ceil(6.6*pi/tr_width) + 1;
n=0:1:M-1;
wc = (ws+wp)/2; % Ideal LPF cutoff frequency
hd = ideal_lp(wc,M); w_ham = (hamming(M))'; h = hd .* w_ham;
[db,mag,pha,grd,w] = freqz_m(h,1); delta_w = 2*pi/1000;
Rp = -(min(db(1:1:wp/delta_w+1))); % Actual Passband Ripple
As = -round(max(db(ws/delta_w+1:1:501))); % Min Stopband attenuation
% plots
subplot(2,2,1); stem(n,hd); title('Ideal Impulse Response')
axis([0 M-1 -0.1 0.3]); xlabel('n'); ylabel('hd(n)')
subplot(2,2,2); stem(n,w_ham);title('Hamming Window')
axis([0 M-1 0 1.1]); xlabel('n'); ylabel('w(n)')
subplot(2,2,3); stem(n,h);title('Actual Impulse Response')
axis([0 M-1 -0.1 0.3]); xlabel('n'); ylabel('h(n)')
subplot(2,2,4); plot(w/pi,db);title('Magnitude Response in dB');grid
axis([0 1 -100 10]); xlabel('frequency in pi units'); ylabel('Decibels')
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
