a1 = zeros(1000,1);
a2 = zeros(1000,1);
for i = 1:950
    a1(i) = sin(2*pi*i/10);
    a2(i+50) = sin(2*pi*i/10);
end

for i =1:850
    a1(i) = a1(i) + 0.5*sin(2*pi*i/10);
    a2(i+150) = a2(i+150) + 0.5*sin(2*pi*i/10);
end

t=-999:999;
cc = xcorr(a1,a2);

subplot(3,1,1)
plot(a1)
xlabel('Time/s', fontsize=15)

subplot(3,1,2)
plot(a2)
xlabel('Time/s', fontsize=15)

subplot(3,1,3)
plot(t, cc)
xlabel('Time/s', fontsize=15)