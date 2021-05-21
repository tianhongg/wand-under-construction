a=load('Grids.ini');
figure
plot(a)


for i=2:149
   a(i)=(a(i-1)+a(i+1))/2;
end



figure
plot(a)

figure
plot(b)

for i=2:75
    a(i)=(a(i-1)+a(i+1))/2;
end
for i=149:-1:76
    a(i)=(a(i-1)+a(i+1))/2;
end
% 
% for i=1:75
%     ave=(a(i)+a(151-i))/2;
%     a(i)=ave;
%     a(151-i)=ave;
% end
% 
% figure
% plot(a(1:75)-fliplr(a(76:150)))
% 
