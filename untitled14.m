clear all
close all
p=load('Denn_4_.dg');

figure
pcolor(p)
shading interp
% caxis([0.99999,1.00001])


% figure
% hold on;
% for i=0:8
% p=load(['Parts_',num2str(i),'_.dg']);
% plot3(p(:,1),p(:,2),p(:,3),'.')
% 
% end