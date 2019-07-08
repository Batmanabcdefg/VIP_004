%Generate random numbers with given mean and standard deviation

% clc;
% 
% clear all;
% 
% close all;
% 
% mean_crss=1.12;
% 
% std_dev=0.01*mean_crss;
% 
% dist_crss=mean_crss+std_dev*randn(1200,1);
% 
% mean_dist_crss=mean(dist_crss)
% 
% std_dev_dist_crss=std(dist_crss)
% 
% std_dev_percent=(std_dev_dist_crss/mean_dist_crss)*100

% data=load('cohes.dat');
% 
% for i=1:1:100
% 
% j=8*(i-1)+1:1:8*(i-1)+8
%     
% data1(j)=data(i,:);
% 
% end

data1=load('cohes1.dat');

count=0;

for i=1:8:1200
    
    k=1;
    
    count=count+1;
  
for j=i:1:i+7
    
    data(count,k)=data1(j);
    
    k=k+1;
    
end

end
