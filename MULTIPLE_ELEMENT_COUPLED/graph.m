clc;

clear all;

data=load('fpprts.dat');

timeval=data(:,1);

reacval=data(:,4);

[m,n]=size(data);

counter=0;

for i=1:21:m
    
    counter=counter+1;
    
    time(counter)=timeval(i);
    
    stressval(i)=0;
    
for j=i:1:i+20
    
    stressval(i)=stressval(i) + reacval(j);
    
end

stress(counter)=stressval(i)/8.5e-09;

end

strain=time./1000;

plot(strain,-stress,'k');

xlabel('Nominal Strain');

ylabel('Nominal Stress');

axis([0 0.12 0 2.2e+09]);


    