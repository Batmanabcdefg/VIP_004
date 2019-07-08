clc;

clear all;

data=dlmread('fpprts.dat');

[m,n]=size(data);

count=0;

for i=1:1:m

    count=count+1;
    
timeval(count)=data(i,1);

reacval(count)=data(i,4);

end

counter=0;

for j=1:21:count
    
    counter=counter+1;
    
    time(counter)=timeval(j);
    
    stressval(j)=0;
    
for k=j:1:j+20
    
    stressval(j)=stressval(j) + reacval(k);
    
end

stress(counter)=stressval(j)/40.0e-09;

end

strain=time./1000;

plot(strain,-stress,'k');

xlabel('Nominal Strain');

ylabel('Nominal Stress');

%axis([0 0.12 0 2.2e+09]);


    