quals = csvread('C:\Documents and Settings\hutter\Desktop\log-rep-qualities.csv');
ccs = csvread('C:\Documents and Settings\hutter\Desktop\log-rep-ccs.csv');
lls = csvread('C:\Documents and Settings\hutter\Desktop\log-rep-lls.csv');
rmses = csvread('C:\Documents and Settings\hutter\Desktop\log-rep-rmses.csv');

figure(1)
subplot(2,2,1)
confplot(1:100, mean(lls(:,1:100)), std(lls(:,1:100))); hold on; title('ll')
subplot(2,2,2)
confplot(1:100, mean(ccs(:,1:100)), std(ccs(:,1:100))); hold on; title('cc')
subplot(2,2,3)
confplot(1:100, mean(rmses(:,1:100)), std(rmses(:,1:100))); hold on; title('rmse')
subplot(2,2,4)
confplot(20:200, mean(quals(:,20:200)), std(quals(:,20:200))); hold on; title('incumbent')

%for i=10, figure(i+4); plot(1:100, rmses(i,1:100)); hold on; title('rmse'); end