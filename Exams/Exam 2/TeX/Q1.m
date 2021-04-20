 % Computer Take-Home 2: Question 1
 % Metrics III
 %Oscar Martinez
 
 %Diary
 diary Q2_Output_Oscar_Martinez.txt
 
 %Introduction
 fprintf('--------------------------------------------------------------\n');
 fprintf('Oscar Martinez \t Take-Home 2: Question 1 \t Metrics III\n');
 fprintf('--------------------------------------------------------------\n');
 
 %Exam 2 Question 1:
 data='mlsUp3'; bin='sold'; dep='days'; ind='vac     lap     ';
 theta=[-3; 0; 0; 1];
 theta=weibull(data,dep,bin,ind,theta);
 
 %closing output
 diary off
