clc
clear
data='mls/mlsUp3'; dep='sold'; ind='age     lot     sqft    ';
beta=zeros(4,1);
beta=probit(data,dep,ind,beta);