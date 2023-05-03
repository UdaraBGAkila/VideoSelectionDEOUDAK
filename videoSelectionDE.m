function [] = videoSelectionDE()
clc;
clear;

xx = [1 2 3 4 5 6 7 8 9 10; 800 700 650 750 600 900 950 875 1050 1500;121 95 85 100 78 125 130 128  135 120];
dataMat = xx'
n = 10;
[chromosome] = init_pop(n)

for itr=1:100
  for i = 1:n
    [donorVector] = doMutation(chromosome,i);
    [trialVector] = doCrossover(chromosome,donorVector,i);
    [targetVector] = doSelection(trialVector,chromosome,dataMat,i);
    chr(i,:) = targetVector;
  end
[fitness cdSize] = evaluateFitness(chromosome,dataMat,n)

chromosome = chr;

[bestChromoFitness ind]= max(fitness);
cdSize(ind);
fitness(ind);
chromosome(ind,:);
xplot=itr;
yplot=fitness(ind);
figure(1)
plot(xplot,yplot,'*','MarkerSize',5);
grid on;
hold on

end

[fitness cdSize] = evaluateFitness(chromosome,dataMat,n);
[bestChromoFitness ind]= max(fitness)
cdSize(ind)
chromosome(ind,:)


function [chromosome] = init_pop(n)
for i=1:n
chromosome(i,:) = randi([0, 1], 1,10);
end

function [fitness cdSize] = evaluateFitness(chromosome,dataMat,n)
playTime = dataMat(:,3);
fileSize = dataMat(:,2);
lenPlayTimeData = length(playTime);
for j=1:n
    sum = 0;
    sumFilesize = 0;
    for k=1:lenPlayTimeData
        sum = sum +(chromosome(j,k).*playTime(k));
        sumFilesize = sumFilesize + (chromosome(j,k).*fileSize(k));
    end
    fitness(j) = sum;
    cdSize(j) = sumFilesize;
   
    if cdSize(j)<=4500
        fitness(j) = sum;
    else
        fitness(j) = 0.25*sum;
    end
end


function[trialVector] = doCrossover(chromosome,donorVector,i)
targetVector = chromosome (i,:);
CR = 0.1; % recombination probability
I = randi(10);

  for k=1:10
    r = rand;
    
    if (r > CR) && (I ~= k)
      U(:,k) = targetVector(:,k);
    else
      U(:,k) = donorVector(:,k);
    end
  end
trialVector = U;


function[donorVector] = doMutation(chromosome,i)
F = 0.1;
N = 1:10; % Generate numbers from 1 to 10 assign them to N
N(i) = []; % Removes ith number from N
p = N(randperm(numel(N),3)); % Take 3 unique numbers from N and assign them to p
a = p(1);
b = p(2);
c = p(3);

V = chromosome (a,:) +F *(chromosome(b,:) -chromosome(c,:));

% To make all values either 1 or 0
  %for x=1:10 
      %if V(x) <= 0
        %V(x) = 0;
      %else
        %V(x) = 1;
      %end
  %end
  
% To make all values either 1 or 0  
  for x=1:10 
      if V(x) == 0
        V(x) = V(x);
      elseif V(x) == 1
        V(x) = V(x);  
      else
        V(x) = randi([0,1]);
      end
  end
donorVector = V;  


function [targetVector] = doSelection(trialVector,chromosome,dataMat,i)
targetVector = chromosome (i,:);
playTime = dataMat(:,3);
fileSize = dataMat(:,2);
lenPlayTimeData = length(playTime);

sumU = 0;
sumUFilesize = 0;
sumX = 0;
sumXFilesize = 0;

      for k=1:lenPlayTimeData
          sumU = sumU +(trialVector(:,k).*playTime(k));
          sumX = sumX +(targetVector(:,k).*playTime(k));
          sumUFilesize = sumUFilesize + (trialVector(:,k).*fileSize(k));
          sumXFilesize = sumXFilesize + (targetVector(:,k).*fileSize(k));
      end
      
      cdSizeU = sumUFilesize;
      cdSizeX = sumXFilesize;
     
      if cdSizeU<=4500
          fitnessU = sumU;
      else
          fitnessU = 0.25*sumU;
      end
      
      if cdSizeX<=4500
          fitnessX = sumX;
      else
          fitnessX = 0.25*sumX;
      end
      
      if fitnessU >= fitnessX
          X = trialVector;
      else
          X = targetVector;
      end
targetVector = X;
