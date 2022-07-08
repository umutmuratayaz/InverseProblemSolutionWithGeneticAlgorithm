%% Natural Selection Function for Genetic Algorithm
% Roulette wheel choise  was used for production of new generations. 

function [tempPop] = gaNaturalSelection(Fitn,pop,Chr)


Fitn=1./Fitn;  %We are trying to make the difference minimum

sumFitn=sum(Fitn);
prob=Fitn/sumFitn;

cProb=prob;

        for k=2:Chr   %cumulative sum of probabilities
            cProb(k)=cProb(k-1)+prob(k);
        end

rN=unifrnd(0,1,[Chr,1]);

tempPop=pop;

for kk=1:Chr
    idx=find(rN(kk)<cProb,1);
    tempPop(kk,:)=pop(idx,:);
end



end

