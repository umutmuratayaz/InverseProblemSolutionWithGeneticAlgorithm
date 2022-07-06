
%% Crossover Function for Genetic Algorithm

function [tempPop] = gaCrossover(tempPop,Chr,pCross,Gen)
    pairs=randperm(Chr);
        for l=1:(Chr/2)
            pr1idx=pairs(2*l-1);
            pr2idx=pairs(2*l);
            pr1=tempPop(pr1idx,:);
            pr2=tempPop(pr2idx,:);
            rN=unifrnd(0,1);
            
            %Crossover for just Resistivity parameters(genes) 
            if rN<pCross
                cP=unidrnd((ceil(Gen/2))-1);
                ABC=pr1(cP+1:(ceil(Gen/2)));
                pr1(cP+1:(ceil(Gen/2)))=pr2(cP+1:(ceil(Gen/2)));
                pr2(cP+1:(ceil(Gen/2)))=ABC;
                tempPop(pr1idx,1:(ceil(Gen/2)))=pr1(1:ceil(Gen/2));
                tempPop(pr2idx,1:(ceil(Gen/2)))=pr2(1:ceil(Gen/2));
            end

            %Crossover for just Thickness parameters(genes)
            if rN<pCross
                cP2=unidrnd(floor(Gen/2)-1);
                DEF=pr1((ceil(Gen/2)+1):(ceil(Gen/2)+cP2));
                pr1((ceil(Gen/2)+1):(ceil(Gen/2)+cP2))=...
                    pr2((ceil(Gen/2)+1):(ceil(Gen/2)+cP2));
                pr2((ceil(Gen/2)+1):(ceil(Gen/2)+cP2))=DEF;
                tempPop(pr1idx,((ceil(Gen/2)+1):end))=...
                    pr1((ceil(Gen/2)+1):end);
                tempPop(pr2idx,((ceil(Gen/2)+1):end))=...
                    pr2((ceil(Gen/2)+1):end);
            end
    end
end

