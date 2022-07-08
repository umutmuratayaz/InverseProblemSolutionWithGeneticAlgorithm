%% Mutation Function for Genetic Algorithm

function[pop] = gaMutation(tempPop,pMutation,Chr,Gen,...
    RoUp,RoDown,ThiUp,ThiDown)


nR=0.06; nT=0.03;   %Neighborhood Value. Some special methods are required 
                    %                   for mutation with "value coding". 
                    %The neighborhood value,  
                    %       together with the random number, 
                    %          contributes to the change in the parameter. 
                    %In the event of mutation: it is directly 
                    %         proportional to the amount of change.
  

rN=unifrnd(0,1,[Chr,Gen]);

for m=1:Chr
    
    %Mutation for Resistivity parameters (randomValue x 5% change in search space) 
    for mm=1:(ceil(Gen/2))    
      if (rN(m,mm)<pMutation)
        rN2=unifrnd(-1,1);
        tempPop(m,mm)=tempPop(m,mm)+(rN2*nR*(RoUp-RoDown));
      end
    end
    
    
    %Mutation for Thickness parameters (randomValue x 2% change in search space) 
    for mmm=(ceil(Gen/2)+1):Gen 
       if (rN(m,mmm)<pMutation)
        rN3=unifrnd(-1,1);
        tempPop(m,mmm)=tempPop(m,mmm)+(rN3*nT*(ThiUp-ThiDown));
       end
    end
       
    pop=tempPop;

end

end
