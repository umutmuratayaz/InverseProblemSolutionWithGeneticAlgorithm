
%% Function of Checking whether the population is in search space.

function [pop] = gaCheck(pop,Gen,RoUp,ThiUp,RoDown,ThiDown,Chr)

    for i=1:Chr
        for ij=1:ceil(Gen/2)
            if pop(i,ij)>RoUp || pop(i,ij)<RoDown
                pop(i,ij)=unifrnd(RoDown,RoUp);
            end
        end 
        for ik=ceil(Gen/2)+1:Gen
            if pop(i,ik)>ThiUp || pop(i,ik)<ThiDown
                pop(i,ik)=unifrnd(ThiDown,ThiUp);
            end
        end 
    end
    
end

