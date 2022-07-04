
%% Fitness Value Calculation for Genetic Algorithm

function [Fitn,bRoa] = gaMOD(Chr,ab2,nt,pop,Gen,Fitn,Olculen,bRoa)

    for i=1:Chr   %Determination of Fitness Value for every chromosome
    nv=size(ab2,1);
    
        %Assign value from Genes to Resistivity and Thickness parameters    
        for j=1:nt
        ro(j)=pop(i,j); 
        end
        for k=1:nt-1
        h(k)=pop(i,ceil(Gen/2)+k);
        end

        %Filter Coefficients for Forward Solution (Guptasarma, 1982)
        a(1)=-0.17445;a(2)=0.09672;a(3)=0.36789;
        a(4)=0.63906;a(5)=0.91023;
        a(6)=1.18143;a(7)=1.45257; 
        fi(1)=0.1732;fi(2)=0.2945;
        fi(3)=2.1470;fi(4)=-2.1733;
        fi(5)=0.6646;fi(6)=-0.1215;fi(7)=0.0155;

        for ii=1:nv
        l=ab2(ii); 
        top=0;
            for r=1:7
            lamda(r)=10^(a(r)-log10(l));
                 for jj=nt-1:-1:1
                     if jj==nt-1
                     uu(jj)=exp(-2*h(jj)*lamda(r)); 
                     kk(jj)=(ro(jj)-ro(jj+1))/(ro(jj)+ro(jj+1));
                     t(jj)=ro(jj)*(1-kk(jj)*uu(jj))/(1+kk(jj)*uu(jj));
                     else
                     u(jj)=exp(-2*h(jj)*lamda(r));
                     w(jj)=ro(jj)*(1-u(jj))/(1+u(jj));
                     t(jj)=(w(jj)+t(jj+1))/(1+(w(jj)*t(jj+1))/ro(jj)^2);
                     end
                 end
            t(r)=t(1);
            top=top+fi(r)*t(r);
            end

     %Determined app.resistivity for every ab2 point 
        roa(ii)=top;

     %Square difference of Determined-Observed value for "ii." ab2 point 
        Fitn(i)=Fitn(i)+((roa(ii)-Olculen(ii))^2);

        end

    %The square root average  
    Fitn(i)=sqrt(Fitn(i)/nv);

        while i~=1    %Resistivity Values for Best Candidate
            if(Fitn(i)<Fitn(i-1))
            bRoa=roa;
            end
        break
        end
    end %i:Chr 
end

