%% Function for Genetic Algorithm  

function [roGA,thiGA,TimeGA,bRoa,bV] = MultiGA(RoDown,RoUp,...
    ThiDown,ThiUp,Chr,Gen,maxIteration,pCross,pMutation,pElt,...
    nt,sheet5,cnk,ab2,Olculen,realRo,realH)     

close all;

bV=1000000; iteration=0; %First Fitness Value and Iteration Counter
cnk=20; %%Is the population in search space? A check in each cnk steps.
cn=0;

% First Population Creation
pop=[unifrnd(RoDown,RoUp,[Chr,ceil((Gen/2))]),...
    unifrnd(ThiDown,ThiUp,[Chr,floor(Gen/2)])];

% Calculations
tic;
for ijk=1:maxIteration

    % Checking whether the population is in search space.  
    if mod(iteration,cnk)==0
    [pop] = gaCheck(pop,Gen,RoUp,ThiUp,RoDown,ThiDown,Chr);
    end

%Creating Fitness Value Maxrix
Fitn=zeros(Chr,1);

%Fitness Value Calculations
[Fitn,bRoa] = gaMOD(Chr,ab2,nt,pop,Gen,Fitn,Olculen);
   
%Find The Best Solution
iteration=iteration+1;
    if (min(Fitn)<bV)
    bV=min(Fitn);
    idxBs=find(Fitn==bV);
    BestSolution=pop(idxBs,:);
    end

%The Best Value for current iteration
bVI(iteration)=bV;

%It leave of elite from population before Crossover and Mutation 
EltLeng=floor(length(pop)*pElt);
[ElitFitn idx]=mink(Fitn,EltLeng,1);
ElitPop=pop(idx,:);

%Natural Selection 
[tempPop] = gaNaturalSelection(Fitn,pop,Chr);

%Crossover Operator  
[tempPop] = gaCrossover(tempPop,Chr,pCross,Gen);

%Mutation Operator
[pop] = gaMutation(tempPop,pMutation,Chr,Gen,...
    RoUp,RoDown,ThiUp,ThiDown);

%Return of Elite to population
pop(idx,:)=ElitPop;

end
TimeGA=toc;

% Plotting The Figure 

bRoa=bRoa';
roGA=BestSolution(1:(ceil(Gen/2)));
thiGA=BestSolution((ceil(Gen/2))+1:end);
FroGA=[0 roGA]; FthiGA=[0, cumsum(thiGA),max(thiGA)*10];

subplot(2,2,2); 
stairs(1:maxIteration,bVI,'Marker','.','MarkerSize',6,'LineWidth',2);
legend('Uygunluk Deðerlerinin Ýterasyonla Deðiþimi','FontSize',12)
xlabel('Ýterasyon Adýmý'); ylabel('Uygunluk(%)');
grid on

subplot(2,2,[1 3]);
loglog(ab2,bRoa,'b','Marker','s','MarkerSize',15,'LineWidth',1.2)
ylim([1 1000]); xlim([1 1000]);
ylabel('Apparent Resistivity (Ohm-m)');
xlabel('AB/2 (m)');

dim=[0.00470370370370371 0.682615629984051 0.0851111111111111...
    0.0685805422647529];
str=sprintf('%6s\n %.2f','RMS(%)',bV);
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',10);

dim=[0.00470370370370371 0.768740031897927 0.0851111111111111...
    0.0685805422647529];
str=sprintf('%8s\n %d','Ýterayon',maxIteration);
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',10);

dim=[0.00562962962962963 0.854864433811802 0.0851111111111111...
    0.0685805422647529];
str=sprintf('%13s\n %.4f','Hes.Zaman(sn)',TimeGA);
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',10);

grid on
hold on 
loglog(ab2,Olculen,'k.','Marker','.','MarkerSize',20);
ylim([1 1000]); xlim([1 1000]);
legend('Hesaplanan Deðerler','Gözlenen Deðerler','FontSize',12)
hold off

subplot(2,2,4);
stairs(FroGA,FthiGA,'b-','LineWidth',2);
set(gca,'Ydir','reverse');
set(gca,'Xscale','log');
axis([1 max(FroGA)*1.2 0 sum(thiGA)*1.3]);

xlabel('Resistivity (Ohm)');
ylabel('Depth (m)');
grid on

hold on 
FrealRo=[0,realRo'];FrealH=[0 cumsum(realH'),max(realH')*10];
subplot(2,2,4);
stairs(FrealRo,FrealH,'r--','LineWidth',2);
legend('Hesaplanan Parametreler','Gerçek Parametreler'...
    ,'Location','southwest','FontSize',12);

%Output figure to .jpeg
set(gcf,'Position', [0 0 1280 720]);
pictName=sprintf('%s',sheet5);
print(pictName,'-r200','-djpeg')

end


