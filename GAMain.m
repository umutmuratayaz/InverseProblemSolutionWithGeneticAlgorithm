%% Genetic Algorithm for Vertical Electrical Sounding
% It was developed by Murat Ayaz (2019) within the scope of the...
%                                        undergraduate graduation project.
%                                           Dokuz Eylul University, Ýzmir. 

close all; clear all; clc; format long;

%% Input Parameters and Data

filename = 'InvInOut.xls';
sheet='Forw';
Olculen=xlsread(filename,sheet,'F6:F60');   %Observed Data Input
ab2=xlsread(filename,sheet,'E6:E60');       %           
RoDown=xlsread(filename,sheet,'R6');        %Search Space for 
RoUp=xlsread(filename,sheet,'S6');          %     resistivity and 
ThiDown=xlsread(filename,sheet,'T6');       %         thickness 
ThiUp=xlsread(filename,sheet,'U6');         %               parameters
Chr=xlsread(filename,sheet,'R12');          %Chromosome(candidate solution)
Gen=xlsread(filename,sheet,'T12');          %Gene(parameter on cand. sol.)
maxIteration=xlsread(filename,sheet,'T9');  %Number of Iteration
pCross=xlsread(filename,sheet,'R15');       %Crossover probability  
pMutation=xlsread(filename,sheet,'T15');    %Mutation probability 
pElt=xlsread(filename,sheet,'S18');         %Elite Ratio  
nt=xlsread(filename,sheet,'R9');            %Number of Layers

bV=1000000; iteration=0; %First Fitness Value and Iteration Counter
cnk=20;%%Is the population in search space? A check in each cnk steps.
cn=0; 

%% First Population Creation

pop=[unifrnd(RoDown,RoUp,[Chr,ceil((Gen/2))]),...
    unifrnd(ThiDown,ThiUp,[Chr,floor(Gen/2)])];

%% Calculations

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

%% Plotting The Figure 

bRoa=bRoa';
roGA=BestSolution(1:(ceil(Gen/2)));
thiGA=BestSolution((ceil(Gen/2))+1:end);
FroGA=[0 roGA]; FthiGA=[0 cumsum(thiGA),max(thiGA)*10];

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
legend('Hesaplanan Parametreler','Location','southwest','FontSize',12);
xlabel('Resistivity (Ohm)');
ylabel('Depth (m)');
grid on

%Output figure to .jpeg
set(gcf,'Position', [0 0 1280 720]);
print('-r600','-djpeg','GAInv.jpeg')

%% Output to .xls File

roGA=BestSolution(1:(ceil(Gen/2)));
thiGA=BestSolution((ceil(Gen/2))+1:end);
xlswrite(filename,TimeGA,sheet,'W20');
xlswrite(filename,ab2,sheet,'Y6');
xlswrite(filename,roGA',sheet,'W6');
xlswrite(filename,thiGA',sheet,'X6');
xlswrite(filename,bRoa,sheet,'Z6');
xlswrite(filename,bV,sheet,'U20');
