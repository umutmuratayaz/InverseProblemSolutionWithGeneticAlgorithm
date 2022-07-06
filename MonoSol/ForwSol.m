%% Forward Solution for Vertical Electrical Sounding
% It was developed from Oru� (2012) by Murat Ayaz (2019) ...
%                within the scope of the undergraduate graduation project.

close all; clear all; clc; format long; 

%% Input Parameters and Data

filename = 'InvInOut.xls';             %Input MS Excel File
sheet='Forw';                     %MS Excel Sheet
ab2=xlsread(filename,sheet,'A6:A60');  %AB/2 Values 
ro=xlsread(filename,sheet,'C6:C16');   %Resistivity Values
h=xlsread(filename,sheet,'D6:D16');    %Thickness of Layers

%Filter Coefficients (Guptasarma, 1982)
a(1)=-0.17445;a(2)=0.09672;a(3)=0.36789;a(4)=0.63906;a(5)=0.91023;
a(6)=1.18143;a(7)=1.45257; 
fi(1)=0.1732;fi(2)=0.2945;fi(3)=2.1470;fi(4)=-2.1733;
fi(5)=0.6646;fi(6)=-0.1215;fi(7)=0.0155;

%Number of Data and Layers
nv=length(ab2); nt=length(ro);    

%% Calculations 

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

roa(ii)=top;

end

%% Output to .xls File

xlswrite(filename,ab2,sheet,'E6');
xlswrite(filename,roa',sheet,'F6');

%% Plotting The Figure 

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.407407407407407 0.125996810207337 ...
    0.575852504486866 0.759170653907496]);
hold(axes1,'on');

% Create loglog
loglog(ab2,roa,'DisplayName',' G�zlenen (ohm.m)','Parent',axes1,...
    'MarkerSize',7,...
    'Marker','o',...
    'LineWidth',2);

% Create ylabel
ylabel('roa (ohm.m)');

% Create xlabel
xlabel('ab/2 (m)');

box(axes1,'on');
grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'XMinorTick','on','XScale','log','YMinorTick','on',...
    'YScale','log');
set(axes1,'ylim',[min(roa)*0.4 max(roa)*2]);

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.485516856040417 0.82028980209614 ...
    0.14298978359365 0.033381550333531]);

% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.0518518518518519 0.543859649122807...
    0.289814814814815 0.341307814992026]);
hold(axes2,'on');

% Create bar
bar(ro,'DisplayName','�zdiren� De�erleri (ohm)','Parent',axes2,...
    'FaceColor',[0 0.447058823529412 0.741176470588235]);

% Create ylabel
ylabel('�zdiren� (ohm-m)');
ylim([0 max(ro)*1.4]);

% Create xlabel
xlabel('Tabaka No');

box(axes2,'on');
% Set the remaining axes properties
set(axes2,'XMinorTick','on','XTick',[1 2 3 4],'YGrid','on','YMinorGrid',...
    'on','YMinorTick','on','ZMinorTick','on');
% Create legend
legend2 = legend(axes2,'show');
set(legend2,'Visible','off',...
    'Position',[0.215209884205591 0.71766822242299 ...
    0.190377931966427 0.0329662774804008]);

% Create axes
axes3 = axes('Parent',figure1,...
    'Position',[0.0527777777777778 0.124401913875598 ...
    0.289814814814815 0.301435406698564]);
hold(axes3,'on');

% Create bar
bar(h,'DisplayName','Tabaka Kal�nl�klar� (m)');

% Create ylabel
ylabel('Kal�nl�k (m)');
ylim([0 max(h)*1.5]);

% Create xlabel
xlabel('Tabaka No');

box(axes3,'on');
% Set the remaining axes properties
set(axes3,'XTick',[1 2 3],'YGrid','on','YMinorGrid','on');
% Create legend
legend3 = legend(axes3,'show');
set(legend3,'Visible','off',...
    'Position',[0.179367862664603 0.448185620512542 ....
    0.16946937314803 0.0329662774804012]);

% Figure Output
set(gcf,'Position', [0 0 1280 720]);
picName = sprintf('%sSol',sheet);
print(picName,'-r600','-djpeg');
