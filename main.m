%% Traditional Inversion for Vertical Electrical Sounding 
% It was used "Damped Least-Squares Method" oth.called "Marquart-Levenberg"
% It was developed from Ekinci and Demirci (2018) by Murat Ayaz (2019)...
%                within the scope of the undergraduate graduation project.

close all; clear all; clc; format long;

%% Input Parameters and Data

filename = 'InvInOut.xls';             %Input MS Excel File
sheet='Forw';                     %MS Excel Sheet
ab2=xlsread(filename,sheet,'E6:E60');  %Observed AB/2 Values  
roa=xlsread(filename,sheet,'F6:F60');  %Observed Resistivity Values
r=xlsread(filename,sheet,'I6:I15');    %Prediction Resistivity parameters  
t=xlsread(filename,sheet,'J6:J15');    %Prediction Thickness parameters
maxiteration=xlsread(filename,sheet,'K7');  %Number of Iteration
kr=xlsread(filename,sheet,'K10');          %Acceptable difference
realRo=xlsread(filename,sheet,'C6:C16');%Real Res.Value (Just Forward Sol.)
realH=xlsread(filename,sheet,'D6:D16'); %Real Thickness (Just Forward Sol.)


%% Calculations

x = ab2; iteration = 1; r=r'; t=t';
m = [r t]; rinitial = r; tinitial = t;
lr = length(r); lt = length(t); 
dfit = 1;
tic;                                    % Calculation Time Starter
while iteration<maxiteration
    r = m(1:lr);
    t = m(1+lr:lr+lt);
        for i = 1:length(x)
        s = ab2(i);
        [g] = marqLevInvMOD (r,t,s);
        roa1(i,:) = g;
        end
    e1 = [log(roa)-log(roa1)];
    dd = e1;
    misfit1 = e1'*e1;
    [A] = marqLevInvJAC(ab2,x,r,t,lr,lt,roa,roa1);
    [U S V] = svd(A,0);
    ss = length(S);
    say = 1;
    k = 0;
    while say<ss
    diagS = diag(S);
    beta = S(say)*(dfit^(1/say));
        if beta<10e-5
        beta = 0.001*say;
        end
        for i4 = 1:ss
        SS(i4,i4) = S(i4,i4)/(S(i4,i4)^2+beta);
        end
    dmg = V*SS*U'*dd;
    mg = exp(log(m)+dmg');
    r = mg(1:lr);
    t = mg(1+lr:lr+lt);
        for i5 = 1:length(x)
        s = ab2(i5);
        [g] = marqLevInvMOD (r,t,s);
        roa4(i5,:) = g;
        end
    e2 = [log(roa)-log(roa4)];
    misfit2 = e2'*e2;
        if misfit2>misfit1
            ('Beta control');
            say = say+1;
            k = k+1;
                if k==ss-1
                iteration = maxiteration;
                say = ss+1;
                end
        else
            say = ss+1;
            m = mg;
            dfit = (misfit1-misfit2)/misfit1;
            iteration = iteration+1;
            a = iteration;
                if dfit<kr
                iteration = maxiteration;
                say = say+1;
                end
        end
    end
       
end
TimeMarq=toc;                           %Calculation Time Finisher

observed = roa;
calculated = roa4;
format bank;
rms = norm(roa4-roa) /sqrt(length(roa));


%% Plotting The Figure 

subplot(1,2,1),
loglog(x,roa,'k.','Marker','.','MarkerSize',20);
hold on
loglog(x,roa4,'b','Marker','s','MarkerSize',15,'LineWidth',1.2);
axis([1 1000 1 1000])
xlabel('AB/2 (m)');
ylabel('Apparent Resistivity (Ohm-m)');
pause(0.001)
legend('G�zlenen De�erler','Hesaplanan De�erler','FontSize',12)
dim=[0.00470370370370371 0.682615629984051 0.0851111111111111...
    0.0685805422647529];
str=sprintf('%6s\n %.2f','RMS(%)',rms);
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',10);
dim=[0.00470370370370371 0.768740031897927 0.0851111111111111...
    0.0685805422647529];
str=sprintf('%8s\n %d','�terayon',maxiteration);
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',10);
dim=[0.00562962962962963 0.854864433811802 0.0851111111111111...
    0.0685805422647529];
str=sprintf('%13s\n %.4f','Hes.Zaman(sn)',TimeMarq);
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',10);
hold off
grid on

rr = [0,r];
tt = [0,cumsum(t),max(t)*10];
subplot(1,2,2),
stairs(rr,tt,'b-','LineWidth',2);
grid on
rrr = [0,rinitial];
ttt = [0,cumsum(tinitial),max(tinitial)*10];
hold on;
subplot(1,2,2),
stairs(rrr,ttt,'k--','LineWidth',2);
grid on
set(gca,'Ydir','reverse');
set(gca,'Xscale','log');
xlabel('Resistivity (Ohm)');
ylabel('Depth (m)');
Hor=[rr rrr];
xlim([1 max(Hor)*1.2]);
if max(tinitial)>max(tt)
ylim([0 sum(tinitial)*1.3]);
elseif max(tinitial)<max(tt)
ylim([0 max(tt)]);
end
subplot(1,2,2),
rRo=[0,realRo']; rH=[0,cumsum(realH'),max(realH')*10];
stairs(rRo,rH,'g-.','LineWidth',2);
legend('Hesaplanan Parametreler','�n-Kestirim Parametreleri',...
    'Ger�ek Parametreler','Location','southwest','FontSize',12)

%Output figure to .jpeg file
set(gcf,'Position', [0 0 1280 720]);
picName = sprintf('%s',sheet);
print(picName,'-r600','-djpeg');
grid on


%% Output to .xls File

xlswrite(filename,ab2,sheet,'O6');
xlswrite(filename,r',sheet,'M6');
xlswrite(filename,t',sheet,'N6');
xlswrite(filename,roa4,sheet,'P6');
xlswrite(filename,rms,sheet,'K13');
xlswrite(filename,TimeMarq,sheet,'K16');

