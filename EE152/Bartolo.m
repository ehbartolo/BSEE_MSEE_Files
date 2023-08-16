%%%%%%%%%%%%%%%%%%%%%%  EE 152 Machine Problem Decoupled Harmonic Power Flow %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  Submitted by: Eldion Vincent H. Bartolo  2010-28167  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear                                                                  %clears all variables
fprintf('\n\n                    Decoupled Harmonic Power Flow\n\n\n');
[num,txt,raw]=xlsread('HPFinput.xls','global');            %read the global sheet
basemva = num(1,1);
accuracy = num(1,2);
accel = num(1,3);
maxiter = num(1,4);
maxharm = num(1,5);

[num,txt,raw]=xlsread('HPFinput.xls','busdata');           % read the busdata sheet
busdataold = num;
nbus = length(busdataold(:,1));
Plload = (busdataold(:,5).')/basemva;                                  % power of linear loads in p.u
Qlload = (busdataold(:,6).')/basemva;                                  % reactive power of linear loads in p.u

[num,txt,raw]=xlsread('HPFinput.xls','linedata');          % read the line data
linedataold = num;                                                     % save linedata of fundamental frequency
%----------------------------- Solve for fundamental bus voltage of the n-bus sytem considering the nonlinear load ----------------------------------
busdatanew = busdataold;                                % transfer data
Pn = zeros(1,nbus);                                     % intialize to zeros all the nonlinear power 
Qn = zeros(1,nbus);                                     % initialize to zeros all nonlinear reactive power
[num,txt,raw]=xlsread('HPFinput.xls','nonlinearloads');    % read the nonlinearloads sheet
nlload = num;                                              
nnl = length(nlload(:,1));                              % number of nonlinear loads
for n=1:nnl                                             % Total nonlinear real and reactive power power per bus
    Pn(1,nlload(n,1)) = Pn(1,nlload(n,1)) + nlload(n,3);
    Qn(1,nlload(n,1)) = Qn(1,nlload(n,1)) + nlload(n,4);
end
busdatanew(:,5)=busdatanew(:,5)+ Pn.';                  % Add the real and reactive power of nonlinear loads to the busdata
busdatanew(:,6)=busdatanew(:,6)+ Qn.';
busdata = busdatanew;
linedata = linedataold;

lfybus                                                  % form the bus admittance matrix of fundamental freq.
lfnewton                                                % Load flow solution by Newton-Raphson method
phase = deltad*pi/180;                                  % bus phase angle in radians
V1 = Vm.*(cos(phase)+j*sin(phase));                     % Fundamental bus voltage (1 x nbus) with the nonlinear loads
%------------------ Compute the power loss of fundamental ---------------%
Hlineflow
Sloss = SLT;
%---------------- Fundamental current injection by each nonlinear loads  ---------------------------------------------------------------------------
I1 = zeros(1,nnl);                                      % fundamental current of each nonlinear load
for n =1:nnl
  Snload = (nlload(n,3)+j*nlload(n,4))/basemva;         % Apparent power in p.u of nonlinear load
  I1(1,n) = conj(Snload/V1(1,nlload(n,1)));             % I fundamental
    
end
%-------------------------------------- harmonic power flow  -----------------------------------------------------------------------------
[num,txt,raw]=xlsread('HPFinput.xls','currentinject');     % read the current inject sheet
nho = length(num(:,1));                                 % number of harm order read from .xls file
SumVhsquared = zeros(1,nbus);                           % initialize(|Vh[]|^2) = zero
curr = num;
dynamic = [];                                           % creates a dynamic array for Harmonic voltages in each bus
counter = 1;
Xaxis = [];
for b=2:nho                                             % for(b=2,b<=nho,b++)---> not including the fundamental
    linedatanew = linedataold;
    if num(b,1)<=maxharm                                % if(h<=H)
        Xaxis(1,b-1) = num(b,1);
        ylh = (Plload./(abs(V1).*abs(V1)))-j*(Qlload./((num(b,1)).*(abs(V1).*abs(V1))));   % load admittance due to harmonics
        linedatanew(:,4)= num(b,1).*(linedataold(:,4));                                    % new X = old X *h(harmonic number)
        linedatanew(:,5)= num(b,1).*(linedataold(:,5));                                    % new C = old C*h(harmonic number)
        linedata = linedatanew;
        lfybus                                                                             % create the y bus matrix(not including the load admittance)

        for m=1:nbus
            Ybus(m,m)=Ybus(m,m)+ylh(1,m);               % add load admittance to diagonal elements of Ybus
        end
        
        Ih =zeros(1,nbus);                              % create the Iharmonic (1 by nbus) matrix
        for q=1:nnl
            Chmag = (num(b,nlload(q,2)*2))/100;         % magnitude of C(h)for a specific nonlinear load
            Chphase = (num(b,nlload(q,2)*2+1))*pi/180;  % phase angle of C(h)in radians for a specific nonlinear load
            Ch = Chmag*(cos(Chphase)+j*sin(Chphase));   % C(h)
            Ih(1,nlload(q,1)) = Ih(1,nlload(q,1)) +Ch.*I1(1,q); % harmonic current ---> only buses with nonlinear loads have values 
            
        end
        Ih = Ih.';
        Vh = (inv(Ybus)*Ih).';                                 % Vh (1 by nbus)
        dynamic(counter,:) = abs(Vh);
        counter = counter+1;
        
        %--------- Calculate the power loss of harmonic frequencies ---------------------------------------
        V = Vh;                                                % Update V
        Hlineflow
        Sloss = Sloss+SLT;
        
        SumVhsquared = SumVhsquared+((abs(Vh)).*(abs(Vh)));    % Summation(|Vh|^2)from b= 2 to maxharm
    
    end
   
end
Plosstotal = real(Sloss);                                       % total real power losses
Qlosstotal = imag(Sloss);                                       % total reactive power losses
V1 = abs(V1);                           
THD = (sqrt(SumVhsquared)./abs(V1))*100;
Vrms = sqrt(SumVhsquared+(Vm.*Vm));
% for printing the Vrms
fprintf('\n\n        Output on buses\n');
fprintf('Bus No.   Vrms(in p.u)    VTHD\n')
for n=1:nbus
  fprintf('%4g', n), fprintf('     %9.3f',Vrms(1,n)),fprintf('     %9.3f',THD(1,n));
  fprintf('\n');
end
fprintf('\nTotal Power loss = ');
fprintf('%9.3f',Plosstotal);
fprintf('  Total Reactive Power loss = ');
fprintf('%9.3f',Qlosstotal);
fprintf('\n');

for counter2=1:nbus
if THD(1,counter2) == max(THD)
    busnumber =counter2; 
end
end

% Bar Plot
Yaxis = dynamic(:,busnumber).';
figure;
bar(Xaxis,Yaxis);
title('Bar Plot');
xlabel('Harm order');
ylabel('Harmonic Voltage');


% Time domain plot
Vplot_total = 0;
t=0:0.0001:0.02;

nho = nho-1;
for counter =1:nho
Vplot=Yaxis(1,counter)*sin(curr(counter+1,1)*2*pi*50*t);
Vplot_total=Vplot_total+Vplot;
end
figure;
Vplot_total=abs(V1(1,busnumber))+Vplot_total;
plot(t,Vplot_total)
title('Time Domain Plot')
xlabel('Time(s)');
ylabel('Equivalent Voltage Waveform(V)at highest THD');
grid;


