%%%Fitting a parabola with 1 free parameter 'A' (y = Ax^2) to parabolic 
%%%data with some noise added to it (noise in range -delta -> delta).  Plot
%%%The error in the fit parameter 'A' as a function of the noise range
%%%'delta', the number of fitting points 'N'.  The error plotted here is an
%%%ensemble average of 'T' separate fitting tests, for rigor%%% 

A = 1;
N = 10;
T = 1000;
%sig = 0.1;
Legend=cell(5,1); 
for j = 1:5
    delta = 0.05*j;
    Eav= zeros(50,1);
    for N = 3:53
        E = zeros(T,1);
        for i = 1:T
            x = linspace(0,1,N);
            
            modelfun = @(a,X)(a*X.^2);
            yn = modelfun(A,x) + 2*delta*rand(1,N)-delta;
            opts = statset('nlinfit');
            opts.RobustWgtFun = 'bisquare';
            yf = nlinfit(x,yn,modelfun,A,opts);
            
            %y = A*x.^2;
            %yn = y + 2*sig*rand(1,N)-sig;
            %p = polyfit(x,yn,2);
            %yf = (p(1)*x.^2 + p(2)*x + p(3));
            E(i) = abs(A-yf);%(sum((y - yf).^2)/N).^0.5; %((A-p(1))*2 + p(2).^2 + p(3).^2).^0.5;
        end
        Eav(N-2) = mean(E);  
    end
    
    plot(3:53,100*Eav)
    Legend{j}= ['\delta = ',num2str(delta)];
    hold on;
end
legend(Legend);
xlabel('# of datapoints')
ylabel('abs(A-p)/A')%ylabel('Standard deviation/Max force value (%)')