clc;
clear;
clf;
%%%A simple Monte Carlo simulation of a traffic jam; cars are represented
%%%by '1's in a matrix, emply spots as '0's.  Each car has its own prefered
%%%speed, chosen randomly at initialization from the range Smin->Smax.
%%%Cars will move forward as long as they wont crash into a car in front of them in
%%%the next timestep; If they cannot move forward, they will move either
%%%left or right 1 lane to pass.  Highway has periodic boundaries, so cars
%%%are in a infinite loop.

%Create matrix of N+2 collumns (with boundaries); represent highway
%lanes
%0s represent open spaces, 1s represent cars.  originally at random
%positions, moving at speed random for each car.  car will lane switch(move
%to left, right adjacent position if open (left, right), bottom and top
%corners). forecasts F spaces ahead, decides to slow down if no spaces left
%or right to move.  Will always move left, if given choice (passing lane)
%Matrix is L by N+2 in size

%%%%%%%USER INPUT%%%%%%%%%%%%%%%
L = 50;
N = 3;
T = 5000;
Smin = 6;
Smax = 10;
Cweight = 0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%SETTING UP MATRIX%%%%%%%%

%proportion of initial positions filled by cars, before adding boundaries
%Problem with redistribution of car psoitions; done without replacement so
%> than Cweight percentile of cars
M = zeros(L,N+2);
M(1:round(Cweight*((N+2)*L))) = 1;
M = M(randi([1,(length(M))],L,N+2));
M(:,1) = -1*ones(L,1);
M(:,N+2) = -1*ones(L,1);
%Speed value for each car
M(end,end,2) = 0;
M(:,(2:N+1),2) = randi([Smin,Smax],L,N);
%finding car positions, setting speed of non- car to -1 (invalid speed)
M(:,:,2) = M(:,:,2).*(~(M(:,:,1)<1))-1;

%Matrix of x,y and speed of all cars
[Clist_x,Clist_y] = find(M(:,:,1)==1);
C = [Clist_x,Clist_y,diag(M(Clist_x,Clist_y,2))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%TIME STEP%%%%%%%%%%%%%%%%%
%move each car up by their speed

Avpos_x = zeros(T,1);
Avpos_y = zeros(T,1);
movevec = [];

for k = 1:T
    move = 0;
    for i = 1:length(C(:,1))
        if M((mod((C(i,1)-C(i,3)):(C(i,1)-2),L)+1),C(i,2),1) == zeros((C(i,3))-1,1)
            %M(mod((C(i,1)-C(i,3)):(C(i,1)-1)),C(i,2),1) == zeros((C(i,3)),1)
            Cpos = zeros(L,1);
            Cpos(C(i,1)) = 1;
            Cpos = circshift(Cpos,-(C(i,3)-1));
            M(:,C(i,2),1) = M(:,C(i,2),1) + Cpos;
            M(C(i,1),C(i,2),1) = 0;
            M(C(i,1),C(i,2),2) = -1;
            M(find(Cpos ==1),C(i,2),2) = C(i,3);
            move = move+1;
        elseif M(mod((C(i,1)-2):C(i,1),L)+1,C(i,2)-1,1) == zeros(3,1)
            M(C(i,1),(C(i,2)-1),1) = 1;
            M(C(i,1),(C(i,2)-1),2) = C(i,3);
            M(C(i,1),C(i,2),1) = 0;
            M(C(i,1),C(i,2),2) = -1;
            %move = move+1;
        elseif M(mod((C(i,1)-2):C(i,1),L)+1,C(i,2)+1,1) == zeros(3,1)
            M(C(i,1),(C(i,2)+1),1) = 1;
            M(C(i,1),(C(i,2)+1),2) = C(i,3);
            M(C(i,1),C(i,2),1) = 0;
            M(C(i,1),C(i,2),2) = -1;
            %move = move+1;
        end
    end
    
    %Matrix of x,y and speed of all cars
    [Clist_x,Clist_y] = find(M(:,:,1)==1);
    C = [Clist_x,Clist_y,diag(M(Clist_x,Clist_y,2))];

    %Find weighted average position of traffic bulk, in x and y
    Avpos_vec_x = zeros(1,L);
    Avpos_vec_y = zeros(1,N);
    
    for j = 1:L
        Avpos_vec_x(j) = sum(M(j,2:(N+1),1));
    end
    
    for j = 1:N
        Avpos_vec_y(j) = sum(M(:,j+1,1));
    end
    
    %calculates average position of car bulk in x and y, as well as
    %proportion of cars that moved this step
    Avpos_x(k) = sum(linspace(1,L,L).*Avpos_vec_x)/sum(Avpos_vec_x);
    Avpos_y(k) = sum(linspace(1,N,N).*Avpos_vec_y)/sum(Avpos_vec_y);
    movevec = [movevec move];
    %M(:,:,1)
end

propcars= length(C(:,1))/(L*N)
numcars = length(C(:,1))
%here, -x direction is in opposite direction of movement

%figure(1)
%plot(Avpos_x)
%ylabel('average position (x)')
%xlabel('time')


figure(2)
plot(Avpos_y,'Color','r')
line([0,T],[0.5 0.5],'LineStyle','-')
for i = 1:N
    line([0,T],[i+0.5 i+0.5],'LineStyle',':')
end
line([0,T],[N+0.5 N+0.5],'LineStyle','-')
ylabel('average lane position of cars (y)')
xlabel('time')
ylim([0 N+1])

figure(3)
plot(movevec/length(C(:,1)))
xlabel('time')
ylabel('Number of cars moved at timestep/ Total number of cars')
ylim([-0.05 1.05])

%%%Fast fourier transform
%window function:  .*sin(pi*linspace(0,Avpos_x(end),length(Avpos_x))/Avpos_x(end))
FFT = (fft((Avpos_x-mean(Avpos_x)),2^nextpow2(T))).^2;
f = L*linspace(1,length(FFT),length(FFT))/T;
%Finds positions in FFT corresponding to max and min car speeds, finds max
%of FFT in this range
Avspeed = f(find(FFT ==max(FFT(max(find(f<Smin)):max(find(f<Smax))))))

figure(4)
plot(f,FFT)
xlabel('average speed (spaces/step)')
xlim([0 Smax+2])
ylim([0 1.1*real(max(FFT))])

