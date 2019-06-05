clf;
clc;
clear;
%%%Monte Carlo simulation of a resin curing in 2D.  Initial state of system\
%%%consists of a grid of monomers and holes; monomers can jump to adjacent
%%%holes with some probability (Gamma_0); if adjacent to another monomer 
%%%polymerize with probability (Gamma_p), polymer is immobile.  System 
%%%surrounded by 'boundary' units; if hole comes in contact with boundary,
%%%immediately becomes boundary unit itself.  So, system of holes and 
%%%monomers slowly shrinks and becomes more rigid (polymerized) with time%%%

%%%Fills matrix NxN with random holes, monomers.  Border around NxN matrix%%%
%Monomers are classified with 1
%Polymers are classified with 0
%Holes are classified with -1
%Boundaries are classified with -2

%4 choices/step: polymerize, hop horizontal, up or down. 'gam' establishes 
%preferential up/down

%List 'I' catalogues monomer position and neighbours; has six columns:
%col 1,2 -> x, y pos of monomer
%col 3-6 -> left, right, up down neighbour identity

%List 'I_hole' catalogues hole position and neighbours; has six columns:
%col 1,2 -> x, y pos of hole
%col 3-6 -> left, right, up down neighbour identity
tic
%Path to save text files to (placeholder here)
Path = "path1\path2\path3\"
%solution size (NxN square)
N = 20;
%porportion of initial NxN system composed of monomers
monomweight = 0.5;
%polymerization weighting
Gam_p = 1;
%Hop weighting
Gam_0 = 1;
%Hop direction weighting
gam = 2;

M = ones(N+2,N+2);
M(:,1) = -2;
M(:,N+2) = -2;
M(1,:) = -2;
M(N+2,:) = -2;
monom = zeros(N,N);
monom(1:round(monomweight*N^2)) = 1;
M(2:N+1,2:N+1) = monom(randi([1,(length(monom)^2)],N,N))*2 - 1;
t = 0;

n_p= [];
n_m = [];
n_h= [];
n_edge = [];
tvec = [];
t = 0;
%dlmwrite([Path,'startM.txt'],M);

while (t<50  && length(find(M==1))>0)
  [row_m,col_m] = find(M==1);
  I = zeros(length(row_m),6);
  for i = 1:length(row_m)
    x_i = row_m(i);
    y_i = col_m(i);
    %fills list of index for monomer elements, with left,right,up,down of element entries
    I(i,:) = [x_i,y_i,M(x_i,y_i-1),M(x_i,y_i+1),M(x_i-1,y_i),M(x_i+1,y_i)];
  end

  N_1 = length(find(I(:,3:6) ==1));
  %counts number of holes adjacent left or right to each monomer
  N_2 = length(find(I(:,3:4) ==-1));
  %counts number of holes adjacent above each monomer
  N_3 = length(find(I(:,5) ==-1));
  %counts number of holes adjacent below each monomer
  N_4 = length(find(I(:,6) ==-1));
  Gam = [Gam_p*N_1,Gam_0*N_2,Gam_0*N_3/gam,Gam_0*N_4*gam];
  Gam_tot = sum(Gam);
  %prob of polym, move horizontal, move up, move down
  P = Gam./Gam_tot;
  
  %chooses index 1:4 (polym, left or right, up, down) based on probabilities of each
  %cumsum makes 'bins' of width of specific probability, rand falls into one of the bins
  %bins are integer indexed, 1+ bin index = choice per step
  Ch_dir = (1+sum(rand>cumsum(P)));
  
  if Ch_dir ==1
  %%%POLYMERIZE%%%
    %Chooses monomer in list I and direction of neighbouring monomer to polymerize
    %find(I=0) allows monomers to polymerize with existing polymers as well
    [Ch_p_row,Ch_p_col] = find((I(:,3:6)==1)); % | I(:,3:6)==0
    if length(Ch_p_row)>0
      %Chooses a random indexed monomer in I to polymerize
      Ch_p_index = randi(length(Ch_p_row));
      Ch_p_array = [Ch_p_row,Ch_p_col]; 
      Ch_p = Ch_p_array(Ch_p_index,:);
      %sets chosen monomer to 0 (polymerized)
      M(I(Ch_p(1),1),I(Ch_p(1),2)) = 0;
      %sets neighbour to 0 (ploymerized)
      if Ch_p(2) == 1
        M(I(Ch_p(1),1),I(Ch_p(1),2)-1) = 0;
      elseif Ch_p(2) == 2
        M(I(Ch_p(1),1),I(Ch_p(1),2)+1) = 0;
      elseif Ch_p(2) == 3
        M(I(Ch_p(1),1)-1,I(Ch_p(1),2)) = 0;
      elseif Ch_p(2) == 4
        if (M(I(Ch_p(1),1)+1,I(Ch_p(1),2))== 1)

          originalpos = I(Ch_p(1),1:2);
          neighbourrow1 = find(I(:,1)== (I(Ch_p(1),1)+1));
          neighbourrow2 = find(I(neighbourrow1,2)== I(Ch_p(1),2));
          neighbour_pos = I(neighbourrow1,(1:2));
        end
        M(I(Ch_p(1),1)+1,I(Ch_p(1),2)) = 0;
        
      end 
    end
  elseif Ch_dir ==2
  %%%HOP HORIZ%%%
  %disp('horizontal')
    [Ch_h_row,Ch_h_col] = find((I(:,3:4)==-1));
    if length(Ch_h_row)>0
      Ch_h_index = randi(length(Ch_h_row));
      Ch_h_array = [Ch_h_row,Ch_h_col];
      Ch_h = Ch_h_array(Ch_h_index,:);
      if Ch_h(2)==1
        M(I(Ch_h(1),1),I(Ch_h(1),2)) = -1;
        M(I(Ch_h(1),1),I(Ch_h(1),2)-1) = 1;
      elseif Ch_h(2)==2
        M(I(Ch_h(1),1),I(Ch_h(1),2)) = -1;
        M(I(Ch_h(1),1),I(Ch_h(1),2)+1) = 1;
      end
    end
    %M
  elseif Ch_dir ==3
  %%%HOP UP%%%
  %disp('up')
    Ch_h_row = find((I(:,5)==-1));
    Ch_h = Ch_h_row(randi(length(Ch_h_row)));
    if length(Ch_h_row)>0
      M(I(Ch_h,1),I(Ch_h,2)) = -1;
      M(I(Ch_h,1)-1,I(Ch_h,2)) = 1;
    end
    %M
  elseif Ch_dir ==4
  %%%HOP DOWN%%%
  %disp('down')
    Ch_h_row = find((I(:,6)==-1));
    Ch_h = Ch_h_row(randi(length(Ch_h_row)));
    if length(Ch_h_row)>0
      M(I(Ch_h,1),I(Ch_h,2)) = -1;
      M(I(Ch_h,1)+1,I(Ch_h,2)) = 1;
    end
    %M
  end
  
  %indexes all the current holes
  [row_hole,col_hole] = find(M==-1);
  I_hole = zeros(length(row_hole),6);
  for i = 1:length(row_hole)
    x_i = row_hole(i);
    y_i = col_hole(i);
    %fills list of index for monomer elements, with left,right,up,down of element entries
    I_hole(i,:) = [x_i,y_i,M(x_i,y_i-1),M(x_i,y_i+1),M(x_i-1,y_i),M(x_i+1,y_i)];
  end
  
  %sets any holes adjacent to boundary to boundary themselves
  [Ch_hole_row,Ch_hole_col] = find((I_hole(:,3:6)==-2));
  if length(Ch_hole_row)>0
    for i = 1:length(Ch_hole_row)
      M(I_hole(Ch_hole_row(i),1),I_hole(Ch_hole_row(i),2)) = -2;
    end
  end
  %M
  
  n_p = [n_p length(find(M==0))];
  n_m = [n_m length(find(M==1))];
  n_h = [n_h length(find(M==-1))];
  n_edge = [n_edge ((N+2)*(N+2) - length(find(M==-1)) - length(find(M==0)) - length(find(M==1)))];
  tvec = [tvec t];
  t = t+(-((Gam_tot)^(-1))*log((1-0.1)*rand(1,1) + 0.1));
end

%dlmwrite([Path,'endM.txt'],M);

figure(1)
plot(log10(tvec),n_m/((N+2)*(N+2)),'Color','b')
xlabel('log(time)')
ylabel('Proportion of system composed of monomers')
title(['N = ' num2str(N) ',  \Gamma_p = ' num2str(Gam_p) ', \Gamma_0 = ' num2str(Gam_0) ', \gamma = ' num2str(gam)])
ylim([0 1])
xlim([log10(tvec(1)) log10(tvec(end))])
hold on;

figure(2)
plot(log10(tvec),(((N^2)- n_edge)/(N^2)),'Color','b')
xlabel('log(time)')
ylabel('Proportion of system volume remaining')
title(['N = ' num2str(N) ',  \Gamma_p = ' num2str(Gam_p) ', \Gamma_0 = ' num2str(Gam_0) ', \gamma = ' num2str(gam)])
ylim([(((N)*(N)- n_edge(end))/((N+2)*(N+2))) 1])
xlim([log10(tvec(1)) log10(tvec(end))])
hold on;

figure(3)
plot(log10(tvec(2:length(tvec))),-diff(((N+2)*(N+2)- n_edge)/((N+2)*(N+2))))
xlabel('log(time)')
ylabel('(Differential of) Proportion of system volume remaining')
title(['N = ' num2str(N) ',  \Gamma_p = ' num2str(Gam_p) ', \Gamma_0 = ' num2str(Gam_0) ', \gamma = ' num2str(gam)])
%ylim([0 1])
xlim([log10(tvec(1)) log10(tvec(end))])

%figure(4)
%plot(log10(tvec(2:length(tvec))),diff(n_m+n_p))



toc