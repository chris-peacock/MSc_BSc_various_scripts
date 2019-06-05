clc;
clear;
%%%Just for fun; generates a technicolour mesh of the Mandelbrot set;
%%%Points 'c' in the real/complex plane are included in the set if the function
%%%f_c = y^2 + c does not diverge when iterated from z = 0, for which the
%%%sequence f_c(0), f_c(f_c(0)),... remains bounded in absolute value
%%%-quoted from Wikipedia

%A (max x max) grid of points in the complex plane, in the range 
%(xlimlow->xlimhi, ylimlow->ylimhi)
max = 200
xlimlo= -2;
xlimhi= 2;
ylimlo=-2;
ylimhi=2;
mat = [];
x = linspace(xlimlo,xlimhi,max);
y = i.*linspace(ylimlo,ylimhi,max);
t_0 = cputime;
for r= 1:max
  for j = 1:max
    mat(r,j) = x(r)+y(j);
    x_q = 0;
    y_q = 0;
    for q = 1:100
      z_0 = 0;
      if x_q <100
        z = x_q^2 + 2*i*x_q*y_q - y_q^2 + mat(r,j);
        x_q = real(z);
        y_q = imag(z);
        z_0 = x_q + i*y_q;
      else 
        mat(r,j) = NaN;
      end
    end
  end
end
%Print time elapsed in calculation
t = cputime-t_0
figure(1)
mesh(real(mat));


