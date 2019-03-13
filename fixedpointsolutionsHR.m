% *** FORWORD
% The quantities xbar, ybar and zbar represent the values of the fixed
% points (points for which an ODE is zero) of the Hindmarsh Rose Equations.
% Note that these solutions do not depend on r, since this is cancelled
% during simplification.

a = 1;
b = 3;
c = 1;
d = 5;
I = 0;
x_1 = -1.6;
s = 4;

% The following is the real root of the cubic for xbar
% There are also two complex roots
xbar_real = ((sqrt((27*(a^2)*c + 27*(a^2)*x_1*s + 27*(a^2)*I - 9*a*b*s + 9*a*d*s + 2*(b^3) - 6*(b^2)*d + 6*b*(d^2) - 2*(d^3))^2 + 4*(3*a*s - (d - b)^2)^3)  + 27*(a^2)*c + 27*(a^2)*x_1*s + 27*(a^2)*I - 9*a*b*s + 9*a*d*s + 2*(b^3) - 6*(b^2)*d + 6*b*(d^2) - 2*(d^3))^(1/3))/(3*2^(1/3)*a) - (2^(1/3)*(3*a*s - (d - b)^2))/(3*a*(sqrt((27*a^2*c + 27*a^2*x_1*s + 27*a^2*I - 9*a*b*s + 9*a*d*s + 2*b^3 - 6*b^2*d + 6*b*d^2 - 2*d^3)^2 + 4*(3*a*s - (d - b)^2)^3) + 27*a^2*c + 27*a^2*x_1*s + 27*a^2*I - 9*a*b*s + 9*a*d*s + 2*b^3 - 6*b^2*d + 6*b*d^2 - 2*d^3)^(1/3)) - (d - b)/(3*a);
xbar_complex1 = -((1 - i*sqrt(3))*(sqrt((27*(a^2)*c + 27*(a^2)*x_1*s + 27*(a^2)*I - 9*a*b*s + 9*a*d*s + 2*(b^3) - 6*(b^2)*d + 6*b*(d^2) - 2*(d^3))^2 + 4*(3*a*s - (d - b)^2)^3)  + 27*(a^2)*c + 27*(a^2)*x_1*s + 27*(a^2)*I - 9*a*b*s + 9*a*d*s + 2*(b^3) - 6*(b^2)*d + 6*b*(d^2) - 2*(d^3))^(1/3))/(6*2^(1/3)*a) + ((1 + i*sqrt(3))*(3*a*s - (d - b)^2))/(3*2^(1/3)*a*(sqrt((27*a^2*c + 27*a^2*x_1*s + 27*a^2*I - 9*a*b*s + 9*a*d*s + 2*b^3 - 6*b^2*d + 6*b*d^2 - 2*d^3)^2 + 4*(3*a*s - (d - b)^2)^3) + 27*a^2*c + 27*a^2*x_1*s + 27*a^2*I - 9*a*b*s + 9*a*d*s + 2*b^3 - 6*b^2*d + 6*b*d^2 - 2*d^3)^(1/3)) - (d - b)/(3*a);
xbar_complex2 = -((1 + i*sqrt(3))*(sqrt((27*(a^2)*c + 27*(a^2)*x_1*s + 27*(a^2)*I - 9*a*b*s + 9*a*d*s + 2*(b^3) - 6*(b^2)*d + 6*b*(d^2) - 2*(d^3))^2 + 4*(3*a*s - (d - b)^2)^3)  + 27*(a^2)*c + 27*(a^2)*x_1*s + 27*(a^2)*I - 9*a*b*s + 9*a*d*s + 2*(b^3) - 6*(b^2)*d + 6*b*(d^2) - 2*(d^3))^(1/3))/(6*2^(1/3)*a) + ((1 - i*sqrt(3))*(3*a*s - (d - b)^2))/(3*2^(1/3)*a*(sqrt((27*a^2*c + 27*a^2*x_1*s + 27*a^2*I - 9*a*b*s + 9*a*d*s + 2*b^3 - 6*b^2*d + 6*b*d^2 - 2*d^3)^2 + 4*(3*a*s - (d - b)^2)^3) + 27*a^2*c + 27*a^2*x_1*s + 27*a^2*I - 9*a*b*s + 9*a*d*s + 2*b^3 - 6*b^2*d + 6*b*d^2 - 2*d^3)^(1/3)) - (d - b)/(3*a);

% The corresponding y and z coordinates
ybar = c - d*xbar_real^2;
ybar_complex1 = c - d*xbar_complex1^2;
ybar_complex2 = c - d*xbar_complex2^2;

zbar = s*(xbar_real-x_1);
zbar_complex1 = s*(xbar_complex1-x_1);
zbar_complex2 = s*(xbar_complex2-x_1);

fprintf("xbar = %2.4d, %2.4d, %2.4d\nybar = %2.4d, %2.4d, %2.4d\nzbar = %2.4d, %2.4d, %2.4d\n", xbar_real, xbar_complex1, xbar_complex2, ybar, ybar_complex1, ybar_complex2, zbar, zbar_complex1, zbar_complex2);