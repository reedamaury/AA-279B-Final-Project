% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~% Following line commented out by Barrows 1/2014% function [V1, V2] = curtis_lambert(R1, R2, t, string)% Following line added by Barrows 4/2015 (function renamed; outputs and inputs changed)function [v1_out, v2_out, error_out] = AA279lambert_curtis(mu, r1_in, r2_in, dm, nrev, dtsec)% Following 3 lines of comments added by Barrows 1/2014% dm should be 'pro' (prograde) or 'retro' (retrograde)% nrev is ignored (no multi-revolution solutions)% dtsec is length of transfer in [sec]% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%{  This function solves Lambert's problem.  mu         - gravitational parameter (km^3/s^2)  R1, R2     - initial and final position vectors (km)  r1, r2     - magnitudes of R1 and R2  t          - the time of flight from R1 to R2 (a constant) (s)  V1, V2     - initial and final velocity vectors (km/s)  c12        - cross product of R1 into R2  theta      - angle between R1 and R2  string     - 'pro'   if the orbit is prograde               'retro' if the orbit is retrograde  A          - a constant given by Equation 5.35  z          - alpha*x^2, where alpha is the reciprocal of the               semimajor axis and x is the universal anomaly  y(z)       - a function of z given by Equation 5.38  F(z,t)     - a function of the variable z and constant t,             - given by Equation 5.40  dFdz(z)    - the derivative of F(z,t), given by Equation 5.43  ratio      - F/dFdz  tol        - tolerance on precision of convergence  nmax       - maximum number of iterations of Newton's procedure  f, g       - Lagrange coefficients  gdot       - time derivative of g  C(z), S(z) - Stumpff functions  dum        - a dummy variable  User M-functions required: stumpC and stumpS%}% ----------------------------------------------% Following line commented out by Barrows 1/2014 (mu is now an input)%global mu% Following 4 lines added by Barrows 1/2014R1 = r1_in;R2 = r2_in;t = dtsec;string = sprintf('%s', dm);%...Magnitudes of R1 and R2:r1 = norm(R1);r2 = norm(R2);c12   = cross(R1, R2);theta = acos(dot(R1,R2)/r1/r2);% Comment from Barrows 1/2014: For following 'if' statement, note that% number of arguments has changed from original function%...Determine whether the orbit is prograde or retrograde:if nargin < 4 || (~strcmp(string,'retro') & (~strcmp(string,'pro')))    string = 'pro';    fprintf('\n ** Prograde trajectory assumed.\n')endif strcmp(string,'pro')    if c12(3) <= 0        theta = 2*pi - theta;    endelseif strcmp(string,'retro')    if c12(3) >= 0        theta = 2*pi - theta;    endend%...Equation 5.35:A = sin(theta)*sqrt(r1*r2/(1 - cos(theta)));%...Determine approximately where F(z,t) changes sign, and%...use that value of z as the starting value for Equation 5.45:z = -100;while F(z,t) < 0    z = z + 0.1;end%...Set an error tolerance and a limit on the number of iterations:tol   = 1.e-8;nmax  = 5000;%...Iterate on Equation 5.45 until z is determined to within the%...error tolerance:ratio = 1;n     = 0;while (abs(ratio) > tol) & (n <= nmax)    n     = n + 1;    ratio = F(z,t)/dFdz(z);    z     = z - ratio;end%...Report if the maximum number of iterations is exceeded:if n >= nmax    fprintf('\n\n **Number of iterations exceeds %g in ''lambert'' \n\n ',nmax)end%...Equation 5.46a:f    = 1 - y(z)/r1;%...Equation 5.46b:g    = A*sqrt(y(z)/mu);%...Equation 5.46d:gdot = 1 - y(z)/r2;%...Equation 5.28:V1   = 1/g*(R2 - f*R1);%...Equation 5.29:V2   = 1/g*(gdot*R2 - R1);% Following 3 lines added by Barrows 1/2014v1_out = V1;v2_out = V2;error_out = 'finished';return% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~% Subfunctions used in the main body:% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%...Equation 5.38:function dum = y(z)    dum = r1 + r2 + A*(z*S(z) - 1)/sqrt(C(z));end%...Equation 5.40:function dum = F(z,t)    dum = (y(z)/C(z))^1.5*S(z) + A*sqrt(y(z)) - sqrt(mu)*t;end%...Equation 5.43:function dum = dFdz(z)    if z == 0        dum = sqrt(2)/40*y(0)^1.5 + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0)));    else        dum = (y(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z)) ...               + 3*S(z)^2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(y(z)) ...               + A*sqrt(C(z)/y(z)));    endend%...Stumpff functions:function dum = C(z)    dum = stumpC(z);endfunction dum = S(z)    dum = stumpS(z);endend %lambert% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	% Subfunctions stumpS and stumpC added by Barrows 1/2014% They were retrieved from Curtis's textbook's website% http://booksite.elsevier.com/9780080977478% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~function s = stumpS(z)% ~~~~~~~~~~~~~~~~~~~~~~%{This function evaluates the Stumpff function S(z) accordingto Equation 3.52.z - input arguments - value of S(z)User M-functions required: none%}% ----------------------------------------------if z > 0s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;elseif z < 0s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;elses = 1/6;endend% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~function c = stumpC(z)% ~~~~~~~~~~~~~~~~~~~~~~%{This function evaluates the Stumpff function C(z) accordingto Equation 3.53.z - input argumentc - value of C(z)User M-functions required: none%}% ----------------------------------------------if z > 0c = (1 - cos(sqrt(z)))/z;elseif z < 0c = (cosh(sqrt(-z)) - 1)/(-z);elsec = 1/2;endend% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~