function [H,pValue] = kupiertest(x,F,alpha)
%KUPIERTEST Single Sample Kupier Goodness-Of-Fit Hypothesis Test.  
%   Similar to Kolmogorov-Smirnov (K-S) test, but (K-S) test tend
%   to be most sensitive around median value of the distribution
%   and less sensitive at the distribution tails. 
%   Doesn't require any Toolbox.
%   
%   Performs a Kupier test to determine if the random sample X could have 
%   the hypothesized continuous cumulative distribution function F. 
%   Null Hypothesis Ho: "The Sample is taken form a population with cumulative 
%   distribution F(x) for all x". 
%   This is a 2-sided test. Alternative: "The population from which the sample is taken 
%   has a cumulative distribution different form F "
%
%   Test statistics is K = max(S(x) - F(x)) + max(F(x) - S(x))
%   where S(x) is the empirical cumulative distribution function
%
%   Calculation are based on the asynthotic approximation of pValue
%   presented in:
%   %  Press, W.H., et. al., "Numerical Recipes in C", 
%         Cambridge University Press, 1992.
%
%   OUTPUT: 
%   H is the Boolean indicating the result of the hypothesis test:
%      H = 0 => Do not reject the null hypothesis at significance level ALPHA.
%      H = 1 => Reject the null hypothesis at significance level ALPHA.
%   PVALUE is the associated pValue
%
%   INPUT:
%   X must be a row or column vector representing a random sample from some
%   underlying distribution. NaN's observations in X are ignored. X don't need 
%   to be sorted. 
%   F is the c.d.f. under the null hypothesis. If specified, it must be a
%   raw or column vector with the same number of elements of X. It is intended 
%   to be the vector of values of F(x) for each x in the samples (in the
%   same order). If F omitted or set to an empty matrix [], the hypothetical
%   c.d.f is assumed to be a standard normal, N(0,1).    
%
%   For Example:
%   randn('state',0)
%   x=randn(10000,1);
%   [H,pValue] = kupiertest(x)
%   --> H =0
%   --> pValue=0.4051
%   
%  Author:  Daniele Mossi, 13th of November 2005
%           Nextra Investment Management SGR
% 

%%%%% Work with X
if nargin<1 | min(size(x))>1
    error('X must be a row or column vector!');
end

x=x(:);  % always column vector
nOR=length(x); % original Number of elements in X
NotNan=~isnan(x); % Logical Indices of valid elements 
if sum(NotNan)==0
    error('X must have some valid elements !');
end

x=x(NotNan);    % Here we take only valid elements
[x,i]=sort(x);  % And sort them
n=length(x);    % Number of Valid elements in X

% Calculates the sample cdf S(x), (n points without the initial zero)
S= (1:n)'/n;   % Once X has been sorted, these are S(x) values

% Here we look for duplicated Values as in "cdfcalc" function of
% Statistics Toolbox
NotZeroDiff = ([diff(x);1] > 0); % Look for Duplicated values in X
x=x(NotZeroDiff);   % Remove them from X 
S=S(NotZeroDiff);   % But in S duplicated values are weighted properly
S=[0;S];             % Now we put the initial zero (for the left/right approach)

%% For Example x=[1, 2, 3, 4] ---> S=[0, 0.25, 0.50, 0.75, 1]
%% If instead  x=[1, 2, 2, 3, 4]-> S=[0, 0.20, 0.60, 0.80, 1]

%%%%% Work with F
if nargin>=2 & ~isempty(F)
    if min(size(F))>1 | length(F)~=nOR;
        error('F must be a row or column vector with the same number of elements as X');
    end    

    F=F(:);   % Always columns vector
    F=F(NotNan); % Here we take only elements correspondent to valid elements in X 
    if any(isnan(F)) % But these n values must be all valid !
        error('Some elements of F are not valid numbers');
    end

    F=F(i);   % F data are placed in the same order as those in X
    F=F(NotZeroDiff);   % Remove elements corrisponding to duplicated valus in X
    if any(diff(F))<=0  % This is when F is not a strictly incresing function of x
       error('Data Provided for F are not those of a continuos distribution');
    end    
else  % When input F is not provided
   F=normcdf(x,0,1);  % Standard Normal cdf
end

% ALPHA, if provided must be a scalar between 0 and 1
if (nargin >= 3) & ~isempty(alpha)
   if prod(size(alpha)) > 1 | (alpha <= 0 | alpha >= 1)
      error('Alpha must be a scalar in (0,1).');
   end
else
   alpha  =  0.05;  % Deafult 5% first type error
end

% max( S(x) - F(x) )
delta1=S(1:end-1)-F;   % difference from the LEFT.
delta2=S(2:end)-F;   % difference from the RIGHT.
deltaF=[delta1 ; delta2];  % S(x) - F(x) vector
Kplus=max(deltaF);  

% max ( F(x) - S(x) )
delta1=-delta1;   % difference from the LEFT.
delta2=-delta2;   % difference from the RIGHT.
deltaF=[delta1 ; delta2]; % F(x) - S(x) vector
Kminus= max(deltaF);

KupierStat=Kplus+Kminus;

% Now we compute the Asyntotic p-Value 
lambda  =  max((sqrt(n) + 0.155 + 0.24/sqrt(n)) * KupierStat,0); % This max is useless if F data in input are correct...
if lambda<0.4  % Useless to compute pValue (pValue equals 1 at the 7th decimal
               % For small value of lambda some problems may also arise with sum convergence)
    pValue=1;  % KupierStat very small for this sample length: never reject Ho
    H=0;
    return
end

j=(1:100)';  % 100 addends are abundant for converge 
             % I think this is preferable instead of 
             % testing for convergence condition in a for-loop   
pare=4*lambda*lambda*(j.^2)-1;
expo=exp(-2*lambda*lambda*(j.^2));
argo=pare.*expo;
pValue=2*sum(argo);

if pValue < 0 , pValue = 0; end
if pValue > 1 , pValue = 1; end

if pValue<alpha
    H=1;  % Reject Ho
else
    H=0;  % Accept Ho
end

%%% END OF FUNCTION

