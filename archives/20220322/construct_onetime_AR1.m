function path = construct_onetime_AR1(persistence, shockSize, length)
% constructOneTimeAR1 ... This m-file constructs one-time multiplicative AR(1) shock.
% The equation is given by 
%          y(t) = 1 + x(t)
% where x(t+1) = rho * x(t) and x(1) is the size of initial shock.
%
% Syntax:  
%   path = constructOneTimeAR1(persistence, shockSize, length)
%
% Inputs:
%    persistence (positive real number) 
%       - AR(1) prameter
%    shockSize   (real number)
%       - size of initial shocks where no shock = 0
%    length (positive integer)     
%       - length of the mean-reverting process
%
% Outputs:
%    path  (array of real numbers)      
%       - the mean-reverting path of the multiplicative factors
%
% Example: 
%    constructOnetimeAR1(0.5, 0.1, 4) 
%       = [1.1; 1.05; 1.025; 1.0125]
%
%    constructOnetimeAR1(0.5, -0.1, 4)  
%       = [0.9; 0.95; 0.975; 0.9875]
%
% No other files are required.

% Return errors if inputs are not compatible
if ~isreal(persistence) | persistence < 0 | ~isnumeric(persistence)
    error("MATLAB:NotPositiveReal","Persistence should be a positive real number")
end
if ~isreal(shockSize) | ~isnumeric(shockSize)
    error("MATLAB:NotReal","The size of shock should be a real number")
end
if ~(mod(length,1) == 0) | length < 0 | ~isnumeric(length)
   error("MATLAB:NotPositiveInteger", "The length should be a positive integer") 
end


path = zeros(length,1);

temporaryValue  =   shockSize;

for t = 1 : length
    path(t)         =   temporaryValue;
    temporaryValue  =   temporaryValue * persistence;
end

path = path + 1;
