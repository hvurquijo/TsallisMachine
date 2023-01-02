function[varargout]=qexp(varargin)
q = varargin{1}; %valor de q
x = varargin{2}; %valor de x matriz ou não
if nargin == 3
    opt = varargin{3}; %tipo de qexponencial
else
    opt = 1;
end

if q == 1
    varargout{1} = exp(x);
    return;
end
if  opt ==2
    
    varargout{1}=(1+(1-q)*x).^(1/(1-q));
    
elseif opt ==1
    
    varargout{1}=(1+(q-1)*x).^(1/(q-1));
    
end