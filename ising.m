function[varargout]=ising(varargin)
J = varargin{1};
h = varargin{2};
si = varargin{3};%[-1,1]
type = varargin{4};
if nargin==5
 q = varargin{5};
end

sij = sij_gen(si); %[-1,1]

switch type
    case '++'
        H = sum(h.*si)+sum(sum(J.*sij));
    case '--'
        H = -sum(h.*si)-sum(sum(J.*sij));%classic
    case 's+'%it works for J,h = [0,1] otherwise its equivalent to 's-'
        H = -sum(h.*(si-1)/2)-sum(sum(J.*(sij-1)/2));
    case 's-'%it works for J,h = [0,1] otherwise its equivalent to 's+'
        H = -sum(h.*(si+1)/2)-sum(sum(J.*(sij+1)/2));
    case '+-'
        if (~isempty(q))
            H = -sum(h.*(si-sign((1-q)*h))) - sum(sum(J.*(sij-sign((1-q)*J))));
%             if q==1
%                 H =2*H;
%             end
        else
            disp("You must provide the \'q \' value for use this option")
        end
    case '-+'
        if (~isempty(q))
            H = -sum(h.*(si-sign((1-q)*h))) - sum(sum(J.*(sij-sign((1-q)*J))));
%             if q==1
%                 H =2*H;
%             end
        else
            disp("You must provide the \'q \' value for use this option")
        end
    case '+-c'
        if (~isempty(q))
            H = -sum(h.*(si-sign((q-1)*h))) - sum(sum(J.*(sij-sign((q-1)*J))));
        else
            disp("You must provide the \'q \' value for use this option")
        end
    otherwise
        H=0;
end
varargout{1}=H;

        