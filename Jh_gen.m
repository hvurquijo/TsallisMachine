function[varargout]=Jh_gen(varargin)
n = varargin{1};
opt = varargin{2};
J = zeros(n);
h = zeros(1,n);
if nargin==3
    opt2 = varargin{3};
    if opt2 == 1
        if opt == 1 % J,h in [-1,1]
            for i = 1:n-1
                for j = i+1:n
                    J(i,j) = normrnd(0,0.25);
                end
            end

            h = rand(1,n);
            h = 2*h-1;

        elseif opt ==3 %J,h in [0,2]

            for i = 1:n-1
                for j = i+1:n
                    J(i,j) = normrnd(1,0.25);
                end
            end
            %J = (J+1)/2;

            h = rand(1,n);
            h = 2*h-1;

        elseif opt ==4 %J,h in [-2,0]

            for i = 1:n-1
                for j = i+1:n
                    J(i,j) = normrnd(-1,0.25);
                end
            end
            %J = (J-1)/2;

            h = rand(1,n);
            h = 2*h-1;

        end
    end
    
else
    
    if opt == 1 % J,h in [-1,1]
    
        J = randi([-1000,1000],n);
        J=J/max(max(abs(J)));
        J = triu(J,1);

        h = rand(1,n);
        h = 2*h-1;
        
    elseif opt ==2 % J,h=(-1,1)

        J = randi([0,1],n);
        J = (2*J-1);
        J = triu(J,1);

        h = rand(1,n);

    elseif opt ==3 %J,h in [0,1]

        J = rand(n);
        J = triu(J,1);

        h = rand(1,n);
        h = 2*h-1;

    elseif opt ==4 %J,h in [-1,0]

        J = -rand(n);
        J = triu(J,1);

        h = rand(1,n);

    end
end




varargout{1} = J;
varargout{2} = h;