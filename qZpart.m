function[varargout]=qZpart(varargin)
J = varargin{1};
h = varargin{2};
n = varargin{3};
q = varargin{4};
dist_type = varargin{5};
ham_type = varargin{6};

H_tmp = 0;
nstates = 2^n;
si = (dec2bin(0:nstates-1,n)-'0')*2-1; %[-1,1]
Z_mat = zeros(1,nstates);
Z = 0;

parfor i=1:nstates
    H_tmp = ising(J,h,si(i,:),ham_type,q);
    Z_mat(i) = qexp(q,-H_tmp,dist_type);
end
Z = sum(Z_mat);
varargout{1}=Z;

