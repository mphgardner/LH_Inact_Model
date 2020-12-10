function results = linearTDRL_Lambda_LH_Inact_Mac(X,r,varargin)
    
    % TD lambda RL algorithm employing linear function approximation. 
    % An additional Mackintosh (1975) salience mechanism is included
    % for adjusting elemental associability (independent alphas for each
    % of the features are adjusted based on how predictive a feature is of
    % reward)
    %
    %
    % USAGE: results = linearTDRL_Lambda_Inact_Mac(X,r,varargin)
    %
    % INPUTS:
    %   X - [N x D] stimulus sequence (N = # trials, D = # features)
    %   r - [N x 1] reward sequence
    %   V (optional) - initial values of the features; vector (D X 1) (default: all zeros)
    %   lambda (optional) - persistence of eligibility trace; scalar (default: zero)
    %   gamma (optional) - temporal discounting factor; scalar (default: 0.9)
    %   alpha (optional) - initial learning rate; scalar (default: 0.05)
    %   opto (optional) - sequence of optogenetic perturbations (default: all zeros)
    %   opto_update (optional) - inactivation of the eligibility trace
    %       'none'     : no inactivation (default)
    %       'update'     : disrupts the establishment of an eligibity trace (0: no effect)
    %   theta (optional) - the magnitude of the incremental update to alpha 
    %                       according to Mackintosh 1975
    %   Mac_Trials (optional) - states at which to update alpha (default: zero)
    %  
    %
    % OUTPUTS:
    %   results - [N x 1] structure with the following fields:
    %               .V - value estimate
    %               .dt - prediction error
    %               .W - weight matrix
    
    %
    % Matt Gardner, 2020
    %See Seijen and Sutton 2014 for TD Lambda with linear function approximation
    %Currently this model uses replacing eligibility traces 
    %The error is normalized by the factor (1-lambda)/lambda). This
    %normalization results in an equivalent total error if all features are
    %present at every time step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [N,D] = size(X);
    
    %Default settings:
    W = zeros(D,1); lambda = 0; alpha = 0.05*ones(N+1,D); opto = zeros(N+1,1); opto_update = 'none'; gamma = 0.9; theta = 0; mac_trials = 0;
       
    %This parses inputs
    if ~isempty(varargin)
        
        for i = 1:2:nargin-2
            
            type = varargin{i};
            
            switch type
                
                case 'Vinitial'
                    
                    if numel(varargin{i + 1}) ~= D
                        fprintf('Vinitial should be a vector with the length of the number of columns of X \n')
                        return
                    else
                       W(:,1) = varargin{i + 1}; 
                    end
                case 'lambda'
                    lambda = varargin{i + 1};
                
                case 'gamma'
                     gamma = varargin{i + 1};   
                    
                case 'alpha'
                    if numel(varargin{i + 1}) == 1
                        alpha = varargin{i + 1}.*ones(N+1,D);
                    
                    elseif all(size(varargin{i + 1}) == [1 D])
                        alpha = repmat(varargin{i + 1},N+1,1);
                        
                    elseif numel(varargin{i + 1}) == N
                        alpha = [alpha(1,:); varargin{i + 1}];
                    end    
                
                case 'opto' 
                    opto = [0; varargin{i + 1}]; 
                
                case 'opto_update'
                    opto_update = varargin{i + 1};
                case 'theta'
                    theta = varargin{i + 1};
                case 'Mac_Trials'
                    mac_trials = varargin{i + 1};    
                otherwise
                    fprintf('%s%s\n','Unrecognized input: ', varargin{i})
                    return
            end
        end
    end     
    
    %Initialize the eligibility trace vector
    E = zeros(N+1,1);
    if lambda > 0
        E(1:2) = lambda*gamma*1;
    else
        E(2) = gamma*1;
    end    
    r = [0;r];
    X = [zeros(1,size(X,2));X];
    %W = WM(1,:)';
    norm = (1 - lambda)/lambda;
    norm(norm == Inf) = 1;
    
    %preallocate structure:
    results.dt = zeros(N,1);
    results.V = zeros(N,1);
    results.W = zeros(N,D);
    
    for n = 2:N %The vector has been padded to add an extra state at the beginning
        
        dt = r(n+1) + (gamma*X(n+1,:) - X(n,:))*W;

        % store results
        results.dt(n-1,:) = dt;
        results.V(n-1,:) = X(n,:)*W; 
        results.W(n-1,:) = W;
        
        % update using an eligibility trace and opto manipulation
        switch opto_update

            case 'update' %in the case of TD(0) this is equivalent to a gain modulation of the
                        %error term
            W = W + norm*sum((1 + opto(1:n)).*alpha(1:n,:).*E(1:n).*X(1:n,:)*dt)';
            
            case 'none'

            W = W + norm*sum(alpha(1:n,:).*E(1:n).*X(1:n,:)*dt)';

        end
        
        %update the eligibility trace 
        if lambda > 0
            E(n+1) = 1;
            E = lambda*gamma*E;
        else
            E = zeros(size(E));  
        end  
        
        %This implements a change in the salience of the features according to 
        %Mackintosh 1975
        if theta > 0 && rem(n,mac_trials) == 0
            alpha(n+1:end,W == min(W)) = alpha(n+1:end,W == min(W)) + theta;
            alpha(n+1:end,W > min(W)) = alpha(n+1:end,W > min(W)) - theta;
            alpha(alpha > 1) = 1;
            alpha(alpha < 0) = 0;
        end    
            
    end

    
end
  

