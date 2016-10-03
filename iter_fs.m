function [ R ] = iter_fs( data, k,  MAX_ITER, ZERO)
% ITER_FS Chooses a set of k columns from the given data matrix that provide a good rank-k approximation
%
% Author: Bruno Ordozgoiti
% Date: September 28, 2016
%
% If you use this code, please cite the following paper:
%
% Bruno Ordozgoiti, Sandra Gomez Canaval, Alberto Mozo, "A Fast Iterative 
% Algorithm for Improved Unsupervised Feature Selection". IEEE International 
% Conference on Data Mining, 2016.
%
%
% Arguments:
%    data: The input data matrix
%    k: The number of columns to choose
%    MAX_ITER: (Optional) The maximum number of iterations
%    ZERO: (Optional) The threshold below which scores are set to zero, to avoid numerical issues
%   
%   Output:
%    R: A set of column indices indicating the ones that were chosen
  
    if nargin < 4
        ZERO = 1e-5;
    end
    if nargin < 3
        MAX_ITER = Inf;
    end

    M = size(data,1);
    N = size(data,2);

    converged = false;

    previous_norm = 0;
    iterations = 0;

    %Transform the matrix into a more compact form (Theorem 2)
    if M > N
        [~,S,V]=svd(data, 'econ');
        %[~,S,V]=svdecon(data);
        sigma_vt = S*V';
        A = sigma_vt(1:N, :);    
    else
        A = data;
    end

    current_subset = randsample(N, k, false);
    C = A(:,current_subset);    
    cpinv = pinv(C); %TODO Sometimes unstable
    E = A - C*(cpinv*A);    
    %This is in order not to lose the previous subset when converged
    previous_subset = current_subset;
    
    % (Theorem 3)
    if (N > M) 
        [~,SE,VE]=svd(E, 'econ');
        %[~,SE,VE]=svdecon(E);
        compact_E = SE*VE';
        F = SE*compact_E;
        f = sum(F.^2);
        g = sum(compact_E.^2);
    else
        F = E'*E;
        f = sum(F.^2);
        g = diag(F)';
    end

    current_iter = 0;   
    while not(converged) && current_iter < MAX_ITER
        for i = 1:k
            j = current_subset(i);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % Update pseudoinverse
            Ct = A(:,current_subset);        
            Ct(:,i) = 0; 

            cPinvAi = cpinv(i,:)*A;            
            CpCpi = cpinv * cpinv(i,:)';
            nn = CpCpi(i,:);
            cpinv = cpinv - 1/nn * (CpCpi * cpinv(i,:));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Update residual
            CtGA = -1/nn * Ct*(CpCpi*(cPinvAi));
            pre_residual = E;
            E = E + A(:,j) * cPinvAi - CtGA;        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Update scores to account for missing one         
            delta = E(:,j)' * E;
            delta_norm = norm(delta)^2; 
            gamma = pre_residual' * (pre_residual * delta');
            delta_sq = delta.*delta;
            deltaj_sq = delta(j)^2;

            f = f + delta_norm*(delta_sq/deltaj_sq) + 2.0*(gamma' .* (delta/delta(j)));        
            g = g + delta_sq/delta(j);         
            scores = f./g;

            % Set the scores of already present columns to 0
            for h = current_subset'
                if h ~= current_subset(i)
                     scores(h) = 0;
                end
            end  
            % To avoid numerical issues
            scores(find(g < ZERO)) = 0;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Choose the winning feature
            [~,winner] = max(scores);
            current_subset(i) = winner;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Update scores to account for winner
            delta = E(:,winner)' * E;
            delta_norm = norm(delta)^2; 

            gamma = E' * (E * delta');   
            delta_sq = delta.*delta;
            deltaj_sq = delta(winner)^2;
            f = f + delta_norm*(delta_sq/deltaj_sq) - 2.0*(gamma' .* (delta/delta(winner)));
            g = g - delta_sq/delta(winner); 

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Update pseudoinverse
            v = cpinv * A(:,winner);
            w = E(:,winner);
            nw = w' * w;
            third = 1/nw *  w';
            G = -1/nw * v*w';
            G(i,:) = third;
            cpinv = cpinv + G;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            E = E - (E(:,winner) * (E(:,winner)' * E) / (E(:,winner)'*E(:,winner)));
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check convergence     
        if MAX_ITER > 1
            C = A(:,current_subset);
            cnorm = norm(C*(pinv(C)*A), 'fro');
            if cnorm > previous_norm               
                previous_norm = cnorm;
                previous_subset = current_subset;
            else
                converged = true;
                current_subset = previous_subset;
            end
            iterations = iterations + 1;
    
        end
        current_iter = current_iter + 1;
    end    
    R = current_subset;
end

