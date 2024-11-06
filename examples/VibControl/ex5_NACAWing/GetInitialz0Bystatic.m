function z0 = GetInitialz0Bystatic(DS,tol,max_iter)

    % use Gaussian Elimination to get inverse of K
    u  = [10000; 2000];
   
    G  = DS.K^-1;
    % G  = zeros(size(DS.K)); % inverse of K, called flexibility
    % for i = 1:size(DS.K,2)
    %     ei     = zeros(size(DS.K,1),1);
    %     ei(i,1)= 1;
    %     G(:,i) = Gaussian_Elim(DS.K,ei,max_iter);
    % end
    x0  = DS.epsilon*G*DS.D*u; % initial of Newton's method
    xk  = Inf(size(x0));
    xk1 = x0;
    dxk = zeros(size(xk));
    % tol = 0.0001;
    % use Newton's method to calculate z0
    i      = 1;
    while (max(abs(xk-xk1)) > tol) && (i<=max_iter)
        xk     =  xk1;
        % X(K+1)=X(K)-g(X(K))/g'(X(K+1))
        gxk    =  DS.K*xk + DS.compute_fnl(xk,dxk) - DS.epsilon*DS.D*u;
        dfnl   =  DS.compute_dfnldx(xk,dxk);         % built-in method
        dgxkdx =  DS.K + dfnl;
        % g'_k * s_k = g_k
        % sk     =  Gram_Schmidt(dgxkdx,gxk);
        sk     = dgxkdx^-1*gxk;

        xk1   =  xk - sk;  % hard to calculate large matrix inverse, use Gaussian method
        i       = i + 1;
    end
    disp("iteration number")
    disp(i)
    disp("max abs error=")
    disp(max(abs(xk-xk1)))
    z0 = xk1;
end

function x = Gaussian_Elim(A,b)
    n = size(A,1);
    L_save = cell(n-1,1);    % Save L{i}
    M_save = cell(n-1,1);
    P_save = cell(n-1,1);    % Save the permutation matrix, P{i}
    % Outside loop over colums
    for k = 1 : n-1
        %% pivoting
        % Find the pivot entry in current colums, searching rows by rows
        p = k;  % the index of the row with max diagonal entry
        P = eye(n,n);
        M = eye(n,n); 
        L = eye(n,n);
        for i = k+1 : n
            if A(p,k) < A(i,k)
                % Find larger entry, change the index
                p = i;
            end
        end
        if p ~= k
            % interchange
            A([k p],:) = A([p k],:);
            P([k p],:) = P([p k],:);
        end
            P_save{k} = P;
        %% if diagonal entry equals to 0, next row
        if A(k,k) == 0
            continue;
        end
    %% Gaussian Elimination
        % Update M matrix, finally get L matrix
        for i = k+1 : n
            L(i,k) = A(i,k)/A(k,k);
            M(i,k) = -A(i,k)/A(k,k);
        end
    
        % Elimination
        for j = k : n
            for i = k+1 : n
                A(i,j) = A(i,j) + M(i,k)*A(k,j);
            end
        end
        
        %% Gaussian Elimination, use vector operation 
        % M(k+1 : n,k) = -A(k+1 : n,k)./A(k,k);
        % L(k+1 : n,k) = A(k+1 : n,k)./A(k,k);
        % A(k+1 : n, k : n) = A(k+1 : n, k : n) + M(k+1 : n,k) * A(k,k : n);
    
        L_save{k} = L;
        M_save{k} = M;
        b = M*P*b;
    end
    %% Obtain L and U
    U = A;
    L = eye(n,n);
    for i = 1:n-1
        L =  L * L_save{i};
    end
    % L = eye(n,n);
    % for i = 1:n-1
    %     L = L * P_save{i}' * L_save{i};
    % end
    
    %% Forward and Backward Substitution
    % % forward
    % y = zeros(n,1);
    % for i = 1 : n
    %     if L(i,i)==0
    %         break
    %     end
    %     y(i) = b(i);
    %     for j = 1:i-1
    %         y(i) = y(i) - L(i,j) * y(j);  
    %     end
    %     y(i) = y(i) / L(i,i);
    % end
    % 
    
    
    % backward
    x = zeros(n,1);
    for i = n : -1 : 1
        if U(i,i)==0
            break
        end
        % x(i) = y(i);
        x(i) = b(i);
        for j = n : -1 : i+1
            x(i) = x(i) - U(i,j) * x(j);  
        end
        x(i) = x(i) / U(i,i);
    end
    % disp(x)

end

%% Gram-Schmidt, Q R facterization
function  x  =  Gram_Schmidt(A,b)
    Q = zeros(size(A,1));
    R = zeros(size(A,2));
    for k = 1:size(A,2)
        Q(:,k) = A(:,k);
        for j = 1:k-1
            R(j,k) = Q(:,j)' * A(:,k);
            Q(:,k) = Q(:,k) - R(j,k)*Q(:,j);
        end
        R(k,k) = norm(Q(:,k),2);
        Q(:,k) = Q(:,k)/R(k,k);
    end
    %% such that Q^T*Q = diag(I,0),
    b = Q'*b;
    x = Backforward_sub(R,b); % upper triangle matrix R
end

function x = Backforward_sub(R,b)
    for i = size(R,2):-1:2
        x(i) = b(i) / R(i,i);       
        b(1:i-1) = b(1:i-1) - R(1:i-1,i)*x(i); % update b
    end
    x(1) = b(1) / R(1,1); 
    x = x';
end
