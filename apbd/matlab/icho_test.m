clear;
rng(12);

% Dense test
%A = randn(5,5);
%A = A*A';

% Sparse test
%A = sprandsym(10,0.5,1e-1);

% Load simulation data
load("Arch_states.mat");
A = A + 1e-6 * eye(size(A,1));
b =  b + 1e-4 * randn(size(A,1),1);

% Matlab's ilu() and ichol() only work with sparse, so convert to full for
% debugging purposes
options.udiag = true;
options.droptol = 1e-4;
[L,U] = ilu(sparse(A), options);
L = full(L);
U = full(U);

% For SPD A, dividing each row of U with the diag gives us L
fprintf('%e\n',norm((U./diag(U))' - L));

% Fpr SPD A, dividing each row of U with the sqrt of the diag gives us C
% C = full(ichol(sparse(A)));
% norm((U./sqrt(diag(U)))' - C)

[L0,U0] = ilu0(A);
fprintf('%e\n',norm(L - L0));
fprintf('%e\n',norm(U - U0));

[L0b,U0b] = ilu0_block(collisionList, size(A,1), 1e-6);
U0b = U0b + 1e-6*eye(size(U0b));
fprintf('%e\n',norm(L0b*U0b - A));

% These won't necessarily be zero for sparse matrices
fprintf('%e\n',norm(full(L*U-A)));
fprintf('%e\n',norm(full(L0*U0-A)));

freeIndices = ones(size(A,1),1);
selectM = zeros(size(A,1),length(find(freeIndices)));
selecti = 1;
for i = 1:3:size(A,1)
    if(freeIndices(i)&&freeIndices(i+1))
        selectM(i:i+2,selecti:selecti+2) = eye(3);
        selecti = selecti + 3;
    end
end

clf;
[x,~,relres,iters,resvec] = minres(selectM'*A*selectM,selectM'*b,1e-6,100);
semilogy(1:length(resvec), resvec, 'DisplayName', 'No Preconditioner','linewidth',2);
hold on;

[x,~,relres,iters,resvec] = minres(selectM'*A*selectM,selectM'*b,1e-6,100, selectM'*L0*selectM, selectM'*U0*selectM);
semilogy(1:length(resvec), resvec, 'DisplayName', 'Removed Incomplete LU Preconditioner','linewidth',2);
hold on;

[x,~,relres,iters,resvec] = minres(selectM'*A*selectM,selectM'*b,1e-6,100, selectM'*L0b*selectM, selectM'*U0b*selectM);
semilogy(1:length(resvec), resvec, 'DisplayName', 'Removed Block Incomplete LU Preconditioner','linewidth',2);
hold on;
legend
xlabel('Iteration number') 
ylabel('Residual') 
title('Convergence rate of different preconditioners for full matrix')


for i = [22 28 46 58 70 100 112 118 133 148 157 184 190]
    freeIndices(i:i+2) = 0;
end

selectM = zeros(size(A,1),length(find(freeIndices)));
selecti = 1;
for i = 1:3:size(A,1)
    if(freeIndices(i)&&freeIndices(i+1))
        selectM(i:i+2,selecti:selecti+2) = eye(3);
        selecti = selecti + 3;
    end
end

clf;
[x,~,relres,iters,resvec] = minres(selectM'*A*selectM,selectM'*b,1e-6,100);
semilogy(1:length(resvec), resvec, 'DisplayName', 'No Preconditioner','linewidth',2);
hold on;

[x,~,relres,iters,resvec] = minres(selectM'*A*selectM,selectM'*b,1e-6,100, selectM'*L0*selectM, selectM'*U0*selectM);
semilogy(1:length(resvec), resvec, 'DisplayName', 'Removed Incomplete LU Preconditioner','linewidth',2);
hold on;

[x,~,relres,iters,resvec] = minres(selectM'*A*selectM,selectM'*b,1e-6,100, selectM'*L0b*selectM, selectM'*U0b*selectM);
semilogy(1:length(resvec), resvec, 'DisplayName', 'Removed Block Incomplete LU Preconditioner','linewidth',2);
hold on;
legend
xlabel('Iteration number') 
ylabel('Residual') 
title('Convergence rate of different preconditioners for sub matrix')

%%
function [L,U] = ilu0(A,e)
if nargin < 2
    e = 0;
end
n = size(A,1);

L = zeros(n);
U = zeros(n);

% For each diagonal
for i = 1 : n
    % Compute row U(i,i:end)
    for j = i : n % NOTE: look for neighbors rather than looping through
        if A(i,j) ~= 0
            U(i,j) = A(i,j);
            % Subtract dot product of L(i,1:i-1) and U(1:i-1,j)
            for k = 1 : i-1 % NOTE: look for neighbors rather than looping through
                if U(k,j) ~= 0 % same as checking L(i,k)
                    U(i,j) = U(i,j) - L(i,k) * U(k,j);
                end
            end
        end
    end

    % Compute column L(i+1:end,i)
    L(i,i) = 1;
    for j = i + 1 : n % NOTE: look for neighbors rather than looping through
        if A(j,i) ~= 0
            L(j,i) = A(j,i);
            % Subtract dot product of L(j,1:i-1) and U(1:i-1,i)
            for k = 1 : i-1 % NOTE: look for neighbors rather than looping through
                if L(j,k) ~= 0 % same as checking U(k,i)
                    L(j,i) = L(j,i) - L(j,k) * U(k,i);
                end
            end
            % Divide by diagonal
            L(j,i) = L(j,i) / (U(i,i) + e);
        end
    end
end
end

%%
function [L,U] = ilu0_block(collision, n, e)
    nc = length(collision);
    %Compute neighbors for each contact
    for i = 1:nc
        collision{i}.neighbors = [];
        collision{i}.collisionTypes = [];
        for j = 1:nc
            if(collision{i}.body1 == collision{j}.body1)
                collision{i}.neighbors(end+1) = j;
                collision{i}.collisionTypes(end+1) = 1;
            elseif(collision{i}.body1 == collision{j}.body2)
                collision{i}.neighbors(end+1) = j;
                collision{i}.collisionTypes(end+1) = 2;
            elseif(collision{i}.body2 ~=0 && collision{i}.body2 == collision{j}.body1)
                collision{i}.neighbors(end+1) = j;
                collision{i}.collisionTypes(end+1) = 3;
            elseif(collision{i}.body2 ~=0 && collision{i}.body2 == collision{j}.body2)
                collision{i}.neighbors(end+1) = j;
                collision{i}.collisionTypes(end+1) = 4;
            end
        end
    end

    L = zeros(n);
    U = zeros(n);

    for i = 1 : nc
        % Compute row U(i,i:end)
        for j = collision{i}.neighbors
            if(j >= i)
                switch (collision{i}.collisionTypes(find(collision{i}.neighbors == j)))
                    case 1
                        U(collision{i}.mIndices,collision{j}.mIndices) = collision{i}.J1I * collision{j}.J1I';
                    case 2
                        U(collision{i}.mIndices,collision{j}.mIndices) = collision{i}.J1I * collision{j}.J2I';
                    case 3
                        U(collision{i}.mIndices,collision{j}.mIndices) = collision{i}.J2I * collision{j}.J1I';
                    case 4
                        U(collision{i}.mIndices,collision{j}.mIndices) = collision{i}.J2I * collision{j}.J2I';
                end
                if(i == j)
                    U(collision{i}.mIndices,collision{j}.mIndices) = collision{i}.J1I* collision{j}.J1I' + collision{i}.J2I* collision{j}.J2I';
                end
                
                for k = collision{i}.neighbors
                    if(k < i)
                        U(collision{i}.mIndices,collision{j}.mIndices) = U(collision{i}.mIndices,collision{j}.mIndices) ...
                            - L(collision{i}.mIndices,collision{k}.mIndices) * U(collision{k}.mIndices,collision{j}.mIndices);
                    end
                end
            end
        end

        % Compute column L(i+1:end,i)
        L(collision{i}.mIndices,collision{i}.mIndices) = eye(length(collision{i}.mIndices));
        for j = collision{i}.neighbors
            if(j > i)
                switch (collision{i}.collisionTypes(find(collision{i}.neighbors == j)))
                    case 1
                        L(collision{j}.mIndices,collision{i}.mIndices) = collision{j}.J1I * collision{i}.J1I';
                    case 2
                        L(collision{j}.mIndices,collision{i}.mIndices) = collision{j}.J2I * collision{i}.J1I';
                    case 3
                        L(collision{j}.mIndices,collision{i}.mIndices) = collision{j}.J1I * collision{i}.J2I';
                    case 4
                        L(collision{j}.mIndices,collision{i}.mIndices) = collision{j}.J2I * collision{i}.J2I';
                end
                
                for k = collision{i}.neighbors
                    if(k < i)
                        L(collision{j}.mIndices,collision{i}.mIndices) = L(collision{j}.mIndices,collision{i}.mIndices) ...
                            - L(collision{j}.mIndices,collision{k}.mIndices) * U(collision{k}.mIndices,collision{i}.mIndices);
                    end
                end
                % Divide by diagonal
                L(collision{j}.mIndices,collision{i}.mIndices) = ((U(collision{i}.mIndices,collision{i}.mIndices) + e * eye(length(collision{i}.mIndices)))\L(collision{j}.mIndices,collision{i}.mIndices)')';
            end
        end
    end
end

%% Simple version for debugging
% function [L,U] = ilu0(A,e)
% if nargin < 2
%     e = 0;
% end
% n = size(A,1);
% for k = 1 : n
%     for i = k + 1 : n
%         A(i,k) = A(i,k)/(A(k,k) + e);
%         for j = k + 1 : n
%             A(i,j) = A(i,j) - A(i,k)*A(k,j);
%         end
%     end
% end
% L = tril(A);
% U = triu(A);
% end
