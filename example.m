%A simple example
A = load('toy_matrix.txt')
k = 3

%Run the algorithm and compute quadratic loss
R = iter_fs(A, k)
C = A(:,R)
loss = norm(A-C*pinv(C)*A, 'fro')^2

%Compare to random choice
R = randsample(size(A,2), k)
C = A(:,R)
loss = norm(A-C*pinv(C)*A, 'fro')^2
