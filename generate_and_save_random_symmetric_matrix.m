function generate_and_save_random_symmetric_matrix(n, folder_path)
    
    % Generate a random symmetric matrix
   % A = randi([-1, 1], n, n);  % Create a random matrix with integer values between 1 and 10
    A = rand(n, n);  
   symmetric_matrix_1 = (A + A') / 2;  % Make it symmetric by averaging with its transpose
   %const=10; 
   
 
% finding eigenvalue of the matrix H
[V,E]=eig(A);
%E_1=abs(E); % finding the absolute value of eigenvalues of H
lambda_max_1=max(sort(diag(A)))% find the exact max eigenvalue of B
lambda_min_1=min(sort(diag(A)))
cond_1=cond(A);
norm_1=norm(A)

% % chrusing the matrix with positive const g (B=A./g), so  norm(B,p)<1 is
% true
b=vecnorm(A,1,2); % first norm on rows
%g=const*max(b); % a const g that smash matrix A
g=max(b);
symmetric_matrix_2=symmetric_matrix_1./g;

% finding eigenvalue of the matrix H
[V1,E1]=eig(symmetric_matrix_2);
%E_1=abs(E); % finding the absolute value of eigenvalues of H
lambda_max_2=max(sort(diag(symmetric_matrix_2)))% find the exact max eigenvalue of B
lambda_min_2=min(sort(diag(symmetric_matrix_2)))
cond_2=cond(symmetric_matrix_2);
norm_2=norm(symmetric_matrix_2)

% Generate Circulant Symmetric Matrix
    first_row = randi([1, 10], 1, n);  % Random values for the first row
  % first_row = rand(n, n);
    circulant_matrix = toeplitz([first_row(1), first_row(end:-1:2)], first_row);
    
% finding eigenvalue of the matrix Circulant matrix
[V2,E2]=eig(circulant_matrix);
%E_1=abs(E); % finding the absolute value of eigenvalues of H
lambda_max_Circulant=max(sort(diag(circulant_matrix)))% find the exact max eigenvalue of B
lambda_min_Circulant=min(sort(diag(circulant_matrix)))
cond_Circulant=cond(circulant_matrix);
norm_Circulant=norm(circulant_matrix)
    
    % Generate Toeplitz Symmetric Matrix
    first_row_toeplitz = randi([1, 10],1, n);  % Random values for example
    toeplitz_matrix = toeplitz(first_row_toeplitz);
    
    % finding eigenvalue of the matrix Toeplitz matrix
[V3,E3]=eig(toeplitz_matrix);
%E_1=abs(E); % finding the absolute value of eigenvalues of H
lambda_max_Toeplitz=max(sort(diag(toeplitz_matrix)))% find the exact max eigenvalue of B
lambda_min_Toeplitz=min(sort(diag(toeplitz_matrix)))
cond_Toeplitz=cond(toeplitz_matrix);
norm_Toeplitz=norm(toeplitz_matrix)
    
% %     % Generate Tridiagonal Symmetric Matrix
%    main_diag = randi([1, 10],n,n);  % Main diagonal elements
%   off_diag = randi([1, 10], n-1, 1); % Off diagonal elements
%  tridiagonal_matrix = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);
% %     
% %     
% %         % finding eigenvalue of the matrix Tridiagonal matrix
%  [V4,E4]=eig( tridiagonal_matrix);
% E_1=abs(E); % finding the absolute value of eigenvalues of H
%  lambda_max_Tr=max(sort(diag(tridiagonal_matrix)))% find the exact max eigenvalue of B
%  lambda_min_Tr=min(sort(diag(tridiagonal_matrix)))
% cond_Tr=cond(tridiagonal_matrix);
% norm_Tr=norm(tridiagonal_matrix)
% %     
       
%     % Generate Centrosymmetric Matrix
%     half_matrix = randi([-1, 1], ceil(n/2), n);  % Random values for the first half
%     centrosymmetric_matrix = [half_matrix; flipud(half_matrix)];
%     if mod(n, 2) == 1  % If n is odd, remove duplicated middle row
%         centrosymmetric_matrix = centrosymmetric_matrix(1:end-1, :);
%     end
%     
%         % finding eigenvalue of the matrix Centrosymmetric matrix
% [V5,E5]=eig(centrosymmetric_matrix);
% %E_1=abs(E); % finding the absolute value of eigenvalues of H
% lambda_max_Ce=max(sort(diag(centrosymmetric_matrix)))% find the exact max eigenvalue of B
% lambda_min_Ce=min(sort(diag(centrosymmetric_matrix)))
% cond_Ce=cond(centrosymmetric_matrix);
% norm_Ce=norm(centrosymmetric_matrix)
%     
    
%     % Generate Coxeter Matrix (for a simple Coxeter group, e.g., A_n)
%     coxeter_matrix = 2 * eye(n) - diag(ones(n-1, 1), 1) - diag(ones(n-1, 1), -1);
%     
%             % finding eigenvalue of the matrix Coxeter matrix
% [V6,E6]=eig(coxeter_matrix);
% %E_1=abs(E); % finding the absolute value of eigenvalues of H
% lambda_max_Cox=max(sort(diag(coxeter_matrix)))% find the exact max eigenvalue of B
% lambda_min_Cox=min(sort(diag(coxeter_matrix)))
% cond_Cox=cond(coxeter_matrix);
% norm_Cox=norm(coxeter_matrix)
%     
%     
%     % Generate GCD Matrix
%     gcd_matrix = zeros(n, n);
%     for i = 1:n
%         for j = 1:n
%             gcd_matrix(i, j) = gcd(i, j);  % Calculate GCD of indices
%         end
%     end
%     
%     
%             % finding eigenvalue of the matrix GCD matrix
% [V7,E7]=eig( gcd_matrix );
% %E_1=abs(E); % finding the absolute value of eigenvalues of H
% lambda_max_GCD=max(sort(diag( gcd_matrix )))% find the exact max eigenvalue of B
% lambda_min_GCD=min(sort(diag( gcd_matrix )))
% cond_GCD=cond( gcd_matrix );
% norm_GCD=norm( gcd_matrix )
%     
    
    % Generate Hankel Matrix
    first_row_hankel = randi([1, 10],1, n);
    last_column_hankel = [first_row_hankel(2:end), randi([1, 10])];
    hankel_matrix = hankel(first_row_hankel, last_column_hankel);
    
           % finding eigenvalue of the Hankel matrix
[V8,E8]=eig( hankel_matrix);
%E_1=abs(E); % finding the absolute value of eigenvalues of H
lambda_max_Hankel=max(sort(diag( hankel_matrix)))% find the exact max eigenvalue of B
lambda_min_Hankel=min(sort(diag( hankel_matrix)))
cond_Hankel=cond( hankel_matrix);
norm_Hanekl=norm( hankel_matrix)
     
    
% %     % Generate Hilbert Matrix
% %     hilbert_matrix = hilb(n);
% 
%     % finding eigenvalue of the Hilbert matrix
% [V9,E9]=eig( hilbert_matrix);
% %E_1=abs(E); % finding the absolute value of eigenvalues of H
% lambda_max_Hi=max(sort(diag( hilbert_matrix)))% find the exact max eigenvalue of B
% lambda_min_Hi=min(sort(diag( hilbert_matrix)))
% cond_Hi=cond( hilbert_matrix);
% norm_Hi=norm( hilbert_matrix)
    
    
%     
%     % Generate Persymmetric Matrix
%     base_values = rand( n, n);  % Random values for base
%     persymmetric_matrix = base_values + flipud(fliplr(base_values));
% 
%         % finding eigenvalue of the matrix Centrosymmetric matrix
% [V10,E10]=eig( persymmetric_matrix);
% %E_1=abs(E); % finding the absolute value of eigenvalues of H
% lambda_max_Per=max(sort(diag( persymmetric_matrix)))% find the exact max eigenvalue of B
% lambda_min_Per=min(sort(diag( persymmetric_matrix)))
% cond_Per=cond( persymmetric_matrix);
% norm_Per=norm( persymmetric_matrix);
%     


    % Construct the full file path for the Excel file
    filename_1 = fullfile('D:\Documents_for_Ph_D\Disertaciq-versions\chapter1', 'random_symmetric_matrix_A.xlsx');
    filename_2 = fullfile('D:\Documents_for_Ph_D\Disertaciq-versions\chapter1', 'random_symmetric_matrix_A_smash.xlsx');
    
    filename_3 = fullfile('D:\Documents_for_Ph_D\Disertaciq-versions\chapter1', 'random_symmetric_matrix_circulant.xlsx');
    filename_4 = fullfile('D:\Documents_for_Ph_D\Disertaciq-versions\chapter1', 'random_symmetric_matrix_toeplitz.xlsx');
    
    filename_5 = fullfile('D:\Documents_for_Ph_D\Disertaciq-versions\chapter1', 'random_symmetric_matrix_tridiagonal.xlsx');
    %filename_6 = fullfile('D:\Documents_for_Ph_D\Disertaciq-versions\chapter1', 'random_symmetric_matrix_centrosymmetric.xlsx');
    
    %filename_7 = fullfile('D:\Documents_for_Ph_D\Disertaciq-versions\chapter1', 'random_symmetric_matrix_coxeter.xlsx');
    %filename_8 = fullfile('D:\Documents_for_Ph_D\Disertaciq-versions\chapter1', 'random_symmetric_matrix_gcd.xlsx');
    
    filename_9 = fullfile('D:\Documents_for_Ph_D\Disertaciq-versions\chapter1', 'random_symmetric_matrix_hankel.xlsx');
   % filename_10= fullfile('D:\Documents_for_Ph_D\Disertaciq-versions\chapter1', 'random_symmetric_matrix_hilbert.xlsx');
    %filename_11 = fullfile('D:\Documents_for_Ph_D\Disertaciq-versions\chapter1', 'random_symmetric_matrix_persymmetric.xlsx');
  
    % Save the symmetric matrix to the specified folder in Excel format
    writematrix(symmetric_matrix_1, filename_1, 'Sheet', 'SymmetricMatrix');
    writematrix(symmetric_matrix_2, filename_2, 'Sheet', 'SymmetricMatrix');
    writematrix(circulant_matrix, filename_3, 'Sheet', 'Circulant');
    writematrix(toeplitz_matrix, filename_4, 'Sheet', 'Toeplitz');
  %  writematrix(tridiagonal_matrix, filename_5, 'Sheet', 'Tridiagonal');
   % writematrix(centrosymmetric_matrix, filename_6, 'Sheet', 'Centrosymmetric');
   % writematrix(coxeter_matrix, filename_7, 'Sheet', 'Coxeter');
   % writematrix(gcd_matrix, filename_8, 'Sheet', 'GCD');
    writematrix(hankel_matrix, filename_9, 'Sheet', 'Hankel');
    %writematrix(hilbert_matrix, filename_10, 'Sheet', 'Hilbert');
    %writematrix(persymmetric_matrix, filename_11, 'Sheet', 'Persymmetric');
    
    fprintf('Random symmetric matrix has been saved to %s\n', filename);
end
