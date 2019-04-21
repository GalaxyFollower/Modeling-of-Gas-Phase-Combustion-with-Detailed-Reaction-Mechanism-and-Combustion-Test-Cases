function x = gauss_elimination ( A, b )

%GAUSS_ELIM   solve the linear system Ax = b using Gaussian elimination
%             with back substitution
%
%     calling sequences:
%             x = gauss_elim ( A, b )
%             gauss_elim ( A, b )
%
%     inputs:
%             A       coefficient matrix for linear system
%                     (matrix must be square)
%             b       right-hand side vector
%
%     output:
%             x       solution vector (i.e., vector for which Ax = b)
%
%     NOTE:
%             this is intended as a demonstration routine - no pivoting
%             strategy is implemented to reduce the effects of roundoff
%             error
%

[nrow ncol] = size ( A );
if ( nrow ~= ncol )
   disp ( 'gauss_elim error: Square coefficient matrix required' );
   return;
end;
nb = length ( b );
if ( nrow ~= nb )
   disp ( 'gauss_elim error: Size of b-vector not compatible with matrix dimension' )
   return;
end;

x = zeros ( 1, nrow );

%
%    Gaussian elimination
%

for i = 1 : nrow - 1
    if ( A(i,i) == 0 )
	   t =  min ( find ( A(i+1:nrow,i) ~= 0 ) + i );
	   if ( isempty(t) ) 
	      disp ( 'gauss_elim error: A matrix is singular' );
		  return
	   end;
	   temp = A(i,:);     tb = b(i);
	   A(i,:) = A(t,:);   b(i) = b(t);
	   A(t,:) = temp;     b(t) = tb;
	end;
    for j = i+1 : nrow
	    m = -A(j,i) / A(i,i);
		A(j,i) = 0;
		A(j, i+1:nrow) = A(j, i+1:nrow) + m * A(i, i+1:nrow);
		b(j) = b(j) + m * b(i);
    end;
end;

%
%    back substitution
%

x(nrow) = b(nrow) / A(nrow, nrow);
for i = nrow - 1 : -1 : 1
    x(i) = ( b(i) - sum ( x(i+1:nrow) .* A(i, i+1:nrow) ) ) / A(i,i);
    
end;