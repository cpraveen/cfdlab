use LinearAlgebra;

proc main()
{
   var A = Matrix([2.0, 1.0],
                  [1.0, 2.0]);
   writeln("A =\n",A);
   var B = A+A;
   var C = dot(A,A);
   writeln("A*A =\n",C);
   var (evals, evecs) = eig(A, right=true);
   for i in 0..1 do
      writeln("eval = ", evals[i], " evec = ", evecs[..,i]);
}
