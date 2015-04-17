function M = stima3 ( vertices )

  d = size ( vertices, 2 );

  G = [ ones(1,d+1); vertices' ] \ [ zeros(1,d); eye(d) ];

  M = det ( [ ones(1,d+1); vertices' ] ) * G * G' / prod ( 1:d );

end
