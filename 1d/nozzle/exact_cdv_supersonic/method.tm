<TeXmacs|1.0.7>

<style|article>

<\body>
  <doc-data|<doc-title|Quasi 1-D nozzle solution>|<doc-author-data|<author-name|Praveen.
  c>>>

  Consider a quasi 1-D nozzle with area <math|A(x)>. Let the subscripts
  <math|i,l,r,o> denote conditions at inflow, to the left of the shock, to
  the right of the shock and the outflow. All conditions at inflow are given
  and outflow pressure <math|p<rsub|o>> is also given. The shock is at some
  unknown location <math|x<rsub|s>> which must be determined.

  <\enumerate-numeric>
    <item>State <math|i> and <math|l> are connected adiabatically

    <\eqnarray*>
      <tformat|<twith|table-halign|l>|<table|<row|<cell|<frac|p<rsub|i>|\<rho\><rsub|i><rsup|\<gamma\>>>-<frac|p<rsub|l>|\<rho\><rsub|l><rsup|\<gamma\>>>=0>|<cell|>|<cell|>>|<row|<cell|<left|[><frac|\<gamma\>p<rsub|i>|(\<gamma\>-1)\<rho\><rsub|i>>+<frac|1|2>u<rsub|i><rsup|2><right|]>-<mid|[><frac|\<gamma\>p<rsub|l>|(\<gamma\>-1)\<rho\><rsub|l>>+<frac|1|2>u<rsub|l><rsup|2><right|]>=0>|<cell|>|<cell|>>|<row|<cell|\<rho\><rsub|i>u<rsub|i>A<rsub|i>-\<rho\><rsub|l>u<rsub|l>A(x<rsub|s>)=0>|<cell|>|<cell|>>>>
    </eqnarray*>

    <item>State <math|l> and <math|r> are connected by the Rankine-Hugoniot
    relations

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<rho\><rsub|l>u<rsub|l>-\<rho\><rsub|r>u<rsub|r>=0>|<cell|>|<cell|>>|<row|<cell|(p<rsub|l>+\<rho\><rsub|l>u<rsub|l><rsup|2>)-(p<rsub|r>+\<rho\><rsub|r>u<rsub|r><rsup|2>)=0>|<cell|>|<cell|>>|<row|<cell|<left|[><frac|\<gamma\>p<rsub|l>|(\<gamma\>-1)\<rho\><rsub|l>>+<frac|1|2>u<rsub|l><rsup|2><right|]>-<mid|[><frac|\<gamma\>p<rsub|r>|(\<gamma\>-1)\<rho\><rsub|r>>+<frac|1|2>u<rsub|r><rsup|2><right|]>=0>|<cell|>|<cell|>>>>
    </eqnarray*>

    <item>State <math|r> and <math|o> are connected adiabatically

    <\eqnarray*>
      <tformat|<table|<row|<cell|<frac|p<rsub|r>|\<rho\><rsub|r><rsup|\<gamma\>>>-<frac|p<rsub|o>|\<rho\><rsub|o><rsup|\<gamma\>>>=0>|<cell|>|<cell|>>|<row|<cell|<left|[><frac|\<gamma\>p<rsub|r>|(\<gamma\>-1)\<rho\><rsub|r>>+<frac|1|2>u<rsub|r><rsup|2><right|]>-<mid|[><frac|\<gamma\>p<rsub|o>|(\<gamma\>-1)\<rho\><rsub|o>>+<frac|1|2>u<rsub|o><rsup|2><right|]>=0>|<cell|>|<cell|>>|<row|<cell|\<rho\><rsub|r>u<rsub|r>A(x<rsub|s>)-\<rho\><rsub|o>u<rsub|o>A<rsub|o>=0>|<cell|>|<cell|>>>>
    </eqnarray*>
  </enumerate-numeric>

  We have nine equations for the nine unknowns
  <math|Q=(\<rho\><rsub|l>,u<rsub|l>,p<rsub|l>,\<rho\><rsub|r>,u<rsub|r>,p<rsub|r>,\<rho\><rsub|o>,u<rsub|o>,x<rsub|s>)>.
  These are nine coupled, non-linear equations which can be solved using a
  Newton algorithm. Symbolically we can write

  <\equation*>
    N(Q,M<rsub|i>)=0
  </equation*>

  where <math|M<rsub|i>> is the inflow Mach number.\ 
</body>