set ylabel '{/Symbol w}(meV)'
set xlabel 'Gamma - M'
set style fill transparent solid 0.25
pl for [i=0:9] 'MoS2_lw_path.d' every 9::i u 1:2 w l lt rgb 'black' title ''
repl for [i=0:9] 'MoS2_lw_path.d' every 9::i u ($1):($2-$3/2):($2+$3/2) w filledc lt rgb 'red' title ''
