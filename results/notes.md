# Notes for gnuplot to generate the Sigma for each method.

```
plot "output_cent.dat" u 1:54 title "Centralized" with line
replot "output_naive.dat" u 1:54 title "Naive" with line
replot "output_ci.dat" u 1:54 title "CI-based" with line
replot "output_proposed.dat" u 1:54 title "Proposed" with line
set yrange [-0.5:4]
```

# Notes for gnuplot to generate the error fig for each method.

```
clear
unset object
plot "output_naive.dat" u 1:11 title "Robot1" with line, \
     "output_naive.dat" u 1:22 title "Robot2" with line, \
     "output_naive.dat" u 1:33 title "Robot3" with line, \
     "output_naive.dat" u 1:44 title "Robot4" with line, \
     "output_naive.dat" u 1:55 title "Robot5" with line, \
     "output_naive.dat" u 1:66 title "Robot6" with line, \
     "output_naive.dat" u 1:77 title "Robot7" with line, \
     "output_naive.dat" u 1:88 title "Robot8" with line

clear
unset object
plot "output_ci.dat" u 1:11 title "Robot1" with line, \
     "output_ci.dat" u 1:22 title "Robot2" with line, \
     "output_ci.dat" u 1:33 title "Robot3" with line, \
     "output_ci.dat" u 1:44 title "Robot4" with line, \
     "output_ci.dat" u 1:55 title "Robot5" with line, \
     "output_ci.dat" u 1:66 title "Robot6" with line, \
     "output_ci.dat" u 1:77 title "Robot7" with line, \
     "output_ci.dat" u 1:88 title "Robot8" with line

clear
unset object
plot "output_proposed.dat" u 1:11 title "Robot1" with line, \
     "output_proposed.dat" u 1:22 title "Robot2" with line, \
     "output_proposed.dat" u 1:33 title "Robot3" with line, \
     "output_proposed.dat" u 1:44 title "Robot4" with line, \
     "output_proposed.dat" u 1:55 title "Robot5" with line, \
     "output_proposed.dat" u 1:66 title "Robot6" with line, \
     "output_proposed.dat" u 1:77 title "Robot7" with line, \
     "output_proposed.dat" u 1:88 title "Robot8" with line

clear
unset object
plot "output_cent.dat" u 1:11 title "Robot1" with line, \
     "output_cent.dat" u 1:22 title "Robot2" with line, \
     "output_cent.dat" u 1:33 title "Robot3" with line, \
     "output_cent.dat" u 1:44 title "Robot4" with line, \
     "output_cent.dat" u 1:55 title "Robot5" with line, \
     "output_cent.dat" u 1:66 title "Robot6" with line, \
     "output_cent.dat" u 1:77 title "Robot7" with line, \
     "output_cent.dat" u 1:88 title "Robot8" with line

clear
unset object
plot "output_pf_ori_naive.dat" u 1:8 title "Robot1" with line, \
     "output_pf_ori_naive.dat" u 1:18 title "Robot2" with line, \
     "output_pf_ori_naive.dat" u 1:28 title "Robot3" with line, \
     "output_pf_ori_naive.dat" u 1:38 title "Robot4" with line, \
     "output_pf_ori_naive.dat" u 1:48 title "Robot5" with line, \
     "output_pf_ori_naive.dat" u 1:58 title "Robot6" with line, \
     "output_pf_ori_naive.dat" u 1:68 title "Robot7" with line, \
     "output_pf_ori_naive.dat" u 1:78 title "Robot8" with line

clear
unset object
plot "output_pf_ori_proposed.dat" u 1:8 title "Robot1" with line, \
     "output_pf_ori_proposed.dat" u 1:18 title "Robot2" with line, \
     "output_pf_ori_proposed.dat" u 1:28 title "Robot3" with line, \
     "output_pf_ori_proposed.dat" u 1:38 title "Robot4" with line, \
     "output_pf_ori_proposed.dat" u 1:48 title "Robot5" with line, \
     "output_pf_ori_proposed.dat" u 1:58 title "Robot6" with line, \
     "output_pf_ori_proposed.dat" u 1:68 title "Robot7" with line, \
     "output_pf_ori_proposed.dat" u 1:78 title "Robot8" with line
```

```
clear
unset object
plot "output_pf_ori_naive.dat" u 1:10 title "Robot1" with line, \
     "output_pf_ori_naive.dat" u 1:20 title "Robot2" with line, \
     "output_pf_ori_naive.dat" u 1:30 title "Robot3" with line, \
     "output_pf_ori_naive.dat" u 1:40 title "Robot4" with line, \
     "output_pf_ori_naive.dat" u 1:50 title "Robot5" with line, \
     "output_pf_ori_naive.dat" u 1:60 title "Robot6" with line, \
     "output_pf_ori_naive.dat" u 1:70 title "Robot7" with line, \
     "output_pf_ori_naive.dat" u 1:80 title "Robot8" with line

clear
unset object
plot "output_pf_ori_proposed.dat" u 1:10 title "Robot1" with line, \
     "output_pf_ori_proposed.dat" u 1:20 title "Robot2" with line, \
     "output_pf_ori_proposed.dat" u 1:30 title "Robot3" with line, \
     "output_pf_ori_proposed.dat" u 1:40 title "Robot4" with line, \
     "output_pf_ori_proposed.dat" u 1:50 title "Robot5" with line, \
     "output_pf_ori_proposed.dat" u 1:60 title "Robot6" with line, \
     "output_pf_ori_proposed.dat" u 1:70 title "Robot7" with line, \
     "output_pf_ori_proposed.dat" u 1:80 title "Robot8" with line
```

# To move the legends

```
set key left top
```

# individual plots for naive

```
clear
unset object
plot "output_pf_ori_naive.dat" u 1:8 title "Err" with line
replot "output_pf_ori_naive.dat" u 1:9 title "Uncertainty" with line

clear
unset object
plot "output_pf_ori_naive.dat" u 1:18 title "Err" with line
replot "output_pf_ori_naive.dat" u 1:19 title "Uncertainty" with line

clear
unset object
plot "output_pf_ori_naive.dat" u 1:28 title "Err" with line
replot "output_pf_ori_naive.dat" u 1:29 title "Uncertainty" with line

clear
unset object
plot "output_pf_ori_naive.dat" u 1:38 title "Err" with line
replot "output_pf_ori_naive.dat" u 1:39 title "Uncertainty" with line

clear
unset object
plot "output_pf_ori_naive.dat" u 1:48 title "Err" with line
replot "output_pf_ori_naive.dat" u 1:49 title "Uncertainty" with line

clear
unset object
plot "output_pf_ori_naive.dat" u 1:58 title "Err" with line
replot "output_pf_ori_naive.dat" u 1:59 title "Uncertainty" with line

clear
unset object
plot "output_pf_ori_naive.dat" u 1:68 title "Err" with line
replot "output_pf_ori_naive.dat" u 1:69 title "Uncertainty" with line

clear
unset object
plot "output_pf_ori_naive.dat" u 1:78 title "Err" with line
replot "output_pf_ori_naive.dat" u 1:79 title "Uncertainty" with line
```

```
clear
unset object
set yrange [-180:180]
plot "output_pf_ori_naive.dat" u 1:10 with line

clear
unset object
set yrange [-180:180]
plot "output_pf_ori_naive.dat" u 1:20 with line

clear
unset object
set yrange [-180:180]
plot "output_pf_ori_naive.dat" u 1:30 with line

clear
unset object
set yrange [-180:180]
plot "output_pf_ori_naive.dat" u 1:40 with line

clear
unset object
set yrange [-180:180]
plot "output_pf_ori_naive.dat" u 1:50 with line

clear
unset object
set yrange [-180:180]
plot "output_pf_ori_naive.dat" u 1:60 with line

clear
unset object
set yrange [-180:180]
plot "output_pf_ori_naive.dat" u 1:70 with line

clear
unset object
set yrange [-180:180]
plot "output_pf_ori_naive.dat" u 1:80 with line
```

# individual plots for proposed

```
clear
unset object
plot "output_pf_ori_proposed.dat" u 1:8 title "Err" with line
replot "output_pf_ori_proposed.dat" u 1:9 title "Uncertainty" with line

clear
unset object
plot "output_pf_ori_proposed.dat" u 1:18 title "Err" with line
replot "output_pf_ori_proposed.dat" u 1:19 title "Uncertainty" with line

clear
unset object
plot "output_pf_ori_proposed.dat" u 1:28 title "Err" with line
replot "output_pf_ori_proposed.dat" u 1:29 title "Uncertainty" with line

clear
unset object
plot "output_pf_ori_proposed.dat" u 1:38 title "Err" with line
replot "output_pf_ori_proposed.dat" u 1:39 title "Uncertainty" with line

clear
unset object
plot "output_pf_ori_proposed.dat" u 1:48 title "Err" with line
replot "output_pf_ori_proposed.dat" u 1:49 title "Uncertainty" with line

clear
unset object
plot "output_pf_ori_proposed.dat" u 1:58 title "Err" with line
replot "output_pf_ori_proposed.dat" u 1:59 title "Uncertainty" with line

clear
unset object
plot "output_pf_ori_proposed.dat" u 1:68 title "Err" with line
replot "output_pf_ori_proposed.dat" u 1:69 title "Uncertainty" with line

clear
unset object
plot "output_pf_ori_proposed.dat" u 1:78 title "Err" with line
replot "output_pf_ori_proposed.dat" u 1:79 title "Uncertainty" with line
```

```
clear
unset object
set yrange [-180:180]
plot "output_pf_ori_proposed.dat" u 1:10 with line

clear
unset object
set yrange [-180:180]
plot "output_pf_ori_proposed.dat" u 1:20 with line

clear
unset object
set yrange [-180:180]
plot "output_pf_ori_proposed.dat" u 1:30 with line

clear
unset object
set yrange [-180:180]
plot "output_pf_ori_proposed.dat" u 1:40 with line

clear
unset object
set yrange [-180:180]
plot "output_pf_ori_proposed.dat" u 1:50 with line

clear
unset object
set yrange [-180:180]
plot "output_pf_ori_proposed.dat" u 1:60 with line

clear
unset object
set yrange [-180:180]
plot "output_pf_ori_proposed.dat" u 1:70 with line

clear
unset object
set yrange [-180:180]
plot "output_pf_ori_proposed.dat" u 1:80 with line
```

