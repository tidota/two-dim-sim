# Two Dimensional Simulation

This is simulation of cooperative localization in a two dimentional space.
Each robot estimates its own location. They can communicate with their neighbors and exchange their estimates.
The estimates are updated by the exchanged data and the detected distances to them.

- robots are placed in a 2-D space
- each one moves based on the distances to its neighbors (maybe following the next?)
- motion model has noises
- pass its estimates and update its own with those from its neighbors
- only one can perform global localization

# TODO
- [ ] Update README for 2D sim
- [ ] Update yaml for 2D sim
- [ ] Import the update code from 1D sim
- [ ] Update the code for 2D sim

# Pseudocode

init locations
each step
  each robot
    decide control
    execute the control
    predict the location
  each robot
    take measurement (communication)
    update its location estimate

# gnuplot notes

```
gnuplot
plot "output.dat" u 1:2 title "x1", "output.dat" u 1:5 title "x2", "output.dat" u 1:8 title "x3", "output.dat" u 1:3:4 with errorbars title "m1", "output.dat" u 1:6:7 with errorbars title "m2", "output.dat" u 1:9:10 with errorbars title "m3"
```
```
plot "output.dat" u 1:2 title "x1", "output.dat" u 1:5 title "x2", "output.dat" u 1:8 title "x3", "output.dat" u 1:11 title "x4", "output.dat" u 1:14 title "x5", "output.dat" u 1:3:4 with errorbars title "m1", "output.dat" u 1:6:7 with errorbars title "m2", "output.dat" u 1:9:10 with errorbars title "m3", "output.dat" u 1:12:13 with errorbars title "m4", "output.dat" u 1:15:16 with errorbars title "m5"
```

