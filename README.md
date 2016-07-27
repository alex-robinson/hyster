# Hyster: a hysteresis package

This code facilitates running transient simulations
of hysteresis by regulating the rate of change of the 
forcing as a function of the rate of change of the 
result.

# Testing
```
make clean
make test env=[airaki,manto,...]
./test_hyster.x > tmp
```

In R, plot the results here:
```
source("plot_test.r")
```

