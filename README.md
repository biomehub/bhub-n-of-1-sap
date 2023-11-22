# n-of-1

To build the required docker image:

```
docker build -t postprandial-n-of-1:reproducibilityisawesome .
```

To run the simulation:

```
docker run -it -w ${PWD} -v "${PWD}:${PWD}" postprandial-n-of-1:reproducibilityisawesome R -e 'targets::tar_make()'
```

If everything is in order, `{targets}` should skip all steps - unless the code has been altered. To *actually* run the simulation, change the variable `overall_seed` within the `_targets.R` file to any integer. Then, run the line above again.