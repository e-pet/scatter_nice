## scatter_nice

A wrapper for scatter(.) that
1) chooses useful defaults for colors and marker size and shape,
2) allows for easy group-wise scatter plots,
3) implements randomized-order plotting to prevent the last plotted data to (misleadingly) dominate plot appearance,
4) uses nonlinear color scaling to exploit the full color range, and
5) implements some automatic simple statistical annotations to be shown (optionally).

The goal is for this to be a one-stop function to easily create a useful and visually appealing scatter diagram.

Requires the `cbrewer` file exchange function to be on the path, see 
https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab.
The file is included here for convenience; copyright remains with the original author.

The following are two examples comparing the default outputs of `scatter(.)` and `scatter_nice(.)`.

![An example plot.](demo.png)

---
Eike Petersen, 2020-2021
