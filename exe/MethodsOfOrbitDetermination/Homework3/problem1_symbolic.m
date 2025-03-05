clc; clear all; close all;

x = sym("x", ["real"]);
y = sym("y", ["real"]);

vx = sym("vx", ["real"]);
vy = sym("vy", ["real"]);


r = sqrt(x^2 + y^2);
fx = -x / r^3;
fy = -y / r^3;

pretty(simplify(jacobian(fx, [x, y])))
pretty(simplify(jacobian(fy, [x, y])))