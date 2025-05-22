function functionf = alpha_optimal(x)
load lamda;
functionf = abs((lamda + x)/(lamda - x));
end