library(gpuR)
library(sRDA)

# generate data with few highly correlated variahbles
dataXY <- generate_data(nr_LVs = 2,
                           n = 550,
                           nr_correlated_Xs = c(5,20),
                           nr_uncorrelated_Xs = 450,
                           mean_reg_weights_assoc_X =
                             c(0.9,0.5),
                           sd_reg_weights_assoc_X =
                             c(0.05, 0.05),
                           Xnoise_min = -0.3,
                           Xnoise_max = 0.3,
                           nr_correlated_Ys = c(10,15),
                           nr_uncorrelated_Ys = 500,
                           mean_reg_weights_assoc_Y =
                             c(0.9,0.6),
                           sd_reg_weights_assoc_Y =
                             c(0.05, 0.05),
                           Ynoise_min = -0.3,
                           Ynoise_max = 0.3)

# seperate predictor and predicted sets
X <- dataXY$X
Y <- dataXY$Y

dim(X)
dim(Y)

predictor = X
predicted = Y
tolerance = 1*10^-20
max_iterations = 100
