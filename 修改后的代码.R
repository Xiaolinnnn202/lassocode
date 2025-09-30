library(msgps)
library(Matrix)
library(foreach)
library(glmnet)
library(fUnitRoots)
library(splines)
library(rlang)
library(ggplot2)
library(zoo)
library(forecast)
library(tseries)
library(lars)
library(corpcor)
library(longitudinal)
library(fdrtool)
library(GeneNet)
#因为parcor包已经被删掉了，所以只能用网上的2014年安装包
install.packages("https://cran.r-project.org/src/contrib/Archive/ppls/ppls_1.6-1.1.tar.gz",repos = NULL, type = "source")
library(ppls)
library(MASS)
library(Epi)
install.packages("https://cran.r-project.org/src/contrib/Archive/parcor/parcor_0.2-6.tar.gz",repos = NULL, type = "source")
library(parcor)
library(TTR)
library(xts)
library(quantmod)

# 获取股票数据
stock_data <- read.csv("C:/__.csv")#这里填写你所引用的数据的文件路径
stock_data
# 响应变量
y <- as.matrix(stock_data[,-1]);y

## 创建时序图
days <- stock_data[,1];days
# 创建时间序列对象
time_series <- ts(y)
# 绘制时序图
plot(days, time_series, type = "l", xlab = "2019年1月1日至2020年12月31日", ylab = "中证从2019年1月1日至2020年12月31日收盘价")

## 用ADF单位根检验法检验时序是否通过平稳性检验
# 导入ADF单位根检验需要的adf.test()函数
library(tseries)
# 进行ADF检验
result1 <- adf.test(time_series)
# 打印检验结果
print(result1)

# 一阶差分处理
diff_time_series <- diff(time_series)
# 绘制差分后的时序图
plot(diff_time_series, type = "l", xlab = "2019年1月1日至2020年12月31日", ylab = "中证100从2019年1月1日至2020年12月31日收盘价一阶对数差分值")
# 对差分后的时序进行ADF检验
result2 <- adf.test(diff_time_series)
# 打印检验结果
print(result2)

diff_diff_time_series <- diff(time_series)


# 绘制样本自相关图，先将lag.max定为22
acf(diff_time_series,lag.max=22,main="series diff",xlab="Lag",ylab="ACF")
# 绘制样本偏自相关图
pacf(diff_time_series,lag.max=22,main="series diff",xlab="Lag",ylab="PACF")

#模型的自相关图为截尾，因此模型应该为ARMA模型，但是偏自相关图呈现震动状态，因此需要借助MBIC来给模型定阶。

##计算MBIC值
ts_data <- diff_time_series;ts_data
calculate_mbic <- function(ts_data, order_range, c_value, n) {
  mbic_values <- numeric(length(order_range))
  ln_n <- log(n)
  c_ln_ln_n <- c_value * log(ln_n)
  for (i in seq_along(order_range)) {
    fit <- arima(ts_data, order = c(order_range[i], 0, 0), method = "ML")
    log_likelihood <- logLik(fit)  # 模型的最大似然函数的对数值
    k <- length(fit$coef)  # 模型的未知参数个数
    mbic_values[i] <- -2 * log_likelihood + c_ln_ln_n*k
  }
  return(mbic_values)
}

# 设置阶数范围
order_range <- 1:22  # 假设阶数范围为1到22

# 给定的参数值
c_value <- 1.01
n <- 487

# 计算MBIC值
mbic_values <- calculate_mbic(ts_data, order_range, c_value, n)
# 画出MBIC值图
plot(order_range, mbic_values, type = "b", xlab = "滞后阶数", ylab = "MBIC值", main = "MBIC Values for AR Models")



##准备好ARMA(8,3)数据
# 提取MA部分的数据
# 1. 拟合一个适当的模型到原始时间序列上，并获得残差序列
model <- arima(time_series, order = c(8, 0, 3))  # ARMA(8,3)模型示例
residuals <- residuals(model)


# 2. 提取残差序列的三阶滞后数据
p <- 3  # 滞后阶数
residuals_lagged <- residuals[(p + 1):length(residuals)]
for (i in 1:p) {
  residuals_lagged <- cbind(residuals_lagged, residuals[(p - i + 1):(length(residuals) - i)])
}

# 打印残差序列的三阶滞后数据
print(residuals_lagged)


# 假设 residuals 是你要导出的数据
write.csv(residuals, file = "C:\\主修\\大四下\\毕业论文\\table\\111\\residuals11.csv", row.names = FALSE)

# 提取AR部分的数据
# 假设你已经拟合了一个ARMA(8,3)模型到时间序列上，名为model

# 提取模型的滞后8阶数据
lags <- 8
lagged_data <- NULL
for (i in 1:lags) {
  lagged_data <- cbind(lagged_data, lag(residuals(model), i))
}

# 去掉结果中的NA值
lagged_data <- lagged_data[-(1:lags), ]

# 打印滞后8阶数据
print(lagged_data)
lagged_data <- lagged_data[1:479,];lagged_data
residuals_lagged <- residuals_lagged[1:479,];residuals_lagged
X <- cbind(lagged_data, residuals_lagged)
Y <- ts_data[1:479,];Y




# 有两个时间序列数据存储在向量中，名为ts_data和residuals_ts
# 将ts_data滞后8阶，residuals_ts滞后3阶

# 准备数据
n <- min(length(ts_data), length(residuals_ts)) # 取两个时间序列长度的最小值
p_A <- 8 # A 的滞后阶数
p_B <- 3 # B 的滞后阶数

# 构建滞后特征矩阵
X_A <- matrix(0, n - p_A, p_A) # 滞后特征矩阵 A
X_B <- matrix(0, n - p_B, p_B) # 滞后特征矩阵 B
for (i in 1:(n - p_A)) {
  X_A[i, ] <- ts_data[(i + (p_A - 1)):(i)] # 提取 A 的滞后值
}
for (i in 1:(n - p_B)) {
  X_B[i, ] <- residuals_ts[(i + (p_B - 1)):(i)] # 提取 B 的滞后值
}
n_selected <- 479
X_A <- X_A[1:n_selected, ];X_A
X_B <- X_B[1:n_selected, ];X_B

# 合并自变量数据集
X <- cbind(X_A, X_B);X
Y <- ts_data[1:n_selected,];Y

##导入处理好的ARMA数据
data <- read.csv("C:/主修/大四下/毕业论文/table/111/ARMA(8,3).csv")
data
# 响应变量&自变量
Y <- as.matrix(data[,1]);Y
X <- as.matrix(data[,-1]);X



## 拟合Lasso模型
# 安装并加载 glmnet 包
install.packages("glmnet")
library(glmnet)

# 安装并加载 glmnet 包
install.packages("glmnet")
library(glmnet)

# 将第九个自变量的系数固定为1
penalty_factors <- rep(1, ncol(X))
penalty_factors[9] <- 0  # 将第九个自变量的惩罚系数设为0，表示不对其进行惩罚

# 使用 glmnet 包进行 Lasso 回归
fit <- cv.glmnet(X, Y, alpha = 1)  # alpha=1 表示 Lasso 回归

# 打印最佳的 lambda 值
best_lambda <- fit$lambda.min
print(best_lambda)

# 获取拟合的模型系数
coef_values <- coef(fit, s = best_lambda)
print(coef_values)


##自适应Lasso

# 定义自适应权重向量计算函数

calculate_weight <- function(beta_hat) {
  return(1/sqrt(abs(beta_hat)))
}
ones_column <- rep(1, nrow(X))

# 将这列数据与原始的自变量数据框合并起来
X <- cbind(ones_column, X)


# 自适应Lasso模型
adaptive_lasso <- function(X, Y, best_lambda) {
  # 计算样本数和自变量数
  #n <- nrow(X)
  #p <- ncol(X)
  
  # 拟合Lasso模型
  #fit <- glmnet::glmnet(X, Y, alpha = 1, lambda = best_lambda, standardize = FALSE)
  
  # 获取参数估计值的绝对值
  #beta_hat <- as.numeric(coef(fit))
  
  # 计算自适应权重向量
  w <- calculate_weight(beta_hat)
  # 将向量中的 Inf 替换为 0
  w[is.infinite(w)] <- 0
  
  # 定义目标函数
  objective_function <- function(beta, X, Y, w, lambda) {
    # 计算最小二乘估计值
    least_squares <- sum((Y - X %*% beta)^2)
  
    # 计算惩罚项
    abs_beta <- abs(beta)
    abs_beta[is.infinite(abs_beta)] <- 0
    penalty <- lambda * (w %*% abs(beta))
  
    # 返回目标函数值
    return(least_squares + penalty)
  }
  
  # 设置初始参数估计值
  # 创建一个包含13行的数据框，其中一列的值都为0
  initial_beta_guess <- data.frame(Column_of_Zeros = rep(0, 12))
  # 将 initial_beta_guess 转换为数值向量
  initial_beta_guess <- as.numeric(unlist(initial_beta_guess))
  
  # 使用optim函数找到最小值
  result <- optim(par = initial_beta_guess, fn = objective_function, X = X, Y = Y, lambda = best_lambda, w = w, method = "L-BFGS-B", lower = -1, upper = 1)
   
  return(result$par)
}

# 通过交叉验证选择最佳的lambda值
cv_fit <- glmnet::cv.glmnet(X, Y, alpha = 1, standardize = FALSE)
best_lambda <- cv_fit$lambda.min

# 根据最佳lambda值拟合自适应Lasso模型
final_beta <- adaptive_lasso(X, Y, best_lambda)

# 输出最佳lambda值
print(paste("Best lambda value:", best_lambda))

# 输出拟合方程的系数
print("Fitted coefficients:")
print(final_beta)

# 通过交叉验证选择最佳的lambda值
cv_fit <- glmnet::cv.glmnet(X, Y, alpha = 1, standardize = FALSE)
best_lambda <- 5.3918

# 根据最佳lambda值拟合自适应Lasso模型
final_beta <- adaptive_lasso(X, Y, best_lambda)

# 输出最佳lambda值
print(paste("Best lambda value:", best_lambda))

# 输出拟合方程的系数
print("Fitted coefficients:")
print(final_beta)


## 改进后的Adaptive Lasso
# 响应变量&自变量
Y <- as.matrix(data[,1]);Y
X <- as.matrix(data[,-1]);X


# 令lambda=3.2911和gamma2=0.5，需要用5折交叉验证方法求出gamma3的值

# 定义惩罚权重的方程,其中j为1~13,gamma为需要argmin的数
#calculate_weight_MA <- function(j,beta_hat) {
  #return(j^0.1/sqrt_abs_beta_hat)
#}
calculate_weight_MA <- function(j, gamma, sqrt_abs_beta_hat) {
  return(j^gamma/sqrt_abs_beta_hat)
}

# 将这列数据与原始的自变量数据框合并起来
# 创建一列全部为 1 的数据
#ones_column <- rep(1, nrow(X))

# 将这列数据与原始的自变量数据框合并起来
#X <- cbind(ones_column, X)
#用训练集定义ARMA的拟合方程
MA_adaptive_lasso <- function(X, Y) {
  p <- ncol(X)
  
  # 初始化参数估计值
  beta <- rep(0, ncol(X))
  #beta_hat <- as.numeric(c(0,-0.4166,-0.3751,-0.0906,0.0606,-0.0845,0.0056,-0.0918,1,0.0329,0.4037,0.5209))
  
  # 拟合Lasso模型
  fit <- glmnet::glmnet(X, Y, alpha = 1, lambda = lambda, standardize = FALSE)
  
  # 获取参数估计值的绝对值
  beta_hat <- as.numeric(coef(fit))
  j <- c(1,2,3,4,5,6,7,8,1,1,2,3)
  # 计算自适应权重向量
  sqrt_abs_beta_hat <- sqrt(abs(beta_hat))
  w <- calculate_weight_MA(j, gamma, sqrt_abs_beta_hat)
  # 将向量中的 Inf 替换为 0
  w[is.infinite(w)] <- 0
  # 定义目标函数
  objective_function <- function(w ,beta, X, Y) {
    # 计算最小二乘估计值
    least_squares <- sum((Y - X %*% beta)^2)
    #least_squares <- sum((Y - X %*% beta)^2)
  
    # 计算惩罚项
    penalty <- 3.2911 * sum(w * abs(beta))
  
    # 返回目标函数值
    return(least_squares + penalty)
  }
  # 为参数提供初始猜测值
  initial_beta_guess <- rep(0, ncol(X))
  # 设置参数
  gamma <- best_gamma
  # 使用 argmin 函数找到最小值
  # 设置参数取值范围和约束条件
  lower_bound <- -1
  upper_bound <- 1
  result <- optim(par = initial_beta_guess, fn = objective_function, X = X, Y = Y, method = "L-BFGS-B", lower = -1, upper = 1)

  #result <- optim(par = initial_beta_guess, fn = objective_function, X = X, Y = Y, model_params = model_params, method = "L-BFGS-B", lower = -1, upper = 1)
  beta <- result$par  
  return(beta)
}

# 自变量
# 按行分割矩阵成 5 个子集
X <- as.matrix(data[,-1]);X
num_rows <- nrow(X)
subset_size <- num_rows %/% 5  # 每个子集的大小
remainder <- num_rows %% 5  # 余数，用于处理行数无法整除 5 的情况

# 初始化一个列表，用于存储子集
x <- list()

start_row <- 1
for (i in 1:5) {
  end_row <- start_row + subset_size - 1
  if (remainder > 0) {
    end_row <- end_row + 1  # 将余数分配给前面的子集
    remainder <- remainder - 1
  }
  x[[i]] <- X[start_row:end_row, ]
  start_row <- end_row + 1
}

# 输出子集
for (i in 1:5) {
  print(x[[i]])
}

# 因变量
Y <- as.matrix(data[,1]);Y
num_rows <- nrow(Y)
subset_size <- num_rows %/% 5  # 每个子集的大小
remainder <- num_rows %% 5  # 余数，用于处理行数无法整除 5 的情况

# 初始化一个列表，用于存储子集
y <- list()

start_row <- 1
for (i in 1:5) {
  end_row <- start_row + subset_size - 1
  if (remainder > 0) {
    end_row <- end_row + 1  # 将余数分配给前面的子集
    remainder <- remainder - 1
  }
  y[[i]] <- Y[start_row:end_row, ]
  start_row <- end_row + 1
}

# 输出子集
for (i in 1:5) {
  print(y[[i]])
}



# 设置要尝试的 gamma 值范围
gamma_values <- seq(0.1, 10, by = 0.1)

# 初始化交叉验证误差向量
cv_errors <- numeric(length(gamma_values))

# 循环遍历每个 gamma 值进行交叉验证
for (i in seq_along(gamma_values)) {
  gamma <- gamma_values[i]
  fold_errors <- numeric()
  
  # 循环进行 5 折交叉验证
  for (i in 1:5) {
    # 将第 i 个子集作为测试集
    test_set_x <- x[[i]]
    test_set_y <- y[[i]]
  
    # 将其余四个子集合并为训练集
    train_set_x <- do.call(rbind, x[-i])
    train_set_y <- do.call(rbind, y[-i])
    # 使用测试集的解释变量和拟合得到的 beta 值来估计被解释变量
    estimated_y <- test_set_x %*% beta
  
    # 计算测试误差
    fold_errors <- sqrt(mean((estimated_y - test_set_y)^2))
  
    # 将测试误差存储到测试误差列表中
    fold_errors <- c(fold_errors, fold_errors)
  }
  
  # 计算平均交叉验证误差
  cv_errors[i] <- mean(fold_errors)
}

# 找到最小交叉验证误差对应的 gamma 值
best_gamma <- gamma_values[which.min(cv_errors)]

# 输出最优的 gamma 值
print(best_gamma)

# 根据最佳lambda值拟合自适应Lasso模型
final_beta <- MA_adaptive_lasso(X, Y)

# 输出拟合方程的系数
print("Fitted coefficients:")
print(final_beta)





