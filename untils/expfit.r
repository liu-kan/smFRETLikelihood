plot(TauDdf)
# 通过这个散点图确定参数a, b的初值
a_start <- TauDdf[1,2] # 参数a是x = 0时y的取值
b_start<- 2*log(2)/a_start # b 是衰减速率
# 拟合模型
x<-as.vector(TauDdf$x)
y<-as.vector(TauDdf$y)
m <- nls(y ~ a*exp(-b*x), start = list(a = a_start, b = b_start))
print(cor(y, predict(m)))
lines(x, predict(m), col = "red", lty = 2, lwd = 3)
print(1/coef(m)[2])