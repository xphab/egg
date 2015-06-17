setwd("C:/Users/Daliang/Dropbox/ToolDevelop/github/egg")

x=c(1,1.1,2.3,2,3.4,1.5)
y=c(2,2,2,2,2,1.9)

# 1 # Normality test
# install.packages("nortest")
source(file = "Rcode/nor.test.r")
nt.x=nor.test(x)
nt.y=nor.test(y)

# 2 # Variance homogenity test
install.packages("car")
source(file = "Rcode/vh.test.r")
vht=vh.test(x,y)

# 3 # t-test
