#pre-processed matrix
# test
y1<-as.matrix(wines[,1:6])
x1<-scale(y1)/sqrt(11)  #centering and normalizing

#PCA of the Data Tables
svd_x1<-svd(x1)
u1<-svd_x1$u
v1<-svd_x1$v
Gamma1<-diag(svd_x1$d)   #get u,v and singular value
alpha1<-1/(max(svd_x1$d)^2)
G1<-u1%*%Gamma1          #get the factor scores for the first table

#Creating the alpha weight factor 
#the same way as get each alpha, miss this step, we get A
a<-c(rep(0.241,6), rep(0.239,6), rep(0.275,6), rep(0.273,5),rep(0.307,6),rep(0.302,5), rep(0.417,4), rep(0.272,6), rep(0.264,5),rep(0.309,4))
diffa<-c(0.241,0.239,0.275,0.273,0.307,0.302,0.417,0.272,0.264,0.309)

#Generalized PCA of X
M<-diag(rep(1/12,12))
X<-as.matrix(wines[1:53])
X<-scale(X)/sqrt(11) #centering and normalizing
A<-diag(a)
sqrtM<-sqrt(M)
sqrtA<-sqrt(A)
delta<-svd(sqrtM%*%X%*%sqrtA)$d    #get the sigular value of X 
P<-solve(sqrtM)%*%(svd(sqrtM%*%X%*%sqrtA)$u)
Q<-solve(sqrtA)%*%(svd(sqrtM%*%X%*%sqrtA)$v)
F<-P%*%diag(delta)   #factor scores

#Partial Factor Scores
F1<-10*0.241*x1%*%Q[1:6,1:2] #get the partial factor score for table 1

#Determining the Importance of the Tables in the Compromise

#Supplementary Table
supx<-wines[54:57]
supx<-scale(supx)/sqrt(11)  #centering and normalizing the supplementary table
xsup<-supx/max(svd(supx)$d) # X_sup (whose first singular value is now equal to 1)
delta<-diag(delta[1:2])
Qsup<-t(xsup)%*%M%*%P[,1:2]%*%solve(delta)
Fsup<-xsup%*%Qsup

#bootstrap
set.seed(1)
sample<-sample(1:10,10,replace=TRUE)
#（接下来利用for loop可实现，就是多次模拟求均值）

#HMFA: Hierarchical Multiple Factor Analysis
table1<-wines[,1:6]
table1<-scale(table1)/sqrt(11)
table1<-table1/max(svd(table1)$d)

table2<-wines[,7:12]
table2<-scale(table2)/sqrt(11)
table2<-table2/max(svd(table2)$d)

table3<-wines[,13:18]
table3<-scale(table3)/sqrt(11)
table3<-table3/max(svd(table3)$d)

table4<-wines[,19:23]
table4<-scale(table4)/sqrt(11)
table4<-table4/max(svd(table4)$d)

table5<-wines[,24:29]
table5<-scale(table5)/sqrt(11)
table5<-table5/max(svd(table5)$d)   #normalizing and centering the first five table

table_first5<-cbind(table1,table2,table3,table4,table5) #get the "men-table"
number1<-max(svd(table_first5)$d)#the first singular value of this ‘men-table’ 

#then, assume we get the first singular value of the last five table is number2=2.169
number2<-2.169
number<-c(rep(number1^(-2),5),rep(number2^(-2),5))
alpha_HMFA<-diffa*number #get the weight alpha
a_HMFA<-c(rep(alpha_HMFA[1],6),rep(alpha_HMFA[2],6),rep(alpha_HMFA[3],6),rep(alpha_HMFA[4],5),rep(alpha_HMFA[5],6),rep(alpha_HMFA[6],5),rep(alpha_HMFA[7],4),rep(alpha_HMFA[8],6),rep(alpha_HMFA[9],5),rep(alpha_HMFA[10],4))
A_HMFA<-diag(a_HMFA)
Y_HMFA<-sqrt(M)%*%X%*%sqrt(A_HMFA)
delta_HMFA<-svd(Y_HMFA)$d    #get the sigular value  
P_HMFA<-solve(sqrtM)%*%(svd(Y_HMFA)$u)
F_HMFA<-P_HMFA%*%diag(delta_HMFA)   #factor scores of HMFA
