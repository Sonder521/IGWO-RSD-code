'''
author:Li Shijjie
2023/7/14
基于IGWO优化lidar多通道拼接RSD三目标

'''
import random
import numpy
from numbers import Real
from telnetlib import BRK
from tkinter import image_names
from turtle import title
import pandas as pd
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
import matplotlib.animation as animation
import math
from netCDF4 import Dataset
from itertools import combinations
from copy import deepcopy
import time

def function_R(x,y):
    
    A1=pd.Series(x)
    B1=pd.Series(y)
    value=B1.corr(A1,method='pearson')
    return value 
def function_S(x,y,h):

    if (len(x)%2)!=0:
        x=x[0:len(x)-1]
        y=y[0:len(y)-1]
    X1=pd.Series(x)
    Y1=pd.Series(y)
    h =pd.Series(h)

    part = len(X1)/2
    part = int(part) #变成整数
    #下半部分拟合
    X1_down = X1[0:part-1]
    Y1_down = Y1[0:part-1]
    h_down = h[0:part-1]
    #X1_down = sm.add_constant(X1_down) # 向矩阵 X 添加截距列（x0=[1,...1]）
    model = sm.OLS(Y1_down, X1_down) # 建立最小二乘模型（OLS）
    results = model.fit() # 返回模型拟合结果
    
    #print('下',ers[0],'上',ers[1])   #这里求的是k的下限和上限

    #获得残差
    res = results.resid
    np.asarray(res)
    np.asarray(h_down)
    model1 = sm.OLS(res,h_down) # 残差和高度拟合
    results1 = model1.fit() # 返回模型拟合结果

    par_down = results1.params
    k_down = par_down[0]
    #print('k=',k_down)
    ers_down = results1.conf_int(0.05) #系数置信区间 
    dertak_down = (ers_down[1]-k_down)/2
    #得到K_down\dertak_down

    #上半部分拟合
    X1_up = X1[part:len(X1)-1]
    Y1_up = Y1[part:len(X1)-1]
    h_up = h[part:len(X1)-1]
    #X1_up = sm.add_constant(X1_up) # 向矩阵 X 添加截距列（x0=[1,...1]）
    model = sm.OLS(Y1_up, X1_up) # 建立最小二乘模型（OLS）
    results = model.fit() # 返回模型拟合结果
    
    #print('下',ers[0],'上',ers[1])   #这里求的是k的下限和上限

    #获得残差
    res = results.resid
    np.asarray(res)
    np.asarray(h_up)
    model1 = sm.OLS(res,h_up) # 残差和高度拟合
    results1 = model1.fit() # 返回模型拟合结果

    par = results1.params
    k_up = par[0]
    #print('k=',k_up)
    ers = results1.conf_int(0.05) #系数置信区间 
    dertak_up = (ers[1]-k_up)/2

    value = abs(k_down-k_up)/(dertak_down[0]**2+dertak_up[0]**2)**0.5
    #print('S=',value)
    return -value               #能理解
def function_D(X0,Y0,x1,y1):#计算区域总偏差

    np.asarray(X0)
    np.asarray(Y0)
    #做拟合
    model = sm.OLS(Y0, X0) # 建立最小二乘模型（OLS）
    results = model.fit() # 返回模型拟合结果
    k0 = results.params     #拟合斜率
    #print(k0[0])
    all_delta_Y1 = 0
    # k0 = 47
    for i in range(0,len(x1)):
        Y1_fit = x1[i]*k0[0]
        
        delta_Y1 = abs(Y1_fit-y1[i])/y1[i]
        all_delta_Y1 = all_delta_Y1 + delta_Y1
    value = all_delta_Y1/len(y1)
    float(value)
    return -abs(value),k0[0]
def F_objective(min_x,max_x,N):#随机产生拼接上下限

    pop_size = N
    z0 = [random.randint(min_x,max_x) for i in range(0,pop_size)] 
    z1 = [random.randint(min_x,max_x) for i in range(0,pop_size)] 

    for i in range(0,pop_size):             #拼接长度设置为不小于10,且保证z0小于z1
        if z0[i]>z1[i]:
            huan=0
            huan = z0[i]
            z0[i] = z1[i]
            z1[i] = huan+10
        elif z0[i]>z1[i]-10 and z0[i]<=z1[i]:
            z1[i]+=10

    output = np.dstack((z0,z1))
    output = output.squeeze()
    D = 2
    MaxValue = np.ones((1, D))*max_x
    MinValue = np.ones((1, D))*min_x
    Boundary = np.array([MaxValue, MinValue])

    return output
#函数

def F_RSD(solution,CH1data1,CH2data1):

    CH1data = CH1data0[solution[0]:solution[1]]    
    CH2data = CH2data0[solution[0]:solution[1]]
    h_raw1 = h_raw[solution[0]:solution[1]]
            
    function1_values=function_R(CH2data,CH1data)
    function2_values=function_S(CH2data,CH1data,h_raw1)
    function3_values=function_D(CH2data,CH1data,CH2data1,CH1data1)
        
    FV = -(0.3380*function1_values+0.3139*function2_values+0.3481*function3_values[0])
    FunctionValue = [FV,function3_values[1],function1_values,function2_values,function3_values[0]]

    return FunctionValue
def IGWO_linear(objf, lb, ub, dim, SearchAgents_no, Max_iter):

    # 初始化 alpha, beta, and delta_pos
    Alpha_pos = numpy.zeros(dim)  # 位置.形成30的列表
    Alpha_score = float("inf")  # 这个是表示“正负无穷”,所有数都比 +inf 小；正无穷：float("inf"); 负无穷：float("-inf")

    Beta_pos = numpy.zeros(dim)
    Beta_score = float("inf")

    Delta_pos = numpy.zeros(dim)
    Delta_score = float("inf")  # float() 函数用于将整数和字符串转换成浮点数。

    Convergence_curve = numpy.zeros(Max_iter)
    
    Positions = numpy.zeros((SearchAgents_no, dim))
    Positions = F_objective(lb,ub,SearchAgents_no)

    CH1data1 = CH1data0[lb:ub]    #这里应该是z0[i]
    CH2data1 = CH2data0[lb:ub] 
    GWO_op_solution=[]
    #迭代寻优
    for l in range(0, Max_iter):  # 迭代100
        for i in range(0, SearchAgents_no):  # 5
            # 返回超出搜索空间边界的搜索代理

            for j in range(dim):  # 30
                Positions[i, j] = numpy.clip(Positions[i, j], lb, ub)  # clip这个函数将将数组中的元素限制在a_min(-100), a_max(100)之间，大于a_max的就使得它等于 a_max，小于a_min, 的就使得它等于a_min。
            #对positions里的数组排序，确保起点小于终点
            if Positions[i][0]>Positions[i][1]:
                huan=0
                huan = Positions[i][0]
                Positions[i][0] = Positions[i][1]
                Positions[i][1] = huan+10
            elif Positions[i][0]>Positions[i][1]-10 and Positions[i][0]<=Positions[i][1]:
                Positions[i][1]+=10

            # 计算每个搜索代理的目标函数
            fitness_all = objf(Positions[i, :],CH1data1,CH2data1)  # 把某行数据带入函数计算；fitness是一个值,objf这里指F_RSD
            fitness=fitness_all[0]-3
            # print("经过计算得到：",fitness)

            # Update Alpha, Beta, and Delta
            if fitness < Alpha_score:
                Alpha_score = fitness  # Update alpha
                Alpha_pos = Positions[i, :].copy()

            if (fitness > Alpha_score and fitness < Beta_score):
                Beta_score = fitness  # Update beta
                Beta_pos = Positions[i, :].copy()

            if (fitness > Alpha_score and fitness > Beta_score and fitness < Delta_score):
                Delta_score = fitness  # Update delta
                Delta_pos = Positions[i, :].copy()

        # 以上的循环里，Alpha、Beta、Delta

        a = 2 - l * ((2) / Max_iter);  #   a从2线性减少到0
        # a = 2*math.exp(-l/Max_iter) #a改为指数下降，增强全局能力

        for i in range(0, SearchAgents_no):#种群个数
            for j in range(0, dim):
                r1 = random.random()  # r1 is a random number in [0,1]主要生成一个0-1的随机浮点数。
                r2 = random.random()  # r2 is a random number in [0,1]

                A1 = 2 * a * r1 - a;  # Equation (3.3)
                C1 = 2 * r2;  # Equation (3.4)
                # D_alpha表示候选狼与Alpha狼的距离
                D_alpha = abs(C1 * Alpha_pos[j] - Positions[
                    i, j]);  # abs() 函数返回数字的绝对值。Alpha_pos[j]表示Alpha位置，Positions[i,j])候选灰狼所在位置
                X1 = Alpha_pos[j] - A1 * D_alpha;  # X1表示根据alpha得出的下一代灰狼位置向量

                r1 = random.random()
                r2 = random.random()

                A2 = 2 * a * r1 - a;  #
                C2 = 2 * r2;

                D_beta = abs(C2 * Beta_pos[j] - Positions[i, j]);
                X2 = Beta_pos[j] - A2 * D_beta;

                r1 = random.random()
                r2 = random.random()

                A3 = 2 * a * r1 - a;
                C3 = 2 * r2;

                D_delta = abs(C3 * Delta_pos[j] - Positions[i, j]);
                X3 = Delta_pos[j] - A3 * D_delta;

                Positions[i, j] = (X1 + X2 + X3) / 3  # 候选狼的位置更新为根据Alpha、Beta、Delta得出的下一代灰狼地址。
        # print(Positions)
        Convergence_curve[l] = Alpha_score

        Alpha_can = F_RSD(Alpha_pos,CH1data1,CH2data1)

        if (l % 1 == 0):
            print(['迭代次数为' + str(l) + ' 的迭代结果' + str(Alpha_can)+' 位置是'+str(Alpha_pos)]);  # 每一次的迭代结果
        op_solution=[Alpha_can[0],Alpha_can[1],Alpha_can[2],Alpha_can[3],Alpha_can[4],Alpha_pos[0],Alpha_pos[1]]
        GWO_op_solution.append(op_solution)
    # final_solution=[Alpha_can[0],Alpha_can[1],Alpha_can[2],Alpha_can[3],Alpha_can[4],Alpha_pos[0],Alpha_pos[1]]
    return GWO_op_solution

def get_ub(data):
    bg = np.mean(data[1800:2000])
    SNR = []
    for i in range(len(data)):
        SNRi = (data[i]-bg)/math.sqrt(data[i])
        SNR.append(SNRi)
    SNR = SNR[10:i]
    up = next(k for k, value in enumerate(SNR)if value < 10)
    return up+10           #返回拼接高度上限
#主程序

nc_obj = Dataset('ahsrl_day_20131016T1900_20131016T2000_3600s_120m.nc')
CH1data0 = (nc_obj.variables['combined_counts_hi'][0][1:2000])
#print(nc_obj.variables['combined_counts_hi'])
#print(len(CH1data0))
CH2data0 = (nc_obj.variables['combined_counts_lo'][0][1:2000])
h_raw = (nc_obj.variables['raw_range'][1:2000])
ub = get_ub(CH2data0)
func_details = ['F_RSD', 150, ub, 2]
function_name = func_details[0]
Max_iter = 100#迭代次数
lb = func_details[1]#下界
ub = func_details[2]#上届
dim = func_details[3]#狼的寻值范围
SearchAgents_no = 20#寻值的狼的数量

# GWO_op_solution = []
# # writer = pd.ExcelWriter('D:/A正事专用文件夹/硕士/瑞利雷达/多通道拼接技术/GWO-lidar/2.20GWO尖峰最优pareto面.xls')
# G = globals()
# for i in range(0,3):
#     print('迭代次数',i+1)
#     op_solution = GWO(F_RSD, lb, ub, dim, SearchAgents_no, Max_iter)
#     GWO_op_solution.append(op_solution)  
start = time.perf_counter()
IGWO_op_solution = IGWO_linear(F_RSD, lb, ub, dim, SearchAgents_no, Max_iter)
end = time.perf_counter()
print('Running time: %s Seconds'%(end-start))
# GWO_op_solution1 = GWO_exp(F_RSD, lb, ub, dim, SearchAgents_no, Max_iter)
#     G[str(i)] = pd.DataFrame(GWO_op_solution)
#     G[str(i)].to_excel(writer, sheet_name=str(i),index=0)
# writer.save()
# writer.close()
# x=range(len(GWO_op_solution))
# y=[x[0] for x in GWO_op_solution]
# y1=[x[0] for x in GWO_op_solution1]
# plt.plot(x,y)
# plt.plot(x,y1)
# plt.xlabel('迭代次数')
# plt.ylabel('目标')
# plt.rcParams['font.sans-serif']=['SimHei']
# plt.rcParams['axes.unicode_minus'] = False
# plt.legend(['linear','exp'])
# plt.show()

# book = xlwt.Workbook(encoding='utf-8',style_compression=0)
# sheet = book.add_sheet('GWO',cell_overwrite_ok=True)

# col = ('次数','K','拟合系数','R','S','D','拼接下限','拼接上限')
# for i in range(0,8):
# 	sheet.write(0,i,col[i])

# for i in range(0,len(GWO_op_solution)):
#     sheet.write(i+1,0,i+1)
#     for j in range(0,len(GWO_op_solution[0])):
#         sheet.write(i+1,j+1,str(GWO_op_solution[i][j]))

# savepath = 'D:/A正事专用文件夹/硕士/瑞利雷达/多通道拼接技术/GWO-lidar/3.6-GWO迭代100次-改权重.xls'
# book.save(savepath)
