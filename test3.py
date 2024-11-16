from matplotlib import pyplot as plt

# ------------------------------- PI制御関連 -------------------------------
# 関数の定義 ------------------------
# PI制御
def PI(kp, ki, theta_goal, theta_current, error_sum):
    error = theta_goal - theta_current# 偏差（error）を計算
    error_sum += error # 偏差の総和（積分）を計算
    m = (kp * error) + (ki * error_sum) # 操作量を計算
    return m, error_sum

# 係数などの設定 --------------------
kp = 0.1 # 比例ゲイン
ki = 0.5 # 積分ゲインの値を大きくして，意図的に振動を発生させる
theta_start = 0.0 # 初期角度
theta_goal = 90.0 # 目標角度
time_length = 200 # 計測時間 
theta_current = theta_start # 目標角度
error_sum = 0.0 # 偏差の総和（積分）
time_list = [0] # 時刻のリスト（描画用）
theta_list = [theta_start] # 現在地のリスト（描画用）

# 制御 ------------------------
for time in range(time_length):
    m, error_sum = PI(kp, ki, theta_goal, theta_current, error_sum) # 操作量を計算
    theta_current += m # 現在角度に操作量を足す（実際は，この操作量をもとにモータを動かす）
    time_list.append(time) # 描画用
    theta_list.append(theta_current) # 描画用
    
# 描画 --------------------------
plt.hlines([theta_goal], 0, time_length, "red", linestyles='dashed') #ゴールを赤色の点線で表示
plt.plot(time_list, theta_list, label="PI", color="black") # PI制御のグラフを描画

# ------------------------------- PID制御関連 -------------------------------
# 関数の定義 ---------------------
# PID制御
def PID(kp, ki, kd, theta_goal, theta_current, error_sum, error_pre):
    error = theta_goal - theta_current# 偏差（error）を計算
    error_sum += error # 偏差の総和（積分）を計算
    error_diff = error-error_pre # PI制御からの追加：1時刻前の偏差と現在の偏差の差分（微分）を計算
    m = (kp * error) + (ki * error_sum) + (kd*error_diff) # 操作量を計算
    return m, error_sum, error

# 係数などの設定 ------------------
kd = 0.5 # 微分ゲイン：急激な変化を抑える
error_pre = 0.0 # 1時刻前の偏差
# PI制御の時の数値を初期化
theta_start = 0.0; theta_current = theta_start; error_sum = 0.0; time_list = [0]; theta_list = [theta_start] 

# PID制御 -----------------------
for time in range(time_length):
    m, error_sum, error = PID(kp, ki, kd, theta_goal, theta_current, error_sum, error_pre) # 操作量を計算
    theta_current += m # 現在角度に操作量を足す（実際は，この操作量をもとにモータを動かす）
    error_pre = error # 一時刻前の偏差として保存しておく（D制御用）
    time_list.append(time) # 描画用
    theta_list.append(theta_current) # 描画用

# 描画
plt.plot(time_list, theta_list, label="PID", color="blue") # PID制御のグラフを描画
plt.xlabel(r'$t$') 
plt.ylabel(r'$\theta$') 
plt.legend(loc='lower right') # 凡例を表示
plt.show() # グラフの表示
