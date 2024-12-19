import math
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
from matplotlib import pyplot as plt
#%matplotlib inline

# ------------------------------- PID制御関連 -------------------------------
# 関数の定義 ---------------------
# PID制御
def PID(kp, ki, kd, theta_goal, theta_current, error_sum, error_pre):
    error = theta_goal - theta_current# 偏差（error）を計算
    error_sum += error # 偏差の総和（積分）を計算
    error_diff = error-error_pre # PI制御からの追加：1時刻前の偏差と現在の偏差の差分（微分）を計算
    m = (kp * error) + (ki * error_sum) + (kd*error_diff) # 操作量を計算
    return m, error_sum, error

# 係数などの設定 --------------------
# 以下の変数の数値を変えると結果が変わる
kp = 0.1 # 比例ゲイン
ki = 0.01 #0.5 # 積分ゲインの値を大きくして，意図的に振動を発生させる
kd = 0.5 #0.5 # 微分ゲイン：急激な変化を抑える
theta_start = 0.0 # 初期角度
theta_goal = 90.0 # 目標角度
time_length = 150 #150 # 計測時間 
offset = 1.0 # 定常偏差

# 以下の変数は変えなくて良い
theta_current = theta_start # 目標角度
error_sum = 0.0 # 偏差の総和（積分）
error_pre = 0.0 # 1時刻前の偏差
time_list = [0] # 時刻のリスト（描画用）
theta_list = [theta_start] # 現在地のリスト（描画用）
animation_time_list = [time_list.copy()]
animation_theta_list = [theta_list.copy()]

# PID制御 -----------------------
for time in range(time_length):
    m, error_sum, error = PID(kp, ki, kd, theta_goal, theta_current, error_sum, error_pre) # 操作量を計算
    theta_current += m # 現在角度に操作量を足す（実際は，この操作量をもとにモータを動かす）
    theta_current -= offset
    error_pre = error # 一時刻前の偏差として保存しておく（D制御用）
    time_list.append(time) # 描画用
    theta_list.append(theta_current) # 描画用
    animation_time_list.append(time_list.copy())
    animation_theta_list.append(theta_list.copy())

# ------------------------- アニメーション関連 -------------------------------------
# PID制御のグラフ描画
def plot_pid_graph(ax, time, theta_goal, animation_time_list, animation_theta_list):
    ax.hlines([theta_goal], 0, max(animation_time_list[-1]), "red", linestyles='dashed') #ゴールを赤色の点線で表示
    ax.plot(animation_time_list[time], animation_theta_list[time], label="PID", color="blue") # PID制御のグラフを描画
    ax.set_xlim(-1, max(animation_time_list[-1])) # min=0の場合，グラフの左端が切れるので，min=-1に設定
    if max(animation_theta_list[-1]) < theta_goal: # 定常偏差によりtheta_goalよりもグラフが下回ってしまったらtheta_goalの赤い点線が見えるように範囲を設定
        ax.set_ylim(0, theta_goal+1) # 赤い点線が見えるように +1 している
    else:
        ax.set_ylim(0, max(animation_theta_list[-1])+1)# 赤い点線が見えるように +1 している
    ax.set_xlabel(r'$t$') 
    ax.set_ylabel(r'$\theta$') 
    ax.legend(loc='lower right') # 凡例を表示

# ロボットアームの描画
def plot_robot_arm(ax, time, animation_theta_list):
    ax.plot([0,0], [0,-1], color="black") # 固定リンク
    rad = math.radians(animation_theta_list[time][-1]-90) # 真下方向を0°とする
    x = math.cos(rad) # 順運動学
    y = math.sin(rad) # 順運動学
    ax.plot([0,x], [0,y], color="orange") # 稼働リンク
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    
# figureを作成
fig = plt.figure()
ax_pid = fig.add_subplot(1, 2, 1)
ax_arm = fig.add_subplot(1, 2, 2)
ax_pid.set_aspect("equal")#画像の比率を同じにする
ax_arm.set_aspect("equal")#画像の比率を同じにする
ax_arm.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False) # 軸のメモリを消す（ロボットアーム側はメモリの情報は不要なので）

# 各フレーム毎の描画処理
def update(time):
    ax_pid.cla()
    ax_arm.cla()
    plot_pid_graph(ax_pid, time, theta_goal, animation_time_list, animation_theta_list)
    plot_robot_arm(ax_arm, time, animation_theta_list)
    plt.suptitle('kp={}, ki={}, kd={}, offset={} \n \n t={}, theta={:.3g}'.format(kp, ki, kd, offset, animation_time_list[time][-1], animation_theta_list[time][-1]), x=0.5, y=0.90) # タイトルに時間と角度を表示
 
# アニメーション化
ani = FuncAnimation(fig, update, interval=50, frames=len(time_list))
HTML(ani.to_jshtml()) # HTMLに
#ani.save('pid.mp4', writer="ffmpeg") # mp4で保存．これを実行すると処理時間が増加します
