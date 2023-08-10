# Default = 3, 
# 25 is slow to wake
sudo pmset -a hibernatemode 3
#sudo pmset -c hibernatemode 3
#sudo pmset -b hibernatemode 25

# Default = 1
sudo pmset -a standby       1
sudo pmset -a standbydelay  1800

# Dont wake for network access
sudo pmset -a tcpkeepalive  0

sudo pmset -a powernap      0

sudo pmset -b lowpowermode  1

# For intel machines which have discrete gpu
#sudo pmset -c gpuswitch     2
#sudo pmset -b gpuswitch     0
