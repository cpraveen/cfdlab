# hibernatemode
#   3, safe sleep, copy to disk, RAM is on
#   25, copy to disk, RAM is off, slow to wake, good for battery
sudo pmset -c hibernatemode 3
sudo pmset -b hibernatemode 25

# Default = 1
sudo pmset -a standby       1
sudo pmset -a standbydelay  1800

# Dont wake for network access
sudo pmset -a tcpkeepalive  0

sudo pmset -a powernap      0

sudo pmset -b lowpowermode  1

sudo pmset -c displaysleep  5
sudo pmset -b displaysleep  2

sudo pmset -c disksleep     5
sudo pmset -b disksleep     2

sudo pmset -a womp          0
sudo pmset -a proximitywake 0

sudo pmset -a sleep         1
sudo pmset -a disablesleep  0

# For intel machines which have discrete gpu
#sudo pmset -c gpuswitch     2
#sudo pmset -b gpuswitch     0
