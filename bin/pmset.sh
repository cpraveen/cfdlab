# To restore default values
#    sudo pmset restoredefaults
#
# hibernatemode
#   3 : copy to disk, RAM is on; fast wake
#   25: copy to disk, RAM is off; slow to wake, good for battery
sudo pmset -c hibernatemode 3
sudo pmset -b hibernatemode 3

# Default = 1
sudo pmset -a standby       1
sudo pmset -c standbydelay  1800 # seconds
sudo pmset -b standbydelay  900  # seconds

# Dont wake for network access
sudo pmset -a tcpkeepalive  0

# On arm cpus, disabling powernap can be bad, see
# https://www.bravolt.com/post/why-won-t-my-computer-sleep
sudo pmset -a powernap      0

sudo pmset -b lowpowermode  1

sudo pmset -c displaysleep  5 # minutes
sudo pmset -b displaysleep  2 # minutes

#sudo pmset -c disksleep     10 # minutes
#sudo pmset -b disksleep     5  # minutes

sudo pmset -a womp          0
sudo pmset -a proximitywake 0

sudo pmset -c sleep         1 # minutes
sudo pmset -b sleep         1 # minutes

sudo pmset -a disablesleep  0

# For intel machines which have discrete gpu
#sudo pmset -c gpuswitch     2
#sudo pmset -b gpuswitch     0
