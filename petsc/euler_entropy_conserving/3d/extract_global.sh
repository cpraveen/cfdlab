echo "Reading log.txt and writing into global.txt"
grep "it,t,ke,ent=" log.txt |awk '{print $2, $3, $4, $5}' > global.txt
