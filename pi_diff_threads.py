import subprocess

if __name__ == '__main__':
    for i in range(1,999999999):
        thread_set_file = open('host','w')
        thread_set_file.write('75.pi slots=' + str(i//2) +'\n76.pi slots='+ str(i-i//2)+'\n')
        thread_set_file.close()
        subprocess.call('(time mpirun -np 1 -hostfile host ./a.out 1000000000 ' + str(i) + ') 22>> ./log.data', shell=True)
     
