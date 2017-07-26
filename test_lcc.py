import os

os.chdir('/home/owht/KIZ/software/tophat-2.1.1.Linux_x86_64')
os.system('./tophat -o ./test_out_lcc_2 ./test_data/test_ref ./test_data/reads_1.fq ./test_data/reads_2.fq')
