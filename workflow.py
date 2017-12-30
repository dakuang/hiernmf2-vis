import sys
import os

data_file = sys.argv[1]
num_leafs = int(sys.argv[2])

print '[INFO] Generating term-document matrix ......'
os.system('python txt2mtx.py {0:s}'.format(data_file))

print '[INFO] Applying TF-IDF ......'
os.system('matlab -nojvm -nodesktop -nosplash -r "mtx2mat(\'{0:s}\'); exit"'.format(data_file))

print '[INFO] Run HierNMF2 ......'
os.system('matlab -nojvm -nodesktop -nosplash -r "run_hiernmf2(\'{0:s}\', {1:d}); exit"'.format(data_file, num_leafs))

print '[INFO] Generating hierarchy display ......'
#os.system('source ~/software/python/bin/activate')
os.system('python generate_qtree_priority_simplified.py --alg=1 --max_cluster={0:d} --dataset={1:s} --width=44in --height=34in hier8_data_{2:s}_1.mat'.format(num_leafs, data_file, data_file))
os.system('mv hier8_data_{0:s}_1.pdf {1:s}_{2:d}.pdf'.format(data_file, data_file, num_leafs))
#os.system('deactivate')
