"""Downloads the public data (about 210 MB) into the data/ folder next to the scripts."""
import os, urllib.request
HERE=os.path.dirname(os.path.abspath(__file__)); D=os.path.join(HERE,'data'); os.makedirs(D,exist_ok=True)
PP='https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/'
DE='https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/'
for name,url in [('pp.dat',PP+'Pantheon%2BSH0ES.dat'),('pp.cov',PP+'Pantheon%2BSH0ES_STAT%2BSYS.cov'),
                 ('2MPP_redshift.txt','https://media.githubusercontent.com/media/KSaid-1/pvhub/main/data/2MPP_redshift.txt'),
                 ('desi_mean.txt',DE+'desi_gaussian_bao_ALL_GCcomb_mean.txt'),('desi_cov.txt',DE+'desi_gaussian_bao_ALL_GCcomb_cov.txt')]:
    if os.path.exists(os.path.join(D,name)): print('already present:',name); continue
    print('downloading',name); urllib.request.urlretrieve(url,os.path.join(D,name))
# other supernova samples (Union3, DES 5-year, recalibrated DES), for SNSET=union3|desy5|dovekie
SN='https://raw.githubusercontent.com/CobayaSampler/sn_data/master/'
for f in ('Union3/lcparam_full.txt','Union3/mag_covmat.txt','DESY5/DES-SN5YR_HD.csv','DESY5/covsys_000.txt','DES-Dovekie/DES-Dovekie_HD.csv','DES-Dovekie/covtot_inv_000.npz'):
    p=os.path.join(D,'sn',*f.split('/')); os.makedirs(os.path.dirname(p),exist_ok=True)
    if os.path.exists(p): print('already present:',f); continue
    print('downloading',f); urllib.request.urlretrieve(SN+f,p)
# Planck 'plik-lite' likelihood in Python (Prince & Dunkley), for cmb_full.py and mcmc_cmb.py
import zipfile, io, shutil
if not os.path.isdir(os.path.join(HERE,'planck-lite-py')):
    print('downloading planck-lite-py')
    z=zipfile.ZipFile(io.BytesIO(urllib.request.urlopen('https://github.com/heatherprince/planck-lite-py/archive/refs/heads/master.zip').read()))
    z.extractall(HERE); shutil.move(os.path.join(HERE,'planck-lite-py-master'),os.path.join(HERE,'planck-lite-py'))
print('done')
