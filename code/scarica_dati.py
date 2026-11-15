"""Scarica i dati pubblici (circa 210 MB) nella cartella data/ accanto agli script."""
import os, urllib.request
HERE=os.path.dirname(os.path.abspath(__file__)); D=os.path.join(HERE,'data'); os.makedirs(D,exist_ok=True)
PP='https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/'
DE='https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/'
for name,url in [('pp.dat',PP+'Pantheon%2BSH0ES.dat'),('pp.cov',PP+'Pantheon%2BSH0ES_STAT%2BSYS.cov'),
                 ('2MPP_redshift.txt','https://media.githubusercontent.com/media/KSaid-1/pvhub/main/data/2MPP_redshift.txt'),
                 ('desi_mean.txt',DE+'desi_gaussian_bao_ALL_GCcomb_mean.txt'),('desi_cov.txt',DE+'desi_gaussian_bao_ALL_GCcomb_cov.txt')]:
    if os.path.exists(os.path.join(D,name)): print('gia\' presente:',name); continue
    print('scarico',name); urllib.request.urlretrieve(url,os.path.join(D,name))
# altri campioni di supernovae (Union3, DES 5 anni, DES ricalibrato), per SNSET=union3|desy5|dovekie
SN='https://raw.githubusercontent.com/CobayaSampler/sn_data/master/'
for f in ('Union3/lcparam_full.txt','Union3/mag_covmat.txt','DESY5/DES-SN5YR_HD.csv','DESY5/covsys_000.txt','DES-Dovekie/DES-Dovekie_HD.csv','DES-Dovekie/covtot_inv_000.npz'):
    p=os.path.join(D,'sn',*f.split('/')); os.makedirs(os.path.dirname(p),exist_ok=True)
    if os.path.exists(p): print('gia\' presente:',f); continue
    print('scarico',f); urllib.request.urlretrieve(SN+f,p)
# verosimiglianza Planck 'plik-lite' in Python (Prince & Dunkley), per cmb_full.py e mcmc_cmb.py
import zipfile, io, shutil
if not os.path.isdir(os.path.join(HERE,'planck-lite-py')):
    print('scarico planck-lite-py')
    z=zipfile.ZipFile(io.BytesIO(urllib.request.urlopen('https://github.com/heatherprince/planck-lite-py/archive/refs/heads/master.zip').read()))
    z.extractall(HERE); shutil.move(os.path.join(HERE,'planck-lite-py-master'),os.path.join(HERE,'planck-lite-py'))
print('fatto')
