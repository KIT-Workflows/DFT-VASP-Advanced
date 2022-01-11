# coding: utf-8

import os, yaml, subprocess, math, shutil
import numpy as np
import pymatgen
import pymatgen.symmetry.bandstructure

from ast import literal_eval

from pymatgen.core.structure import Structure, Lattice
from pymatgen.io.vasp import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.outputs import Xdatcar

def line(array):
    l='   '
    for el in array: l+='% .14f   '%(el)
    return l.rstrip()

def bool2vasp(b):
    if b == True or b == 'True': return '.TRUE.'
    elif b == False or b == 'False': return '.FALSE.'

def minimal_incar_dict():
    global wano_file
    incar_dict={'NELM' : 1}

    scf_dict = wano_file['Electronic structure settings']

    incar_dict['ENCUT'] = math.ceil(0.6*scf_dict['Energy cutoff [eV]'])
    incar_dict['ALGO'] = scf_dict['SCF algorithm']
    incar_dict['LREAL'] = scf_dict['LREAL']
    incar_dict['ISPIN'] = scf_dict['ISPIN']
    incar_dict['ISMEAR'] = scf_dict['ISMEAR']
    incar_dict['SIGMA'] = scf_dict['SIGMA']

    if scf_dict['SOC']: 
        incar_dict['ISYM'] = 0

    return incar_dict

def electronic_incar_dict():
    global wano_file
    incar_dict={'LWAVE': bool2vasp(wano_file['Output files']['Keep WAVECAR']), 'LCHARG': bool2vasp(wano_file['Output files']['Keep CHGCAR'])}

    scf_dict = wano_file['Electronic structure settings']

    if scf_dict['Type of functional'] == 'VdW':
        func=scf_dict['VdW functional']
        incar_dict['LASPH'] = '.TRUE.'
        incar_dict['LUSE_VDW'] = '.TRUE.'

        if func == 'SCAN+rVV10':
            incar_dict['METAGGA'] = 'SCAN'
            incar_dict['BPARAM'] = 6.3
            incar_dict['CPARAM'] = 0.0093

        else:
            incar_dict['AGGAC'] = 0.0
            if func == 'optPBE-vdW':
                incar_dict['GGA'] = 'OR'

            elif func == 'optB88-vdW':
                incar_dict['GGA'] = 'BO'
                incar_dict['PARAM1'] = 0.1833333333
                incar_dict['PARAM2'] = 0.2200000000

            elif func == 'optB86b-vdW':
                incar_dict['GGA'] = 'MK'
                incar_dict['PARAM1'] = 0.1234
                incar_dict['PARAM2'] = 1.0000

            elif func == 'vdW-DF2':
                incar_dict['GGA'] = 'ML'
                incar_dict['Zab_vdW'] = -1.8867

            elif func == 'rev-vdW-DF2':
                incar_dict['GGA'] = 'ML'
                incar_dict['Zab_vdW'] = -1.8867
                incar_dict['PARAM1'] = 0.1234
                incar_dict['PARAM2'] = 0.711357

    else:
        if scf_dict['Type of functional'] == 'Hybrid':
            hybrid_func = scf_dict['Hybrid functional']

            func = 'B3' if hybrid_func == 'B3LYP' else 'PE'

            incar_dict['LHFCALC'] = '.TRUE.'

            if hybrid_func == 'HSE03':
                incar_dict['HFSCREEN'] = 0.3

            elif hybrid_func == 'HSE06':
                incar_dict['HFSCREEN'] = 0.2

            elif hybrid_func == 'B3LYP':
                incar_dict['AEXX'] = 0.20
                incar_dict['AGGAX'] = 0.72
                incar_dict['AGGAC'] = 0.81
                incar_dict['ALDAC'] = 0.20
        else:
            if scf_dict['Functional'] != 'LDA':
                func_dict = {'PW91': '91', 'PBE': 'PE', 'rPBE': 'RP', 'revPBE': 'RE', 'PBEsol': 'PS', 'VWN': 'VW'}
                func = func_dict[scf_dict['Functional']]

        incar_dict['GGA'] = func

        if scf_dict['Use vdW correction']:
            incar_dict['ADDGRID'] = '.FALSE.'
            incar_dict['LASPH'] = '.FALSE.'
            disp_dict={'D2':1,'D3':11,'D3BJ':12,'TS':2,'TS-SCS':2,'TS-HP':'21','MBD@rsSCS':202,'dDsC':4}
            incar_dict['IVDW'] = disp_dict[scf_dict['Type of correction']]

            if scf_dict['Type of correction'] == 'TS-SCS':
                incar_dict['LVDWSCS'] = '.TRUE.'
                incar_dict['LSCSGRAD'] = '.TRUE.'
                incar_dict['LSCALER0'] = '.TRUE.'

        else:
            incar_dict['LASPH'] = '.TRUE.'

    incar_dict['ENCUT'] = scf_dict['Energy cutoff [eV]']
    incar_dict['ALGO'] = scf_dict['SCF algorithm']
    incar_dict['LREAL'] = scf_dict['LREAL']
    incar_dict['PREC'] = scf_dict['Precision mode']
    incar_dict['NELMIN'] = scf_dict['Minimum number of SCF cycles']
    incar_dict['NELM'] = scf_dict['Maximum number of SCF cycles']
    incar_dict['EDIFF'] = scf_dict['EDIFF']
    incar_dict['ISPIN'] = scf_dict['ISPIN']
    incar_dict['ISMEAR'] = scf_dict['ISMEAR']
    incar_dict['SIGMA'] = scf_dict['SIGMA']

    if scf_dict['SOC']: 
        incar_dict['LASPH'] = '.TRUE.'
        incar_dict['LSORBIT'] = '.TRUE.'
        incar_dict['GGA_COMPAT'] = '.TRUE.'
        incar_dict['SAXIS'] = '0 0 1'
        incar_dict['ISYM'] = 0

    if scf_dict['Dipole correction to forces/potentials']:
        incar_dict['DIPOL'] = '.TRUE.'

    if scf_dict['Dipole correction to energy']:
        idipol_dict = {'a': 1, 'b': 2, 'c': 3, 'all': 4}
        incar_dict['IDIPOL'] = idipol_dict[scf_dict['Direction']]

    incar_dict['IBRION'] = -1

    return incar_dict

def dos_bs_incar_dict(incar_dict):
    global wano_file

    dos_dict=wano_file['DOS settings']
    incar_dict['NEDOS'] = dos_dict['Gridpoints']

    if dos_dict['Calculate lm-decomposed DOS']:
        incar_dict['LORBIT'] = 11

    else:
        incar_dict['LORBIT'] = 10

    if not dos_dict['Default energy range']:
        incar_dict['EMIN'] = dos_settings['Lower energy limit [eV]']
        incar_dict['EMAX'] = dos_settings['Upper energy limit [eV]']

    if wano_file['SP calculation'] == 'Band structure calculation':
        if wano_file['Electronic structure settings']['Type of functional'] != 'Hybrid' and wano_file['Band structure settings']['Non-SC calculation']:
            incar_dict['ICHARG'] = 11 

    return incar_dict

def struct_incar_dict(incar_dict):
    global wano_file

    if wano_file['Type of calculation'] == 'MD calculation':
        incar_dict['IBRION'] = 0
        incar_dict['ISYM'] = 0
        md_dict = wano_file['MD settings']

        if md_dict['Ensemble'] == 'NVE':
            incar_dict['MDALGO'] = 1
            incar_dict['ANDERSEN_PROB'] = 0.0

        elif md_dict['Ensemble'] == 'NVT':
            incar_dict['ISIF'] = 2
            mdalgo_dict={'Andersen': 1, 'Nosé-Hoover': 2, 'Langevin': 3, 'Multiple Andersen': 13}
            incar_dict['MDALGO'] = mdalgo_dict[md_dict['Thermostat']]

        incar_dict['NSW'] = md_dict['Number of steps']
        incar_dict['POTIM'] = md_dict['Time step [fs]']
        incar_dict['TEBEG'] = md_dict['Temperature [K]']

        if md_dict['Temperature gradient']:
            incar_dict['TEEND'] = md_dict['Final temperature [K]']

        if md_dict['Save snapshots']:
            incar_dict['NBLOCK'] = md_dict['Snapshot after n steps']

    else:
        if wano_file['Type of calculation'] == 'Structure optimisation':
            opt_dict = wano_file['Structure optimisation']

        elif wano_file['Type of calculation'] == 'NEB calculation':
            opt_dict = wano_file['NEB settings']
            incar_dict['IMAGES'] = opt_dict['Image interpolation']['Number of images']
            incar_dict['SPRING'] = -5
            incar_dict['ISYM'] = 0

        if opt_dict['Optimisation algorithm'] == 'RMM-DIIS':
            incar_dict['IBRION'] = 1
            incar_dict['POTIM'] = opt_dict['Step width scaling constant']

        elif opt_dict['Optimisation algorithm'] == 'Conjugate gradient':
            incar_dict['IBRION'] = 2

        elif opt_dict['Optimisation algorithm'] == 'Damped MD':
            incar_dict['IBRION'] = 3
            incar_dict['SMASS'] = opt_dict['Damping factor'] 
            incar_dict['POTIM'] = opt_dict['Time step [fs]']

        incar_dict['NSW'] = opt_dict['Maximum optimisation steps']

        if not opt_dict['Default convergence criteria']:
            incar_dict['EDIFFG'] = opt_dict['EDIFFG']

        if opt_dict['Optimise cell parameters']:
            incar_dict['ISIF'] = 4 if opt_dict['Constant cell volume'] else 3

        else:
            isif_dict = {'No': 0, 'Trace only': 1, 'Yes': 2}
            incar_dict['ISIF'] = isif_dict[opt_dict['Calculate stress tensor']]

        return incar_dict
    
def write_kpoints(kpoints_dict):
    global incar_dict

    if kpoints_dict['KPOINTS type'] == 'Gamma point only':
        kpoints = Kpoints.gamma_automatic([1,1,1])
        kpoints.write_file('KPOINTS')
        gamma_only = True

    elif kpoints_dict['KPOINTS type'] == 'Set minimum k-spacing':
        k_resol = kpoints_dict['K-resolution [1/Å]']
        incar_dict['KSPACING'] = '%.6f'%round(2*np.pi/k_resol,6)
        incar_dict['KGAMMA'] = bool2vasp(kpoints_dict['Gamma centering'])
        gamma_only = not True in [i/(2*np.pi)>1/k_resol for i in structure.lattice.reciprocal_lattice.abc]

    elif kpoints_dict['KPOINTS type'] == 'Set mesh manually':
        mesh = [int(i) for i in literal_eval(kpoints_dict['Mesh'])[0]]

        if kpoints_dict['Gamma centering']:
            kpoints = Kpoints.gamma_automatic(mesh)

        else: 
            kpoints = Kpoints.monkhorst_automatic(mesh)

        kpoints.write_file('KPOINTS')
        gamma_only = not True in [i > 1 for i in mesh]

    else:
        os.rename('KPOINTS_regular','KPOINTS')
        gamma_only = len(Kpoints.from_file('KPOINTS').kpts[0]) == 1

    return gamma_only

def write_bs_kpoints():

    bs_dict = wano_file['Band structure settings']
    k_path = pymatgen.symmetry.bandstructure.HighSymmKpath(structure).kpath
    path = k_path['path']
    k_points = k_path['kpoints']

    if wano_file['Electronic structure settings']['Type of functional'] != 'Hybrid' and bs_dict['Non-SC calculation']:
        if bs_dict['Non-SC']['Auto-generate k-path']:
            with open('KPOINTS','w') as outfile:
                for i in range(len(path)-1):
                    for j in range(len(path[i])-1): outfile.write(path[i][j]+' - ')
                    outfile.write(path[i][-1]+' | ')
                for i in range(len(path[-1])-1): outfile.write(path[-1][i]+' - ')
                outfile.write(path[-1][-1]+'\n')
                outfile.write('21\nline\nreciprocal\n')
                for i in range(len(path)):
                    for j in range(len(path[i])-1):
                        outfile.write(line(k_points[path[i][j]])+' '+path[i][j]+'\n')
                        outfile.write(line(k_points[path[i][j+1]])+' '+path[i][j+1]+'\n')
                        outfile.write('\n')

        else: os.rename('KPOINT_bs','KPOINTS')

    else:
        if bs_dict['0-weight k-points']['Auto-generate k-points']:
            k_steps = 10

            with open('KPOINTS_0weight','w') as outfile:
                for i in range(len(path)):
                    for j in range(len(path[i])-1):
                        outfile.write(line(k_points[path[i][j]])+'             0 '+path[i][j]+'\n')
                        for n in range(1,k_steps):
                            outfile.write(line(k_points[path[i][j]]+(n/k_steps)*(k_points[path[i][j+1]]-k_points[path[i][j]]))+'             0 \n')
                        outfile.write(line(k_points[path[i][j+1]])+'             0 '+path[i][j+1]+'\n')

        with open('IBZKPT') as infile:
            ibzkpt = infile.readlines()

        nkpt = int(ibzkpt[1].rstrip())

        with open('KPOINTS_0weight') as infile:
            kpt_0w = infile.readlines()

        nkpt_0w = len(kpt_0w)

        with open('KPOINTS','w') as outfile:
            outfile.write(ibzkpt[0])
            outfile.write('    %s\n'%(nkpt+nkpt_0w))
            for i in range(2,nkpt+3):
                outfile.write('%s\n'%ibzkpt[i].rstrip())
            for kpt in kpt_0w:
                outfile.write(kpt)

            if len(ibzkpt)>nkpt+3:
                for i in range(nkpt+3,len(ibzkpt)):
                    outfile.write(ibzkpt[i])

def write_potcar(idx):
    global wano_file,structure

    if 'GGA' in incar_dict or 'METAGGA' in incar_dict:
        potcar_path=os.getenv('VASP_GGA_PATH')

    else:
        potcar_path=os.getenv('VASP_LDA_PATH')

    with open('POSCAR') as infile:
        element_list=infile.readlines()[5].replace('\n','').split()

    extensions = ['','_pv','sv']
    for element in element_list:
        basedir = '%s/%s'%(potcar_path,element)

        while not os.path.isdir('%s%s'%(basedir,extensions[idx])):
            idx=(idx+1)%len(extensions)

        os.system ('cat %s%s/POTCAR >> POTCAR'%(basedir,extensions[idx]))

def run_cmd(incar_dict,cmd,outfile):
    global nkpt

    images = int(incar_dict['IMAGES']) if 'IMAGES' in incar_dict else 1
    cores_per_calc=int(os.getenv('UC_TOTAL_PROCESSORS'))/images
    if round(cores_per_calc%1,5) != 0:
        print('The total number of total cores is no multiple of the number of NEB images. Please adjust your settings and rerun')
        exit(0)
    kpar = math.gcd(nkpt, int(os.getenv('UC_PROCESSORS_PER_NODE')))
    incar_dict['KPAR'] = kpar
    cores_per_kpt=cores_per_calc/kpar
    ncore=math.floor(cores_per_kpt**0.5)
    while not cores_per_kpt%ncore == 0: ncore-=1
    incar_dict['NCORE'] = ncore

    Incar(incar_dict).write_file('INCAR')

    with open(outfile,'w') as outfile:
        process=subprocess.Popen(cmd,stdout=outfile,stderr=subprocess.PIPE)
        out, err = process.communicate()

if __name__ == '__main__':

    outfile='vasp.out'
    nkpt = 1

    with open('rendered_wano.yml') as infile:
        wano_file = yaml.full_load(infile)
    
    neb = wano_file['Type of calculation'] == 'NEB calculation'
    if neb:
        if wano_file['NEB settings']['Image selection'] == 'Load images from tar file':
            os.system('tar -xf images.tar.xz')
            dirnames = []
            for s in os.listdir():
                if os.path.isdir(s):
                    dirnames.append(s)

            dirnames.sort()
            wano_file['NEB settings']['Image interpolation']['Number of images'] = max([int(s) for s in dirnames])-1

        else: 
            n_images = wano_file['NEB settings']['Image interpolation']['Number of images']
            struct_init = Structure.from_file('POSCAR_init')
            atoms = struct_init.species
            sel_dyn = Poscar.from_file('POSCAR_init').selective_dynamics
            pos_init = struct_init.frac_coords
            latt_init = struct_init.lattice.matrix
            struct_final = Structure.from_file('POSCAR_final')
            pos_final = struct_final.frac_coords
            latt_final = struct_final.lattice.matrix
            dirnames = []
            for i in range(n_images+2):
                dirnames.append('%02i'%i)
                os.mkdir(dirnames[i])
                frac_init = 1-(i/(n_images+1))
                frac_final = 1-frac_init
                struct = Structure(species=atoms,lattice=frac_init*latt_init+frac_final*latt_final,coords=frac_init*pos_init+frac_final*pos_final)
                poscar = Poscar(struct,selective_dynamics=sel_dyn)
                poscar.write_file('%s/POSCAR'%dirnames[i])

        shutil.copy('00/POSCAR','.')

    structure=Structure.from_file('POSCAR')
    incar_dict=electronic_incar_dict()
    if wano_file['Type of calculation'] == 'Single point calculation' and not wano_file['SP calculation'] == 'Energy calculation': incar_dict = dos_bs_incar_dict(incar_dict)
    if not (wano_file['Type of calculation'] == 'Single point calculation' and wano_file['SP calculation'] == 'Band structure calculation'): gamma_only=write_kpoints(wano_file['K-sampling'])
    else: 
        write_bs_kpoints()
        gamma_only=False

    if wano_file['Plane waves file']['Auto-generate POTCAR']: 
        paw_type=wano_file['Plane waves file']['PAW type']
        if paw_type == 'Standard': idx=0
        elif paw_type == 'Include semi-core p-states': idx=1
        elif paw_type == 'Include semi-core sp-states': idx=2
        write_potcar(idx)

    if wano_file['Electronic structure settings']['SOC']: suffix='ncl'
    elif gamma_only: suffix='gam'
    else: suffix='std'
    vasp_cmd=('%s %s_%s'%(os.getenv('VASP_PARALLEL'),os.getenv('VASP_BASE'),suffix)).split(' ')


    if suffix != 'gam':
        ibz_dir = 'ibzkpt'
        k_keys = ['KSPACING','KGAMMA']
        minimal_incar = minimal_incar_dict()
        os.mkdir(ibz_dir) 

        if os.path.isfile('KPOINTS'):
            shutil.copy('KPOINTS',ibz_dir)
            shutil.copy('KPOINTS','KPOINTS_0')
        else: 
            with open('KPOINTS_0','w') as kpoints_0:
                for key in k_keys: kpoints_0.write('%s = %s\n'%(key, incar_dict[key]))
            for key in k_keys: minimal_incar[key]=incar_dict[key]
        if not neb: 
            shutil.copy('POSCAR',ibz_dir)

        else:
            shutil.copy('00/POSCAR',ibz_dir)

        os.chdir(ibz_dir)
        Incar(minimal_incar).write_file('INCAR')
        write_potcar(0)
        run_cmd(minimal_incar,vasp_cmd,outfile)
        kpt_file = 'IBZKPT' if os.path.isfile('IBZKPT') else 'KPOINTS'

        with open(kpt_file) as infile:
            nkpt = int(infile.readlines()[1].rstrip())

        for key in k_keys: 
            if key in incar_dict:
                del incar_dict[key]

        os.chdir('..')
        shutil.copy('%s/%s'%(ibz_dir,kpt_file),'KPOINTS')
        shutil.rmtree(ibz_dir)

    if neb:
        os.remove('POSCAR')
        for dirname in ['00',dirnames[-1]]:
            shutil.copy('KPOINTS',dirname)
            shutil.copy('POTCAR',dirname)
            os.chdir(dirname)
            run_cmd(incar_dict,vasp_cmd,outfile)
            os.chdir('..')

    if wano_file['Type of calculation'] != 'Single point calculation':
        incar_dict = struct_incar_dict(incar_dict)
            
    run_cmd(incar_dict,vasp_cmd,outfile)

    output_files = ['POSCAR','POTCAR','INCAR','IBZKPT','KPOINTS','OUTCAR','OSZICAR','XDATCAR','vasp.out','vasprun.xml']

    for filename in output_files:
        if not os.path.isfile(filename): output_files.remove(filename)

    snapshots = []

    if wano_file['Type of calculation'] == 'MD calculation' and wano_file['MD settings']['Save snapshots']:
        poscar = Poscar.from_file('POSCAR')
        sel_dyn = poscar.selective_dynamics
        xdatcar = Xdatcar('XDATCAR')
        structures = xdatcar.structures

        poscar_0 = 'POSCAR_00'
        snapshots.append(poscar_0)
        poscar.write_file(poscar_0)

        for i in range(len(structures)):
            poscar_new = Poscar(structures[i],selective_dynamics=sel_dyn)
            filename = 'POSCAR_%02i'%(i+1)
            snapshots.append(filename)
            poscar_new.write_file(filename)

        output_files += snapshots

    outdict = {'snapshots' : snapshots}
    
    os.system('tar -cf results.tar.xz %s'%(' '.join(output_files)))
    
    with open('output_dict.yml','w') as outfile: yaml.dump(outdict,outfile)

    #for filename in ['CHGCAR','CONTCAR','WAVECAR']:
    #    if not os.path.isfile(filename): open(filename,'w').close()
    if not os.path.isfile('IBZKPT'): shutil.copyfile('KPOINTS','IBZKPT')
