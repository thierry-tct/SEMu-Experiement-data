
import os
import sys
import shutil
import json
from multiprocessing.pool import ThreadPool
from multiprocessing import Pool 

sys.path.insert(0, 'D:\\PhD\\muteria')
try:
    import muteria.common.fs as common_fs
    import muteria.common.matrices as common_matrices
    import muteria.statistics.mutant_quality_indicators as mqi
except ImportError:
    print("# ERROR: Muteria import failed")
    exit(1)


def error_exit(s):
    print("Error: "+s)
    exit(1)
#~ def error_exit()

def merge_matrices(in_init_mat, in_mfi_mat, in_killnon_mat, out_final_mat_file=None, mart_0=True):
    """
        Input 3 matrices and create output matrice and return it
    """
    non_key_cols = set()
    non_key_cols |= set(in_init_mat.get_nonkey_colname_list())
    non_key_cols |= set(in_mfi_mat.get_nonkey_colname_list())
    non_key_cols |= set(in_killnon_mat.get_nonkey_colname_list())

    non_key_cols = list(non_key_cols)
    
    # 1. merge matrices into final
    print("        >> step 1..")
    final_mat_raw = common_matrices.ExecutionMatrix(non_key_col_list=non_key_cols)
    final_mat_raw.update_with_other_matrix(in_init_mat, allow_missing=True)
    print("            - step 1 - updating 2..")
    final_mat_raw.update_with_other_matrix(in_mfi_mat, allow_missing=True)
    print("            - step 1 - updating 3..")
    final_mat_raw.update_with_other_matrix(in_killnon_mat, allow_missing=True)
    ## change uncertain into inactive
    uncertain = final_mat_raw.query_uncertain_columns_of_rows()
    for row, col_list in uncertain.items():
        if len(col_list) == 0:
            continue
        new_v = {c: final_mat_raw.getInactiveCellVal() for c in col_list}
        final_mat_raw.update_cells(row, new_v)
    
    # 2. transform the columns to have mart_0 if set
    print("        >> step 2..")
    final_mat = common_matrices.ExecutionMatrix(filename=out_final_mat_file, non_key_col_list=non_key_cols)
    for k, v in final_mat_raw._get_key_values_dict().items(): 
        if mart_0:
            k = ':'.join(['mart_0', k.split('/')[1]])
        final_mat.add_row_by_key(k, v)

    return final_mat
#~ def merge_matrices()
    
def get_mutants_info(raw_mut_inf, fdupes):
    """ take in_raw_mutant_info and fdupes
        return final mutant info
    """
    res = {}
    dups = set()
    for rem, duplist in fdupes.items():
        dups |= set(duplist)
    remaining = set(raw_mut_inf) - dups
    for m_id in remaining:
        res[m_id] = raw_mut_inf[m_id]
    return res
#~ def get_mutants_info()

def extract_proj(proj_in_topdir, proj_out_topdir):
    out_SM_mat_file = os.path.join(proj_out_topdir, "STRONG_MUTATION.csv")
    out_mut_info_file = os.path.join(proj_out_topdir, "mutantsInfos.json")
    out_raw_mut_info_file = os.path.join(proj_out_topdir, "raw-mutantsInfos.json")
    out_subs_cluster_file = os.path.join(proj_out_topdir, "subsuming-clusters.json")
    #out_test_equivalent_file =  os.path.join(proj_out_topdir, "test_equivalent.json")

    in_raw_mut_info_file = os.path.join(proj_in_topdir, 'inputs', 'mutantsdata', 'mutantsInfos.json')
    in_fduped_file = os.path.join(proj_in_topdir, 'inputs', 'mutantsdata', 'fdupes_duplicates.json')
    in_initial_sm_mat_file = os.path.join(proj_in_topdir, 'inputs', 'matrices', 'SM.dat')
    in_mfi_run_sm_mat_file = os.path.join(proj_in_topdir, 'OUTPUT', 'mfirun_output', 'data', 'matrices', 'matrices', 'SM.dat')
    in_killed_non_mat_file = os.path.join(proj_in_topdir, 'OUTPUT', 'killed_non_mfirun_output', 'data', 'matrices', 'matrices', 'SM.dat')

    # Computation
    # 1. mutantsinfo
    print("    # MutantsInfo ...")
    in_raw_mut_inf = common_fs.loadJSON(in_raw_mut_info_file)
    in_fdupes = common_fs.loadJSON(in_fduped_file)
    mut_info = get_mutants_info(in_raw_mut_inf, in_fdupes)
    common_fs.dumpJSON(mut_info, out_mut_info_file, pretty=True)
    shutil.copy2(in_raw_mut_info_file, out_raw_mut_info_file)

    # 2. matrices
    print("    # Matrices ...")
    tmp_file = os.path.join(proj_out_topdir, 'tmp_file')
    load_mat = [None, None, None]
    for ind, mat_file in enumerate([in_initial_sm_mat_file, in_mfi_run_sm_mat_file, in_killed_non_mat_file]):
        with open(mat_file) as f:
            f_dat = f.read().replace('ktestSM', 'MUTERIA_MATRIX_KEY_COL')
        assert "ktestSM" not in f_dat, "after replace, ktestSM in matrice {}".format(matfile)
        with open(tmp_file, "w") as g:
            g.write(f_dat)
        load_mat[ind] = common_matrices.ExecutionMatrix(filename=tmp_file)
        os.remove(tmp_file)
    in_init_mat, in_mfi_mat, in_killnon_mat = load_mat
    
    if os.path.isfile(out_SM_mat_file):
        os.remove(out_SM_mat_file)
    final_mat = merge_matrices(in_init_mat, in_mfi_mat, in_killnon_mat, out_final_mat_file=out_SM_mat_file)
    final_mat.serialize()

    # 3. subsuming
    print("    # Subsuming ...")
    eq, subs_clust = mqi.getSubsumingMutants (out_SM_mat_file)
    common_fs.dumpJSON({'subsume': [None, subs_clust]}, out_subs_cluster_file, pretty=True)
    
#~ def extract_proj

def parallel_extract_wrapper (arg_tuple):
    p, p_in, p_out = arg_tuple

    print ("# Processing {} ...".format(p))
        
    if not os.path.isdir(p_out):
        os.mkdir(p_out)

    if not os.path.isfile(os.path.join(p_out, "subsuming-clusters.json")):
#       try:
        extract_proj(p_in, p_out)
#       except:
#           error_exit ("\n\n >>>>> {} FAILED".format(p))
    print ("# Done with {}!".format(p))
#~ def parallel_extract_wrapper()
    
def main():
    src_top_dir=os.path.abspath(os.path.dirname(sys.argv[0]))
    dest_dir=os.path.join(os.path.dirname(os.path.dirname(src_top_dir)), "semu_cleaned_data")
    
    if not os.path.isdir(dest_dir):
        os.mkdir(dest_dir)

    tasks = []
    for p in os.listdir(src_top_dir):
        p_in = os.path.join(src_top_dir, p)
        p_out = os.path.join(dest_dir, p)
        if not os.path.isdir(p_in): # or p != 'false':
            continue
        
        tasks.append((p, p_in, p_out))

    thread = False
    if thread:
        pool = ThreadPool(24)
        pool.map(parallel_extract_wrapper, tasks)
        pool.terminate()
        pool.close()
        pool.join()
    else:
        with Pool(6) as pool:
            pool.map(parallel_extract_wrapper, tasks)

    print("# Data Conversion Done!")

    print("# Getting source files ...")
    #repo_tmp = os.path.join(os.path.dirname(dest_dir), "coreutils_repo.tmp")
    repo_tmp = os.path.join(os.path.dirname(src_top_dir), "coreutils-v8.22")
    #git_cmd = "git clone --branch v8.22 --single-branch {} {}".format('https://github.com/coreutils/coreutils', repo_tmp)
    #if os.system(git_cmd) != 0:
    #    error_exit("git clone failed.\nCMD: {}".format(git_cmd))
    for p, p_in, p_out in tasks:
        minf_file = os.path.join(p_out, "mutantsInfos.json")
        srcs = set()
        with open(minf_file) as f:
            minf = json.load(f)
        for mut, m_dat in minf.items():
            if m_dat["SrcLoc"] == "":
                continue
            srcs.add(os.sep.join(m_dat["SrcLoc"].split(':')[0].split('/')))
                
        for src in srcs:
            dest = os.path.join(p_out, os.path.dirname(src))
            if not os.path.isdir(dest):
                os.makedirs(dest)
            dest_src = os.path.join(dest, os.path.basename(src))
            shutil.copy2(os.path.join(repo_tmp, src), dest_src)
            if os.path.basename(dest_src) == "false.c":
                with open(dest_src) as f:
                    f_dat = f.readlines()
                with open(os.path.join(repo_tmp, src.replace('false', 'true'))) as f:
                    true_dat = f.read()
                with open(dest_src, 'w') as f:
                    for line in f_dat:
                        if '#include "true.c"' in line:
                            line = line.replace('#include "true.c"', '/*#include "true.c"*/')
                        f.write(line)
                    f.write(true_dat)
                
    #try:
    #    shutil.rmtree(repo_tmp)
    #except PermissionError:
    #    error_exit("ALL is DONE. Just manually remove the cloned coreutils folder: {}".format(repo_tmp))
    print("# ALL DONE!")
#~ def main()


if __name__ == "__main__":
    main()
