import loompy
import glob
import velocyto as vcy
import numpy as np
from sklearn.manifold import TSNE
import umap
import sys

def main():
    
    if len(sys.argv) != 5:
        print('runs the velocyto workflow on a given hdf5 object')
        print('usage: velocyto_workflow <path to file> <random_seed> <k for kNN smoothing> <tsne/umap>')
        sys.exit()
        
    input_path = sys.argv[1]
    random_seed = sys.argv[2]
    k = sys.argv[3]
    viz = sys.argv[4]
    
    if viz == 'tsne':
        output_path = input_path[:-5] + '_kNN' + k + '_tsne' + random_seed + '.hdf5'
    else:
        output_path = input_path[:-5] + '_kNN' + k + '_umap.hdf5'
    
    print('loading data')
    vlm = vcy.load_velocyto_hdf5(input_path)

    print(len(vlm.ca['CellID']), 'cells')
    print(len(vlm.ra['Gene']), 'genes')

    print('filtering cells')
    vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.5))

    print('filtering genes')
    vlm.score_cv_vs_mean(3000, plot=False, max_expr_avg=35)
    vlm.filter_genes(by_cv_vs_mean=True)

    print(len(vlm.ca['CellID']), 'cells')
    print(len(vlm.ra['Gene']), 'genes')

    #print('setting sample names as clusters')
    #samplenames = list(map(lambda x: x.split(':')[0], vlm.ca['CellID']))
    #vlm.ca['sample_name'] = samplenames
    #vlm.set_clusters(vlm.ca["sample_name"])

    print('normalizing data matrices')
    vlm._normalize_S(relative_size=vlm.S.sum(0), target_size=vlm.S.sum(0).mean())
    vlm._normalize_U(relative_size=vlm.U.sum(0), target_size=vlm.U.sum(0).mean())

    print('running pca')
    vlm.perform_PCA()

    print('knn smoothing')
    #default k = 0.025 * nCells
    vlm.knn_imputation(n_pca_dims=15, k=int(k), balanced=True, b_sight=3000, b_maxl=1500, n_jobs=20)

    print('fit gammas')
    vlm.fit_gammas()

    print('calculate velocity')
    vlm.predict_U()
    vlm.calculate_velocity()
    vlm.calculate_shift(assumption="constant_velocity")
    vlm.extrapolate_cell_at_t(delta_t=1.)

    if viz == 'tsne':
        #print('running tsne')
        #bh_tsne = TSNE(random_state=int(random_seed))
        #vlm.ts = bh_tsne.fit_transform(vlm.pcs[:, :15])
        vlm.ts = np.column_stack([vlm.ca["tsne1"], vlm.ca["tsne2"]])
        print('projection of velocity onto embeddings')
        vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1, n_neighbors=3500, knn_random=True, sampled_fraction=0.5)
    else:
        print('running UMAP')
        vlm.umap = umap.UMAP().fit_transform(vlm.pcs[:, :18])
        print('projection of velocity onto embeddings')
        vlm.estimate_transition_prob(hidim="Sx_sz", embed="umap", transform="sqrt", psc=1, n_neighbors=3500, knn_random=True, sampled_fraction=0.5)
    
    print('calculate embedding shift')
    vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)

    print('calculate grid arrows')
    vlm.calculate_grid_arrows(smooth=0.5, steps=(40, 40), n_neighbors=150)

    print('saving hdf5')
    vlm.to_hdf5(output_path)


if __name__ == '__main__':
    main()
