import numpy as np
import matplotlib.pyplot as plt
from .dropcalling import countMatrix

def bc_rk_plt(raw_mtx:countMatrix, init_idx, eval_idx, fig_path:str='.'):

    my_bc_colors = {'blue_1': '#08336E',
                    'blue_2': '#105C9C',
                    'blue_3': '#3888C0',
                    'blue_3': '#68ACCD',
                    'blue_4': '#AAD7E5',
                    'blue_5': '#D2E3F3',
                    'blue_6': '#F4F9FE'}

    bc_colors = np.array([my_bc_colors['blue_5'] for _ in range(len(raw_mtx.bcs))])
    bc_colors[init_idx] = my_bc_colors['blue_1']
    bc_colors[eval_idx] = my_bc_colors['blue_2']
    ind_plt = np.argsort(raw_mtx.counts_per_bc)[::-1]                       

    fig = plt.figure()
    fig1 = fig.add_subplot(111)

    fig1.scatter(x=range(1, len(raw_mtx.counts_per_bc)+1), y=raw_mtx.counts_per_bc[ind_plt],
                 c=bc_colors[ind_plt], s=7, linewidths=0, zorder=10)
    fig1.loglog()
    fig1.set_xlim(0, len(raw_mtx.counts_per_bc)*1.25)
    fig1.set_xlabel('Barcode index')
    fig1.set_ylabel('Count')
    fig1.spines['bottom'].set_linewidth(0.5)
    fig1.spines['left'].set_linewidth(0.5)
    fig1.spines['top'].set_visible(False)
    fig1.spines['right'].set_visible(False)
    fig1.set_xticks([], minor=True)
    fig1.set_yticks([], minor=True)
    fig1.tick_params(bottom=False, top=False, left=False, right=False)
    fig1.grid(color="#F0F0F2", zorder=0)

    fig.savefig(f"{fig_path}/bcs_rk_plt.png", bbox_inches='tight', dpi=400, transparent=True)


def calcu_part_cellp(bc_count, init_idx, eval_idx, scale:int=100):
    # bc_count 未降序排序
    called_cell_idx = np.append(init_idx, eval_idx)
    bc_count_dec_idx = np.argsort(bc_count)[::-1]
    # baseline_idx_tmp1 = np.where(np.in1d(bc_count_dec_idx, called_cell_idx))
    baseline_idx = max(np.where(np.in1d(bc_count_dec_idx, called_cell_idx))[0])
    # baseline_idx = max(called_cell_idx)
    cellp = list()
    px = list()
    py = list()
    sidx = 0
    eidx = scale - 1
    while eidx <= baseline_idx:
        px.append(range(sidx+1, eidx+1))
        py.append(bc_count[bc_count_dec_idx[sidx:eidx]])
        cellp.append(sum(np.in1d(bc_count_dec_idx[sidx:eidx+1], called_cell_idx)))
        sidx += scale
        eidx += scale        
    ept_drops_sidx = eidx+1
    return px, py, np.array(cellp)/scale, ept_drops_sidx


def bc_curve_plt(raw_mtx:countMatrix, init_idx, eval_idx, fig_path:str='.', scale:int=100):
    
    px, py, cellp, ept_drops_sidx = calcu_part_cellp(bc_count=raw_mtx.counts_per_bc, init_idx=init_idx, eval_idx=eval_idx, scale=scale)
    ind_plt = np.argsort(raw_mtx.counts_per_bc)[::-1]                       

    fig = plt.figure()
    fig1 = fig.add_subplot(111)
    my_color = [102/255,8/255,116/255]                          # Tsinghua purple

    for x, y, c in zip(px, py, cellp):
        pcolor = my_color
        fig1.plot(x, y, color=pcolor, alpha=c*c, linewidth=2)    
    fig1.plot(range(ept_drops_sidx+1, raw_mtx.bcs_dim+1), raw_mtx.counts_per_bc[ind_plt[ept_drops_sidx:]], color=my_color, alpha=0.1, linewidth=2.2)
    
    fig1.loglog()
    fig1.set_xlim(0, len(raw_mtx.counts_per_bc)*1.25)
    fig1.set_xlabel('Barcode index')
    fig1.set_ylabel('Count')
    fig1.spines['bottom'].set_linewidth(0.5)
    fig1.spines['left'].set_linewidth(0.5)
    fig1.spines['top'].set_visible(False)
    fig1.spines['right'].set_visible(False)
    fig1.set_xticks([], minor=True)
    fig1.set_yticks([], minor=True)
    fig1.tick_params(bottom=False, top=False, left=False, right=False)
    fig1.grid(color="#F0F0F2", zorder=0)

    fig.savefig(f"{fig_path}/bcs_rk_plt.png", bbox_inches='tight', dpi=400, transparent=True)
