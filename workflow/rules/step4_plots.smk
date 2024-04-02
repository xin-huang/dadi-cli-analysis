rule all:
    input:
        "results/plots/1KG.two.epoch.gamma.dfes.png",
        "results/plots/1KG.two.epoch.lognormal.dfes.png",


rule plot:
    input:
        res = "results/dfes/all.res.txt",
    output:
        gamma_dfe = "results/plots/1KG.two.epoch.gamma.dfes.png",
        lognormal_dfe = "results/plots/1KG.two.epoch.lognormal.dfes.png"
    run:
        import matplotlib.pyplot as plt
        import pandas as pd
        import matplotlib
        matplotlib.use("Agg")

        data = pd.read_csv(input.res, sep="\t")
        population_list = {
            'AFR': ['ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI'],
            'AMR': ['CLM', 'MXL', 'PEL', 'PUR'],
            'EAS': ['CDX', 'CHB', 'CHS', 'JPT', 'KHV'],
            'EUR': ['CEU', 'FIN', 'GBR', 'IBS', 'TSI'],
            'SAS': ['BEB', 'GIH', 'ITU', 'PJL', 'STU'],
        }
        sorter = [
            'ACB', 'ASW', 'ESN', 'GWD', 'LWK', 'MSL', 'YRI',
            'CLM', 'MXL', 'PEL', 'PUR',
            'CDX', 'CHB', 'CHS', 'JPT', 'KHV',
            'CEU', 'FIN', 'GBR', 'IBS', 'TSI',
            'BEB', 'GIH', 'ITU', 'PJL', 'STU',
        ]
        colors = [
            'black', 'black', 'black', 'black', 'black', 'black', 'black',
            'green', 'green', 'green', 'green',
            'gold', 'gold', 'gold', 'gold', 'gold',
            'blue', 'blue', 'blue', 'blue', 'blue',
            'brown', 'brown', 'brown', 'brown', 'brown',
        ]
        df = data[data['DFE'] == 'gamma']
        gamma_df = df.sort_values(by="Pop", key=lambda column: column.map(lambda e: sorter.index(e)))
        df = data[data['DFE'] == 'lognormal']
        lognormal_df = df.sort_values(by="Pop", key=lambda column: column.map(lambda e: sorter.index(e)))

        fig, axs = plt.subplots(nrows=2, ncols=2, constrained_layout=True, figsize=(7.5,4), dpi=350)
        gridspec = axs[0, 0].get_subplotspec().get_gridspec()
        for a in axs[:,1]: a.remove()
        axs[0, 0].scatter(sorter, gamma_df['DFE_param1'].values, color=colors, zorder=2)
        axs[0, 0].scatter(sorter, gamma_df['DFE_param1_lb'].values, marker='_', color='grey')
        axs[0, 0].scatter(sorter, gamma_df['DFE_param1_ub'].values, marker='_', color='grey')
        for i in list(range(0,26)): 
            axs[0, 0].plot([i, i], [gamma_df['DFE_param1'].values[i], gamma_df['DFE_param1_ub'].values[i]], linestyle='dashed', color='grey', zorder=1)
            axs[0, 0].plot([i, i], [gamma_df['DFE_param1_lb'].values[i], gamma_df['DFE_param1'].values[i]], linestyle='dashed', color='grey', zorder=1)
        axs[0, 0].set_ylim([0,0.5])
        axs[0, 0].set_ylabel('shape')
        axs[0, 0].set_xticks(list(range(0,26)), sorter, rotation=90)
        axs[1, 0].scatter(sorter, gamma_df['DFE_param2'].values, color=colors, zorder=2)
        axs[1, 0].scatter(sorter, gamma_df['DFE_param2_lb'].values, marker='_', color='grey')
        axs[1, 0].scatter(sorter, gamma_df['DFE_param2_ub'].values, marker='_', color='grey')
        for i in list(range(0,26)):
            axs[1, 0].plot([i, i], [gamma_df['DFE_param2'].values[i], gamma_df['DFE_param2_ub'].values[i]], linestyle='dashed', color='grey', zorder=1)
            axs[1, 0].plot([i, i], [gamma_df['DFE_param2_lb'].values[i], gamma_df['DFE_param2'].values[i]], linestyle='dashed', color='grey', zorder=1)
        axs[1, 0].set_xticks(list(range(0,26)), sorter, rotation=90)
        axs[1, 0].set_yticks([0, 100000, 200000, 300000], ['0', '1', '2', '3'])
        axs[1, 0].set_ylim([-10000, 300000])
        axs[1, 0].set_ylabel('scale (× $10^5$)')

        subfig = fig.add_subfigure(gridspec[:,1])
        handles, labels = subfig.gca().get_legend_handles_labels()
        afr = axs[0,1].scatter([0], [0], label='AFR', color='black')
        amr = axs[0,1].scatter([0], [0], label='AMR', color='green')
        eas = axs[0,1].scatter([0], [0], label='EAS', color='gold')
        eur = axs[0,1].scatter([0], [0], label='EUR', color='blue')
        sas = axs[0,1].scatter([0], [0], label='SAS', color='brown')
        handles.extend([afr, amr, eas, eur, sas])
        subfig.legend(handles=handles, fontsize=8, handlelength=1.5, loc='upper left')
        fig.set_constrained_layout_pads(w_pad=4 / 72, h_pad=4 / 72, hspace=0, wspace=0.1)
        plt.savefig(output.gamma_dfe, bbox_inches='tight')

        fig, axs = plt.subplots(nrows=2, ncols=2, constrained_layout=True, figsize=(7.5,4), dpi=350)
        gridspec = axs[0, 0].get_subplotspec().get_gridspec()
        for a in axs[:,1]: a.remove()
        axs[0, 0].scatter(sorter, lognormal_df['DFE_param1'].values, color=colors, zorder=2)
        axs[0, 0].scatter(sorter, lognormal_df['DFE_param1_lb'].values, marker='_', color='grey')
        axs[0, 0].scatter(sorter, lognormal_df['DFE_param1_ub'].values, marker='_', color='grey')
        for i in list(range(0,26)): 
            axs[0, 0].plot([i, i], [lognormal_df['DFE_param1'].values[i], lognormal_df['DFE_param1_ub'].values[i]], linestyle='dashed', color='grey', zorder=1)
            axs[0, 0].plot([i, i], [lognormal_df['DFE_param1_lb'].values[i], lognormal_df['DFE_param1'].values[i]], linestyle='dashed', color='grey', zorder=1)
        axs[0, 0].set_ylim([-1,7])
        axs[0, 0].set_ylabel('$\mu$')
        axs[0, 0].set_xticks(list(range(0,26)), sorter, rotation=90)
        axs[1, 0].scatter(sorter, lognormal_df['DFE_param2'].values, color=colors, zorder=2)
        axs[1, 0].scatter(sorter, lognormal_df['DFE_param2_lb'].values, marker='_', color='grey')
        axs[1, 0].scatter(sorter, lognormal_df['DFE_param2_ub'].values, marker='_', color='grey')
        for i in list(range(0,26)):
            axs[1, 0].plot([i, i], [lognormal_df['DFE_param2'].values[i], lognormal_df['DFE_param2_ub'].values[i]], linestyle='dashed', color='grey', zorder=1)
            axs[1, 0].plot([i, i], [lognormal_df['DFE_param2_lb'].values[i], lognormal_df['DFE_param2'].values[i]], linestyle='dashed', color='grey', zorder=1)
        axs[1, 0].set_xticks(list(range(0,26)), sorter, rotation=90)
        axs[1, 0].set_ylim([0, 40])
        axs[1, 0].set_ylabel('$\sigma$')

        subfig = fig.add_subfigure(gridspec[:,1])
        handles, labels = subfig.gca().get_legend_handles_labels()
        afr = axs[0,1].scatter([0], [0], label='AFR', color='black')
        amr = axs[0,1].scatter([0], [0], label='AMR', color='green')
        eas = axs[0,1].scatter([0], [0], label='EAS', color='gold')
        eur = axs[0,1].scatter([0], [0], label='EUR', color='blue')
        sas = axs[0,1].scatter([0], [0], label='SAS', color='brown')
        handles.extend([afr, amr, eas, eur, sas])
        subfig.legend(handles=handles, fontsize=8, handlelength=1.5, loc='upper left')
        fig.set_constrained_layout_pads(w_pad=4 / 72, h_pad=4 / 72, hspace=0, wspace=0.1)
        plt.savefig(output.lognormal_dfe, bbox_inches='tight')
