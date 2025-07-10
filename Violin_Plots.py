#!/usr/bin/env python3

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import rdMolTransforms, rdmolops, rdMolDescriptors, Descriptors
from rdkit.Chem import AllChem 
import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler, MinMaxScaler
import statistics


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colorbar
import matplotlib.colors
import matplotlib.cm
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.gridspec as gridspec

import Compile_Data

plt.rcParams['figure.figsize'] = [12, 8]

plt.rcParams["font.fantasy"] = "Cambria"



class graph_ensambles:
    def __init__(self, df=None):
        df_j = Compile_Data.curate_data(exclude_irrelavent_features=True, boltz_99=True,  remove_high_energy_confs=False, drop_low_data=False, export_csv=False, outfile=None)
        if df is None:
            self.df = df_j.generate_dataset()
        else:
            self.df = df

        self.weighted_dataset = df_j.boltz_avg_data(self.df, 'G(MeCN)')
        self.bridge_ls = df_j.get_bridge_ls(self.df)
        

        self.scaled_df = df_j.scale_data(self.df, 0, 1)
        properties = []
        numeric_list = ['int', 'float', 'number'] 
        # select columns based on the above list
        numeric_data = self.df.select_dtypes(include=numeric_list)
        
        for key in numeric_data.keys():
            properties.append(key)
        
        self.properties = properties
        
    def get_X_vals(self, data):
        np.random.seed(42)
        x_vals = np.random.normal(0, 0.5, size=len(data))
        adder = 0
        for i in range(len(x_vals)):
            x_vals[i] += adder
            if i % 10 == 0:
                adder += .5
        
        return x_vals
    
    def get_BB_data(self):
        """ Creating phosphine linker tags for plotting """
        
        sort_by_part_dict = {
            'C(SP3)' : 'a',
            'N(SP3)' : 'b',
            'O(SP3)' : 'c',
            'Fe(SP3D)' : 'd',
            'C(SP2)' : 'e',
            'N(SP2)' : 'f',
            'O(SP2)' : 'g',
            'C(SP)' : 'h'
        }
        
        sp3_l = ['a', 'b', 'c']
        sp2_l = ['e', 'f', 'g']
        
        
        
        for i, x in enumerate(self.bridge_ls):
            if 'ferrocene005' in str(x[0]):
                self.bridge_ls[i][1] = '2C(SP2),1Fe(SP3D)'
        
        
        ligs = list(set(list(self.df['Ligand'].values)))
        sorted_bridges = []
        
        
        for i, x in enumerate(self.bridge_ls):
            total_score=0
            if x[0] in ligs:
                try:
                    xs = x[1].split(',')
                    n_tags=[]
                    tags=''
                    for s in xs:
                        tag = sort_by_part_dict[s[1:]]
                        tags += tag
                        tag_n = f'{s[0]}{tag}'
                        n_tags.append([s[0], tag])
                        num = s[0]
                        let = tag
                        if let in sp3_l:
                            total_score += int(num)*1
                        elif let in sp2_l:
                            total_score += int(num)*-1
                        elif let == 'd':
                            total_score += int(num)*-10
                    
                    
                    
                    n_tags.sort(key = lambda row: row[1])
                    n_tags = [f'{x[0]}{x[1]}' for x in n_tags]
                    n_tags= ','.join(n_tags)
                    
                    sorted_bridges.append([x[0], n_tags, total_score])
                except:
                    num = x[1][0]
                    tag = sort_by_part_dict[x[1][1:]]
                    let = tag
                    if let in sp3_l:
                        total_score += int(num)*1
                    elif let in sp2_l:
                        total_score += int(num)*-1
                    elif let == 'd':
                        total_score += int(num)*-10
                    
                    
                    
                    sorted_bridges.append([x[0], f'{num}{tag}', total_score])
                
        
        
        sorted_bridges.sort(key = lambda row: row[2])
        ligands = [x[0] for x in sorted_bridges]
        
        
        tt2_ls = [(x[1], x[2]) for x in sorted_bridges]
        tt2_ls = list(set(tt2_ls))
        tt2_ls.sort(key = lambda row: row[1])
        bb_type = [x[0] for x in tt2_ls]
        
        
        tt1_ls = [(x[0], bb_type.index(x[1])) for x in sorted_bridges]
        tt1_ls.sort(key = lambda row: row[1])
        ligands = [x[0] for x in tt1_ls]
        
        self.bridge_ls = sorted_bridges
        
        
        self.bb_colors = plt.cm.viridis(np.linspace(0, 1, len(bb_type)))
        self.norm = mpl.colors.Normalize(vmin=0, vmax=len(bb_type))
        self.cmap = LinearSegmentedColormap.from_list('viridis', self.bb_colors)
        self.cbar_data = np.random.rand(len(bb_type), len(bb_type))
        
        ref_cmap = {}
        for bb in bb_type:
            c = self.bb_colors[bb_type.index(bb)]
            #c = np.array([c])
            ref_cmap[bb] = c
        
        self.bb_type = bb_type
        
        self.ref_cmap = ref_cmap
        
        
        #no_ferro = [x for x in ligands if 'errocene' not in x]
        #fer = [x for x in ligands if 'errocene' in x]
        #self.ligands = no_ferro + fer
        self.ligands = ligands
        self.x_vals = self.get_X_vals(ligands)
        self.x_vals.sort()

    def set_vplot_BB_colors(self, vplot, ligand):
        for cb in vplot['bodies']: 
            ss = [x[1] for x in self.bridge_ls if x[0] == ligand]
            c = self.bb_colors[bb_type.index(ss[0])]
            c = np.array([c])
            cb.set_facecolor(c)
            cb.set_edgecolor('black')
            cb.set_linewidth(1)
            cb.set_alpha(0.8)
        
        return vplot





    def jitter(self, datasets, colors, energies, sorted_NRG, ax, nx, point_scale=1):
        """Scatter points that may overlap when graphing by randomly offsetting them."""
        for i, p in enumerate(datasets):
            NRG = energies[i]
            idx = sorted_NRG.index(NRG)
            if idx == 0:
                marker = "D"
                color = 'black'
                edgecol = 'gold'
                markeredgewidth = 3
                markersize = 15
            else:
                marker = 'o'
                edgecol = 'black'
                color = colors[idx]
                markeredgewidth = 1.5
                markersize = 10
                
            y = [p]
            x = np.random.normal(nx, 0.020, size=len(y))
            ax.plot(x, y, alpha=1, markersize=round(markersize*point_scale), color=color, marker=marker, markeredgecolor=edgecol,
                    markeredgewidth=markeredgewidth, linestyle='None')


    def get_cb_cmap(self, colors, values):
        
        
        scaled_values = values
        clist = []
        for c, v in zip(colors, scaled_values):
            
            clist.append((v, c))
        
        cm = mpl.colors.LinearSegmentedColormap.from_list('val_cmap', clist, 1024)
        
        return cm 
    
    def vplot_single(self, outfig, properties, ligands):
        
        
        df = drop_low_data(self.df)
        scaled_df = scale_data(df, 0, 1)
        #ligands = list(set(list(df['Ligand'].values)))
        #properties = [x for x in properties if x not in excluded_features]
        weighted_dataset = boltz_avg_data(df)
        x_vals = get_X_vals(properties)
        for ligand in ligands:
            energies = list(df['def2-TZVP-GAS'].loc[df['Ligand'] == ligand].astype(float))
            
            for i, prop in enumerate(properties):
                data = list(df[prop].loc[df['Ligand'] == ligand].astype(float))
                scaled_data = list(scaled_df[prop].loc[df['Ligand'] == ligand].astype(float))
                data = [val for val in data if val != 0]
                scaled_data = [val for val in scaled_data if val != 0]
                if len(data) < 2 or len(scaled_data) == 0:
                    continue
                    
                
                fig,ax = plt.subplots(figsize=(6,10), dpi=self.dpi)
                norm = mpl.colors.Normalize(vmin=0, vmax=len(energies))
                
                
                x = [x_vals[i]]
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.set_facecolor('lavender')
                
                
                vplot = ax.violinplot(data,x,showmeans=False, showmedians=False,
                showextrema=False)
                for cb in vplot['bodies']: 
                    cb.set_facecolor('white')
                    
                    cb.set_edgecolor('darkslategray')
                    cb.set_alpha(1.0)
                    cb.set_linewidth(2)
                    
                boltz_avg = list(weighted_dataset['Boltz_avg'].loc[(weighted_dataset['Ligand'] == ligand) & (weighted_dataset['Property'] == prop)].astype(float))
                boltz_avg = boltz_avg[0]
                    
                colors = plt.cm.Blues(np.linspace(0, 1, len(energies)))
                norm = mpl.colors.Normalize(vmin=0, vmax=len(energies))
                cmap = LinearSegmentedColormap.from_list('Blues', colors)
                
                NRGs = energies
                NRGs.sort()
                 
                plt.set_cmap(cmap)
                
                jitter(data, colors, energies, NRGs, ax, x[0])
                
                L_lim = abs(min(data))-abs((abs(max(data)) - abs(min(data)))*0.10)
                H_lim = abs(max(data))+abs((abs(max(data)) - abs(min(data)))*0.20)
                
                
                
                scaled_rng = max(scaled_data) - min(scaled_data)
                
                
                pro = prop.replace('/', '-')
                fig_name = f'{fig_dir}/{pro}_{ligand}.png'
                
                
                ax.set_ylim(L_lim, H_lim)
                plt.axhline(boltz_avg, 0.25, 0.75, color="black", linestyle="--", lw=3)
                plt.axhline(boltz_avg, 0.40, 0.60, color="black", linestyle="-", lw=5)
                ax.plot(statistics.mean(x), boltz_avg, marker = 's', markeredgecolor='black', color = 'white', markersize = 15, markeredgewidth = 3)
                
                ## SFR = Scaled Feature Range
                fig.text(0.50, 0.92, f'SFR: {round(scaled_rng, 2)}', fontsize=18, fontweight='bold', ha = 'center', va = 'top', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=1'))
                ax.spines['right'].set_visible(False)
                plt.xticks([])
                plt.yticks([])
                plt.grid(axis = 'y', color = 'white', linestyle = '--', linewidth = 1.0)
                plt.ylim(L_lim, H_lim)
                plt.tight_layout()
                plt.savefig(f'{outfig}.png', transparant='True')
                plt.show()
                
                
                plt.close()



    def get_feat_range_plot(self, prop_ranges, outfile, dpi):
        
        fig,ax = plt.subplots(figsize=(25,8), dpi=self.dpi)
        for i, pr in enumerate(prop_ranges):
            plt.bar(i, pr[1], edgecolor = 'black', linewidth=1, color = pr[2], width = 0.6, label=pr[3])
        
        data = [pr[1] for pr in prop_ranges]
        plt.ylabel('Average feature range', fontsize=24, labelpad=10)
        plt.xlabel('Feature', fontsize=24, labelpad=10)
        plt.xticks([])
        plt.yticks([])
        ax.tick_params(axis='y', labelsize=20)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        y_max = round(max(data), 2) + 0.01
        y_half = round((max(data)/2), 2)
        plt.ylim(0, round((max(data)+(max(data)*0.15)), 2))
        
        plt.legend(['Electronic', 'Steric'], loc='upper left', fontsize=22, facecolor='lavender', edgecolor='black', framealpha=0.5)
        leg = ax.get_legend()
        leg.legendHandles[0].set_color('steelblue')
        leg.legendHandles[1].set_color('indianred')
        
        fign = fn.replace('.csv', '.png')
        plt.savefig(f'{outfile}.png')
        plt.show()


    def get_label_color(self, feat):
        steric_feats = ['Bmax', 'Bmin', 'L', 'Vshell', 'distance', 'Vbur', 'Area', 'Volume', '6-5-7', '6-7', '5-7', '5-6', 'anisotropy', 'Pd-P difference']
        steric=False
        for f in steric_feats:
            if f in feat or f == feat:
                steric=True
        if steric:
            c = 'indianred'
            
            l = 'Steric'
        else:
            c = 'steelblue'
            
            l = 'Electronic'
        
        return c, l


    def get_feature_ranges(self, outfile, ferro=False, noferro=False):
        
        
        ## List of features excluded from this analyses due to: 
        # 1. Confomrational independance or 2. Descriptors derived from atoms other than P or Pd
        
        
        ferro_ligs = []
        no_ferro_ligs = []
        all_ligs = []
        for lig in self.ligands:
            
            all_ligs.append(lig)
            if 'erro' in lig:
                ferro_ligs.append(lig)
            else:
                no_ferro_ligs.append(lig)
        
        
        
        
        prop_ranges = []
        label_ranges = []
        for prop in self.properties:
            
            
            
            temp_ranges = []
            
            if ferro:
                ligs = list(set(ferro_ligs))
                
            elif noferro:
                ligs = list(set(no_ferro_ligs))
                
            else:
                ligs = list(set(all_ligs))
                
            
            for lig in ligs:
                data = self.scaled_df[prop].loc[self.scaled_df['Ligand'] == lig]
                vals = [val for val in data if val != 0]
                
                if len(vals) > 1:
                    L_lim = min(vals)
                    H_lim = max(vals)
                    LH_diff = abs(abs(H_lim) - abs(L_lim))
                    temp_ranges.append(LH_diff)
            
            label_ranges.append([prop, temp_ranges])
        
        all_ranges = [len(x[1]) for x in label_ranges]
        min_cutoff = max(all_ranges) - max(all_ranges)*0.25
        for x in label_ranges:
            if len(x[1]) < min_cutoff:
                continue
            c, l = get_label_color(x[0])
            
            lig_range = statistics.mean(x[1])  
            prop_ranges.append([x[0],lig_range,c,l])
        
        
        prop_ranges.sort(key = lambda x: x[1])
        
        with open(f'{outfile}.csv', 'w') as wf:
            for pr in prop_ranges:
                wf.write(f'{pr[0]},{pr[1]}\n')
                
        self.get_feat_range_plot(prop_ranges, outfile)
        

    
    def vplot_BB_per_feature(self, outfig, prop): 
        
        
        plt.rcParams.update({'mathtext.default':  'regular'})
        
        
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20,16), dpi=self.dpi, gridspec_kw={'height_ratios': [12, 1]})
        
        for i, ligand in enumerate(self.ligands):
            x = [self.x_vals[i]+((i+1)*0.05)]
            
            data = self.df[prop].loc[self.df['Ligand'] == ligand]
        
            print(data)
            d = max(data) - min(data)
            
            if d == 0 or len(list(data)) < 2:  
                continue
            
            vplot = ax1.violinplot(data, x, showmeans=False, showmedians=False, showextrema=False)
            
            
            bb = [x[1] for x in self.bridge_ls if x[0] == ligand]
            c = self.ref_cmap[bb[0]]
            c = np.array([c])
                
            
            for cb in vplot['bodies']: 
                cb.set_facecolor(c)
                cb.set_edgecolor('black')
                cb.set_linewidth(1)
                cb.set_alpha(0.8)
            
            
            y = self.weighted_dataset['Boltz_avg'].loc[(self.weighted_dataset['Ligand'] == ligand) & (self.weighted_dataset['Property'] == prop)]
            
            
                
            ax1.scatter(x[0],y,c='black',marker='_',s=[256],label=ligand)
            
            
        
        box_props = dict(boxstyle='round', facecolor='lavender', alpha=0.5, pad=0.5)
        prop = prop.replace('/', '-')
        ax1.set_title(prop, fontsize=30, fontweight='bold', bbox=box_props)
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_ylabel('Feature value', fontsize=30, labelpad=5)
        ax1.set_xlabel('Ligand', fontsize=30, labelpad=5)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        
        plt.set_cmap(self.cmap)
            
        cbar = mpl.colorbar.ColorbarBase(ax2,cmap=self.cmap,norm=self.norm,orientation='horizontal')
        cbar.ax.get_yaxis().set_ticks([])
        cbar.ax.get_xaxis().set_ticks([])
        
        for j, lab in enumerate(self.bb_type):
            x_t = ((2 * j + 1) / 1.94)-0.53
            txt_c = cbar.ax.text(x_t, -0.15, lab, ha='right', va='center', size=20, color='black', rotation=45, rotation_mode='anchor')
            txt_c = cbar.ax.text(x_t, 0, '-', ha='center', va='center', size=30, color='black', rotation=90, rotation_mode='anchor')
            
        
        sp3, sp3d, sp2, sp = "$\it{sp\u00b3}$", "$\it{sp\u00b3}$", "$\it{sp\u00b2}$", "$\it{sp}$"
        
        cbar_title = f'Backbone Type (a=C-{sp3}, b=N-{sp3}, c=O-{sp3}, d=Fe-{sp3}-d, e=C-{sp2}, f=N-{sp2}, g=O-{sp2}, h=C-{sp})'
        cbar.ax.text(17.0, -2.4, cbar_title, ha='center', va='center', size=20, color='black')
        fig.subplots_adjust(bottom=0.15)
        
        plt.savefig(f'{outfig}.png', bbox_inches='tight')
        plt.show()
        
        plt.close()


    def plot_BB_all_ligs(self, properties, dpi):
        self.dpi = dpi
        for prop in self.properties:
            outfig = f'{prop}_BB_plot'
            self.vplot_BB_per_feature(outfig, prop)


    def plot_defined_ligs_props(self, properties, ligands, dpi, outfig):
        self.dpi = dpi
        self.vplot_single(outfig, properties, ligands)
        

    def plot_feat_range(self, dpi, outfile, ferro=False, noferro=False):
        self.get_feature_ranges(outfile, ferro=ferro, noferro=noferro)

"""
dpi = 100
j = graph_ensambles(df=pd.read_csv('Bisphos_all_features.csv'))

## Three plotting options...

## For below plotting, all ligands are used in a single
## plot, one for each specified property 
properties = ['Vbur_5','Atom5','Atom6']
j.plot_BB_all_ligs(properties, dpi)

## For below plotting, define exact properties and ligands
## No backbone key will be given here
properties = ['Vbur_5','Atom5','Atom6']
ligands = ['100_cs', '101_cs']
outfig = 'Sterics_NBO_2_ligands'
j.plot_defined_ligs_props(properties, ligands, dpi, outfig)

## Below plotting is for the feature range bar graph, can
## inlcude/exclude ferrocene structures 
outfile = 'Feature_ranges_all_ligands'
j.plot_feat_range(dpi, outfile, ferro=True, noferro=False)
"""



