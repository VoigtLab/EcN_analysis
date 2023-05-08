import os
import csv
import json
import pandas as pd
import numpy as np
import dnaplotlib as dpl
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import stats
import glob
import pymongo
import datetime

inch_plot_width = 1.022*2*1.18/1.5841 #plot size
lims = [1.0e-4,1.0e4] #plot limits for scatter plots FPKM
profile_plot_margin = 750 #margin for profile plots in bp
plt.rcParams["font.family"] = 'Arial'
plt.rcParams["font.size"] = 6
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
colormap = np.array([[0.89411765,0.10196078,0.10980392,1.],
        [0.21568627,0.49411765,0.72156863,1.],
        [0.30196078,0.627451,0.29019608,1.],
        [1.,0.6803922,0.25,1.],
        [0.65098039,0.3372549,0.15686275,1.],
        [0.96862745 ,0.50588235,0.74901961,1.],
        [0.59607843 ,0.30588235,0.63921569,1.],
        [0.3,0.3,0.3,1.]])

username = 'hdoosth'
# Location for files and outputs
input_folder = '2021_05_23_EcN_input'
output_folder = '2021_05_23_EcN_output'
root = os.getcwd()
default_prefix = "/home/project_source/"

regions_of_int = {'AJT206_EcN_SensorArrayOnly':{'LP3':[4226580,4229433],'LP2':[4307675,4319200],'LP1':[4641508,4644416]},
    'EcN_WT':{'LP3':[4226580,4229081],'LP2':[4307323,4309824],'LP1':[4632132,4634633]}}
    # insertion sites: 4227830-1, 4308573-4, 4633382-3

# Metadata and preCAD FPKM data file
metad_file = 'experiment.ginkgo.42587_QC_and_metadata.csv'
data_file = 'experiment.ginkgo.42587_ReadCountMatrix_preCAD_FPKM.csv'
profile_dataframe_file = '2021_01_08_EcN_hdoost_h_updated_dataframe.csv'

# Data loading (alldata includes metadata and FPKM values)
metadata = pd.read_csv(os.path.join(root,input_folder,metad_file), low_memory=False)
alldata = pd.read_csv(os.path.join(root,input_folder,data_file), low_memory=False)
profile_dataframe = pd.read_csv(os.path.join(root,input_folder,profile_dataframe_file), low_memory=False)

# Separate metadata from main dataframe so it only contains gene FPKM values
data = alldata[~alldata["row_titles"].isin(metadata["row_titles"].tolist())]
metadata.set_index("row_titles", inplace=True)
data.set_index("row_titles", inplace=True)

def mongo_query(experiment_id):
    
    # Mongo settings:
    dbURI = '' # datacatalog entered here
    client = pymongo.MongoClient(dbURI)
    db = client.catalog_staging
        
    # Select a collection
    jobs = db.collection1
    query_sam={}
    query_sam['data.experiment_id'] = {'$in': experiment_id}
    preprocessing_jobs = {}
    alignment_jobs = {}
    
    # Create a list of the query results
    for job in jobs.find(query_sam):
        sample_id = job['data']['sample_id']
        if job["pipeline_id"] == "preprocessing_id":
            preprocessing_jobs[sample_id] = job
        if job["pipeline_id"] == "alignment_id":
            alignment_jobs[sample_id] = job
    print(len(preprocessing_jobs), " pre-processing jobs and ", len(alignment_jobs), " alighment jobs found in mongo query for ", experiment_id)

    # Create dictionary for sam files:
    sam_dict = {}
    all_dict = {}
    for sample_id, job in alignment_jobs.items():

        all_outputs = glob.glob(default_prefix + job['archive_path'] + '/*')
        # For just the sam file we can use the extension:
        sam_file = glob.glob(default_prefix + job['archive_path'] + '/*.sorted.sam')
        if len(sam_file) > 1:
            print("Warning in sample ", sample_id, "- more than one .sorted.sam file found")
        elif len(sam_file) == 0:
            print("Warning in sample", sample_id, "- no .sorted.sam file found")
        else:
            sam_dict[sample_id] = sam_file[0].split(default_prefix)[1]
        all_dict[sample_id] = all_outputs

    return sam_dict, all_dict, alignment_jobs, preprocessing_jobs

def parse_gff3(file_name):
    cols = ['genome','part_type','start','end','strand','Name','label','note','product','protein_id','translation']
    all_parts = {}
    with open(file_name) as file:
        readcsv = csv.reader(file,delimiter='\t')
        for row in readcsv:
            if not row[0].startswith('##'):
                if row[2]=='CDS':
                    info = row[8].split(';')
                    info_dict = {detail.split('=')[0]:detail.split('=')[1] for detail in info}
                    info_dict['genome']=row[0]
                    info_dict['part_type']=row[2]
                    info_dict['start']=float(row[3])
                    info_dict['end']=float(row[4])
                    info_dict['strand']=row[6]
                    all_parts[info_dict['Name']]=[info_dict[x] if x in info_dict.keys() else "NA" for x in cols]
    gff_df = pd.DataFrame.from_dict(all_parts, orient='index', columns=cols)
    return gff_df

def create_profile(file_path, file_name, prefix=os.path.join(root,input_folder), output_path=os.path.join(root,input_folder), sample_id=""):
    SAM_file = prefix + file_path + file_name
    reads = {}
    profile = {}
    genomes = {}
    nline = 0
    mappedline = 0
    min_length = 10
    stats = {}
    max_values = {}
    with open(SAM_file, 'rU') as ins:
        for lines in ins:
            field = lines.strip().split()
            if lines[0] == '@':
                if field[0] == '@SQ':
                    name = field[1].split(':')
                    length = field[2].split(':')
                    genomes[name[1]] = length[1]
                    profile[name[1] + '_fwd'] = np.zeros(int(length[1]))
                    profile[name[1] + '_rev'] = np.zeros(int(length[1]))
                    reads[name[1] + '_fwd'] = {}
                    reads[name[1] + '_rev'] = {}
            if lines[0] != '@':
                if field[2] in genomes and int(field[4]) >= 10:                     
                    start = min(int(field[3]),int(field[7]))
                    readlen = abs(int(field[8]))-1
                    if int(field[1]) == 83 or int(field[1]) == 163:
                        if readlen >= min_length:
                            if start not in reads[field[2]+'_fwd']:
                                reads[field[2]+'_fwd'][start] = []
                            reads[field[2]+'_fwd'][start].append(readlen)
                            profile[field[2]+'_fwd'][start:start+readlen] += np.ones(readlen)
                            mappedline += 1
                    elif int(field[1]) == 99 or int(field[1]) == 147:
                        if readlen >= min_length:
                            if start not in reads[field[2]+'_rev']:
                                reads[field[2]+'_rev'][start] = []
                            reads[field[2]+'_rev'][start].append(readlen)
                            profile[field[2]+'_rev'][start:start+readlen] += np.ones(readlen)
                            mappedline += 1
                nline += 1
    for seq, length in genomes.items():
        max_values[seq+'_fwd'] = max(profile[seq+'_fwd'])
        max_values[seq+'_rev'] = max(profile[seq+'_rev'])
    print(str(nline) + " reads found in SAM file " + str(mappedline) + " met qc to be mapped into profile")
    stats['reads_sam'] = nline
    stats['reads_profile'] = mappedline
    stats['max_values'] = max_values
    stats['genomes'] = genomes
    print("Profile Generation Complete - Saving")
    json_reads = json.dumps(reads)
    with open(output_path+sample_id+"_reads.json","w") as f_reads:
        f_reads.write(json_reads)
    profile_for_save = {}
    for key, length in genomes.items():
        profile_for_save[key + '_fwd'] = profile[key + '_fwd'].tolist()
        profile_for_save[key + '_rev'] = profile[key + '_rev'].tolist()
    json_profile = json.dumps(profile_for_save)
    with open(output_path+sample_id+"_profile.json","w") as f_profile:
        f_profile.write(json_profile)
    print("Save Complete")
    return profile, stats

def analyze_all_samples():
    (sam_dict, all_outputs, alignment_jobs, preprocessing_jobs) = mongo_query(['experiment.ginkgo.42587'])
    output_path = os.path.join(root,input_folder)
    df = pd.read_csv('2021_01_07_EcN_hdoost_h_dataframe.csv', header=3, low_memory=False)
    n = df.shape[0]
    m = df.shape[1]
    now = datetime.datetime.now()
    df.insert(m, 'reads_file', ['']*n, True)
    df.insert(m+1, 'profile_file', ['']*n, True)
    df.insert(m+2, 'stats', ['']*n, True)
    df.insert(m+3, 'profile_flag', [0]*n, True)
    df.insert(m+4, 'chrom', ['']*n, True)
    count_samples = 0

    for sample_id,bwa_archive_path in sam_dict.items():
        print('\n', sample_id)
        if os.path.isfile(output_path+sample_id+"_reads.json") and  os.path.isfile(output_path+sample_id+"_profile.json"):
            a = df.loc[df['sample_id'] == sample_id, 'flag_profile'] = 1
            print("File exists - assuming previously completed successfully")
        else:
            (profile, stats) = create_profile(bwa_archive_path, '',sample_id=sample_id)
            a = df.loc[df['sample_id'] == sample_id, 'stats'] = str(stats)
            a = df.loc[df['sample_id'] == sample_id, 'flag_profile'] = 1
        a = df.loc[df['sample_id'] == sample_id, 'reads_file'] = output_path+sample_id+"_reads.json"
        a = df.loc[df['sample_id'] == sample_id, 'profile_file'] = output_path+sample_id+"_profile.json"
        dict_chroms = {}
        for key in profile.keys():
            if key[-4::] == '_fwd':
                chrom = key[0:-4]
                dict_chroms[chrom] = len(profile[key])
        a = df.loc[df['sample_id'] == sample_id, 'chrom'] = json.dumps(dict_chroms)
        count_samples += 1
        updated_dataframe = output_path+now.strftime("%Y_%m_%d__%H_%M_")+username+'updated_dataframe.csv'
    df.to_csv(updated_dataframe, index=False)
    print(count_samples, " samples analyzed and profiles generated, metadata has been updated with addresses")

def create_normalized_profiles(profile_dataframe):
    ## Demonstration and manual checks
    print(len(profile_dataframe[profile_dataframe['measurement_type']=='RNA_SEQ']['sample_id'].unique()))  # prints number of unique samples
    print(profile_dataframe[profile_dataframe['sample_id']=='sample.ginkgo.30832357.experiment.ginkgo.42587'].iloc[0]) # prints metadata available for a single sample
    print(profile_dataframe[profile_dataframe['sample_id']=='sample.ginkgo.30832357.experiment.ginkgo.42587'].iloc[0]['IPTG_concentration']) # prints IPTG concentrations for said sample
    print(json.loads(profile_dataframe[profile_dataframe['sample_id']=='sample.ginkgo.30832357.experiment.ginkgo.42587'].iloc[0]['chrom'])) # prints dictionary of chromosomes in associated data files

    os.path.join(root,input_folder,metad_file)
    ## Create all the .bed files usable by DNAplotlib. *_profile.json files are raw data and *.bed files are normalized
    ## Warning: the following is only suitable for single chromosome data files (no plasmids).
    for sample_id in profile_dataframe[profile_dataframe['measurement_type']=='RNA_SEQ']['sample_id'].unique():
        if os.path.isfile(os.path.join(root,input_folder,sample_id+"_profile.json")):
            if not os.path.isfile(os.path.join(root,input_folder,sample_id+"_fwd.bed")) or not os.path.isfile(os.path.join(root,input_folder,sample_id+"_rev.bed")):
                with open(os.path.join(root,input_folder,sample_id+"_profile.json"),"rU") as f_profile:
                    profile_load = json.load(f_profile)
                profile_raw={} # raw profile values
                profile = {} # normalized profile values
                profile_for_save = {}
                for key in profile_load.keys():
                    if key[-4::] == '_fwd':
                        genome = key[0:-4]
                        profile_raw[genome + '_fwd'] = np.array(profile_load[genome + '_fwd'])
                        profile_raw[genome + '_rev'] = np.array(profile_load[genome + '_rev'])
                        total_reads = sum(profile_raw[genome + '_fwd']) + sum(profile_raw[genome + '_rev'])
                        profile[genome + '_fwd'] = profile_raw[genome + '_fwd']/total_reads*1e9
                        profile[genome + '_rev'] = profile_raw[genome + '_rev']/total_reads*1e9
                        c=1
                        with open(os.path.join(root,input_folder,sample_id+"_fwd.bed"),"w") as fwd_bed_file:
                            for ii in profile[genome+'_fwd']:
                                fwd_bed_file.write(key+'\t0\t'+str(len(profile[genome + '_fwd']))+'\t'+str(c)+'\t'+str(ii)+'\n') 
                                c = c+1
                        e=1
                        with open(os.path.join(root,input_folder,sample_id+"_rev.bed"),"w") as rev_bed_file:
                            for jj in profile[genome+'_rev']:
                                rev_bed_file.write(key+'\t0\t'+str(len(profile[genome + '_rev']))+'\t'+str(e)+'\t'+str(jj)+'\n') 
                                e = e+1
                        profile_for_save[genome + '_fwd'] = profile[genome + '_fwd'].tolist()
                        profile_for_save[genome + '_rev'] = profile[genome + '_rev'].tolist()
                        json_profile = json.dumps(profile_for_save)
                        with open(os.path.join(root,input_folder,sample_id+"_profile_norm.json"),"w") as f_profile_norm:
                            f_profile_norm.write(json_profile)
                kstr = str(list(profile_load.keys()))
                print("BED files created successfully for sample: " + sample_id)
            else:
                print("BED files for sample: " + sample_id + " already exist")
        else:
            print("No _profile.json file found for sample: " + sample_id)

def plot_annotated_profiles(profile_dataframe,cur_region):

    ## The following will loop through all:
    for sample_id in profile_dataframe[profile_dataframe['measurement_type']=='RNA_SEQ']['sample_id'].unique():
        if os.path.isfile(os.path.join(root,input_folder,sample_id+"_fwd.bed")) and os.path.isfile(os.path.join(root,input_folder,sample_id+"_rev.bed")):
            fwd_bed_filename = os.path.join(root,input_folder,sample_id+"_fwd.bed")
            rev_bed_filename = os.path.join(root,input_folder,sample_id+"_rev.bed")
            
            ## The following would loop through all chromosomes (there's only one in these cases)
            for chrom, chrom_len in json.loads(profile_dataframe[profile_dataframe['sample_id']==sample_id].iloc[0]['chrom']).items(): 
                print(sample_id, chrom, chrom_len)
                gff_file = os.path.join(root,input_folder,chrom+'.gff')

                # Load the design from a GFF file (see DNAplotlib github example for further details)
                design = dpl.load_design_from_gff(gff_file, chrom, region=cur_region)
                profile_fwd = dpl.load_profile_from_bed(fwd_bed_filename, chrom+'_fwd', [0, chrom_len-1])
                profile_rev = dpl.load_profile_from_bed(rev_bed_filename, chrom+'_rev', [0, chrom_len-1])

                # Create the DNAplotlib renderer
                dr = dpl.DNARenderer(scale=10.0)
                part_renderers = dr.trace_part_renderers()

                # Create the figure
                fig = plt.figure(figsize=(3.5,2.0))
                gs = gridspec.GridSpec(2, 1, height_ratios=[1, 0.2])
                ax_dna = plt.subplot(gs[1])

                # Redender the DNA to axis
                start, end = dr.renderDNA(ax_dna, design, part_renderers)
                ax_dna.set_xlim(cur_region)
                ax_dna.set_ylim([-5,8])
                ax_dna.axis('off')

                ax_profile = plt.subplot(gs[0])
                ax_profile.fill_between(list(range(cur_region[0],cur_region[1])), profile_fwd[cur_region[0]:cur_region[1]], np.zeros(cur_region[1]-cur_region[0]), color=(0.5,0.5,0.5), edgecolor=(0.5,0.5,0.5), linewidth=1, zorder=1.5)
                ax_profile.fill_between(list(range(cur_region[0],cur_region[1])), np.array(profile_rev[cur_region[0]:cur_region[1]])*-1.0, np.zeros(cur_region[1]-cur_region[0]), color=(1,0,0), edgecolor=(1,0,0), linewidth=1, zorder=1.5)
                ax_profile.plot(cur_region, [0,0], color=(0,0,0), zorder=10)
                ax_profile.set_xlim(cur_region)
                ax_profile.axis('off')

                # Update subplot spacing
                plt.subplots_adjust(hspace=0.001, left=0.01, right=0.99, top=0.99, bottom=0.01)

                # Save the figure
                fig.savefig(os.path.join(output_folder, sample_id + '_' + str(cur_region[0]) + '_' + str(cur_region[1]) + '.png'), dpi=300)

                # Clear the plotting cache
                plt.close('all')
    
def plot_regions_of_interest(profile_dataframe):

    for strain in np.unique(metadata.loc['Strain'].values).tolist():
        chrom_dict = json.loads(profile_dataframe[profile_dataframe['strain']==strain].iloc[0]['chrom'])
        for chrom in chrom_dict.keys():
            for region in regions_of_int[chrom].keys():
                cur_region = regions_of_int[chrom][region]
                cur_region[0]-=profile_plot_margin
                cur_region[1]+=profile_plot_margin
                for timepoint in np.unique(metadata.loc['Timepoint'].values).tolist():
                    for iptg in np.unique(metadata.loc['IPTG'].values).tolist():
                        
                        print(strain, timepoint, iptg)
                        # Create the figure
                        fig = plt.figure(figsize=(6.0*(cur_region[1]-cur_region[0])/3000,4.0))
                        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 0.2])
                        ax = plt.subplot(gs[0])
                        plt.title(chrom + ' ' + timepoint + 'h IPTG=' + iptg + 'M ' + region)

                        x_conds = metadata.loc['Strain'].eq(strain) & metadata.loc['Timepoint'].eq(timepoint) & metadata.loc['IPTG'].eq(iptg) & metadata.loc['QC_mapped_reads_BOOL'].eq('TRUE')
                        x_meta = x_conds[x_conds].index
                        print(x_meta)
                        avg_profile_fwd = np.zeros(len(range(cur_region[0],cur_region[1])))
                        avg_profile_rev = np.zeros(len(range(cur_region[0],cur_region[1])))
                        samp_count = 0
                        for sample_id in x_meta:
                            if os.path.isfile(os.path.join(root,input_folder,sample_id+"_profile_norm.json")):
                                profile ={}
                                with open(os.path.join(root,input_folder,sample_id+"_profile_norm.json"),"r") as f_profile:
                                    profile_load = json.load(f_profile)
                                for key in profile_load.keys():
                                    if key[-4::] == '_fwd':
                                        genome = key[0:-4]
                                        profile[genome + '_fwd'] = np.array(profile_load[genome + '_fwd'])
                                        profile[genome + '_rev'] = np.array(profile_load[genome + '_rev'])
                            else:
                                print("Warning: Normalized Profiles haven't been created") 
                                print("Use create_normalized_profiles")
                                return
                            for key in profile.keys():
                                if key[-4::] == '_fwd':
                                    genome = key[0:-4]
                                    super_threshold_indices = profile[genome + '_fwd'] < 1
                                    profile[genome + '_fwd'][super_threshold_indices] = 1
                                    super_threshold_indices = profile[genome + '_rev'] < 1
                                    profile[genome + '_rev'][super_threshold_indices] = 1
                                    avg_profile_fwd += np.log10(profile[genome + '_fwd'][cur_region[0]:cur_region[1]])
                                    avg_profile_rev += -np.log10(profile[genome + '_rev'][cur_region[0]:cur_region[1]])
                                    samp_count +=1
                            del profile
                        ax.fill_between(np.array(range(cur_region[0],cur_region[1])), 0.0, avg_profile_fwd/samp_count, facecolor='grey', linewidth=0.0)
                        ax.fill_between(np.array(range(cur_region[0],cur_region[1])), 0.0, avg_profile_rev/samp_count, facecolor='lightpink', linewidth=0.0)
                        ax.axvline(cur_region[0]+profile_plot_margin,-5,5,color='maroon',linestyle='--',linewidth=0.5)
                        ax.axvline(cur_region[1]-profile_plot_margin,-5,5,color='maroon',linestyle='--',linewidth=0.5)

                        ax.plot(cur_region,[0,0],color='k',linewidth=1.0)
                        ax.set_xlim(cur_region)
                        ax.ticklabel_format(style='plain')
                        ax.set_ylim([-5,5])
                        for axis in ['bottom','top','right']:
                            ax.spines[axis].set_linewidth(1.0)
                            ax.spines[axis].set_color('k')
                        ax.spines['left'].set_linewidth(1.0)
                        ax.spines['left'].set_color('k')
                        
                        ax.tick_params(axis='y', which='major', labelsize=7, length=4, width=0.25, direction='out')
                        ax.tick_params(axis='x', which='major', labelsize=7, length=4, width=0.25, direction='out')
                        ax.set_xticks([cur_region[0],cur_region[0]+500,cur_region[0]+1000,cur_region[0]+1500,cur_region[1]-1500,cur_region[1]-1000,cur_region[1]-500,cur_region[1]])
                        ax.set_xticklabels(["-1500","-1000","-500","0","0","500","1000","1500"])
                        ax.set_yticks([-5,0,5])
                        ax.set_yticklabels([-5,0,5])

                        ax_dna = plt.subplot(gs[1])
                        gff_file = os.path.join(root,input_folder,chrom + ".gff")
                        # print(design) 

                        # Create the DNAplotlib renderer
                        design = dpl.load_design_from_gff(gff_file, chrom, region=cur_region)
                        dr = dpl.DNARenderer(scale=10.0)
                        part_renderers = dr.trace_part_renderers()
                        start, end = dr.renderDNA(ax_dna, design, part_renderers)
                        ax_dna.set_xlim(cur_region)
                        ax_dna.set_ylim([-8,8])
                        ax_dna.axis('off')

                        # Update subplot spacing
                        plt.subplots_adjust(hspace=0.2, left=0.15, right=0.85, top=0.85, bottom=0.15)
                        
                        # Save the figure
                        fig.savefig(os.path.join(output_folder, sample_id + '_' + str(cur_region[0]) + '_' + str(cur_region[1]) + '.png'), dpi=300)

                        # Clear the plotting cache
                        plt.close('all')
                        plt.clf()
                        plt.cla()
        del chrom_dict
        
def plot_diff_ex(x_vals, y_vals, title='2D Scatter Plot of Expression Levels', x_label='X', y_label='Y',strain=""):
    gs = gridspec.GridSpec(1,1)
    fig1 = plt.figure(figsize=(inch_plot_width,inch_plot_width))

    ax1 = plt.subplot(gs[0,0])

    thresh_x = np.geomspace(lims[0],lims[1],1000)
    thresh_y = np.maximum(thresh_x*10,thresh_x+0.5)
    indices = thresh_y<=lims[1]
    thresh_x = thresh_x[indices]
    thresh_y = thresh_y[indices]
    for i in range(len(y_vals)):
        ax1.scatter(x_vals,y_vals[i],s=1,color=colormap[i],marker='.',alpha=0.75,linewidth=0,zorder=10)
    ax1.plot(thresh_x,thresh_y,color='darkgrey',linewidth=0.5,zorder=5)
    ax1.plot(thresh_y,thresh_x,color='darkgrey',linewidth=0.5,zorder=5)

    ax1.set_yscale('log')
    ax1.set_ylim(lims)
    ax1.set_xscale('log')
    ax1.set_xlim(lims)
    for axis in ['bottom', 'left','top','right']:
        ax1.spines[axis].set_linewidth(1)
        ax1.spines[axis].set_color('black')
    ax1.tick_params(axis='both', which='major', labelsize=8, length=3, width=0.5, direction='in', pad=2)
    ax1.tick_params(axis='both', which='minor', labelsize=8, length=2, width=0.5, direction='in')
    locmaj = matplotlib.ticker.LogLocator(base=10,numticks=10)
    ax1.yaxis.set_major_locator(locmaj)
    ax1.xaxis.set_major_locator(locmaj)
    locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),numticks=100)
    ax1.yaxis.set_minor_locator(locmin)
    ax1.xaxis.set_minor_locator(locmin)
    ax1.plot(lims,lims,color='k',linewidth=0.5)
    ax1.set_ylabel(y_label)
    ax1.set_xlabel(x_label)
    ax1.set_title(title)

    plt.savefig(os.path.join(output_folder,strain+", DE_"+title+".png"),dpi=300)
    plt.clf()
    plt.close()

def diff_ex_all():
    gff_file = "EcN_WT.gff"
    df = parse_gff3(os.path.join(root,input_folder,gff_file))

    ## Analyze for every timepoint, induction condition, and strain separately using appropriate reference strain
    for timepoint in np.unique(metadata.loc['Timepoint'].values).tolist():
        for iptg in np.unique(metadata.loc['IPTG'].values).tolist():
            for strain in np.unique(metadata.loc['Strain'].values).tolist():
                if strain == 'EcN_SensorArrayOnly':
                    ref_strain = 'Escherichia_coli_Nissle_WT'
                else:
                    ref_strain = 'EcN_SensorArrayOnly'

                ## Select samples (replicates) which pass QC and correspond to the selected conditions
                x_conds = metadata.loc['Strain'].eq(ref_strain) & metadata.loc['Timepoint'].eq(timepoint) & metadata.loc['IPTG'].eq(iptg) & metadata.loc['QC_mapped_reads_BOOL'].eq('TRUE')
                x_meta = metadata[x_conds[x_conds].index]
                x_data = data[x_conds[x_conds].index].astype(float)
                y_conds = metadata.loc['Strain'].eq(strain) & metadata.loc['Timepoint'].eq(timepoint) & metadata.loc['IPTG'].eq(iptg) & metadata.loc['QC_mapped_reads_BOOL'].eq('TRUE')
                y_meta = metadata[y_conds[y_conds].index]
                y_data = data[y_conds[y_conds].index].astype(float)

                ## Calculate statistics and fold change each gene
                stat, pvalue = stats.ttest_ind(x_data.T,y_data.T,equal_var=False)
                pv = pd.Series(pvalue,index=y_data.index)
                x_vals = (x_data.mean(1),x_data.std(1))
                y_vals = (y_data.mean(1),y_data.std(1))
                reps = y_data.columns
                fc=y_vals[0]/x_vals[0]
                sig_indices = np.all([pvalue<0.05,np.any([fc>=10,fc<=0.1],axis=0)],axis=0)

                ## append mean FPKM (_values), fold change (_fc), coefficien of variance (_cv), p values (_pv) for each strain/condition to data frame
                if strain == 'EcN_SensorArrayOnly':
                    df_wt = pd.DataFrame({ref_strain+"_IPTG="+iptg+"M_values": x_vals[0]})
                    df_sens = pd.DataFrame({strain+"_IPTG="+iptg+"M_values": y_vals[0]})
                    df = pd.concat([df, df_wt], axis=1)
                    df = pd.concat([df, df_sens], axis=1)
                df_right = pd.DataFrame({strain+"_IPTG="+iptg+"M_fc": (y_vals[1]/y_vals[0])})
                df_right1 = pd.DataFrame({strain+"_IPTG="+iptg+"M_cv": (y_vals[1]/y_vals[0])})
                df_right2 = pd.DataFrame({strain+"_IPTG="+iptg+"M_pvalues": pv})
                df = pd.concat([df, df_right], axis=1)
                df = pd.concat([df, df_right1], axis=1)
                df = pd.concat([df, df_right2], axis=1)

                ## Replace zeros with 3e-4 for plotting in logspace
                x_vals[0].replace(0,3e-4,inplace=True)
                y_vals[0].replace(0,3e-4,inplace=True)
                plot_diff_ex(x_vals, y_vals, title=timepoint+"h, IPTG="+iptg+"M",x_label="RNAP Flux [FPKM]\n"+ref_strain,y_label=strain+"\n RNAP Flux [FPKM]", strain=strain)

    df.to_csv(os.path.join(output_folder,"2021_04_09_CV_all_23h.csv"))

# analyze_all_samples()
diff_ex_all()
plot_regions_of_interest(metadata)
