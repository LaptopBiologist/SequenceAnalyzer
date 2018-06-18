#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      I am
#
# Created:     22/11/2016
# Copyright:   (c) I am 2016
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import sklearn
from mpl_toolkits import mplot3d
import matplotlib
from matplotlib import pyplot
import matplotlib.gridspec as gridspec

import sklearn.decomposition
import sklearn.neighbors
import sklearn.linear_model
import sklearn.model_selection
import sklearn.ensemble
import statsmodels.sandbox.stats
import statsmodels.stats
import statsmodels.stats.multitest


import seaborn
import numpy
import scipy
import xlsxwriter


import os

from matplotlib.widgets import Lasso
from matplotlib.colors import colorConverter
from matplotlib.collections import RegularPolyCollection
from matplotlib import path

import matplotlib.pyplot as plt
from numpy import nonzero
from numpy.random import rand
try:
    import wpca
except:
    pass
##matplotlib.interactive(True)

colors=numpy.array( ['firebrick', 'firebrick', 'firebrick', 'firebrick', 'firebrick', 'firebrick', 'firebrick', 'firebrick', 'firebrick', 'firebrick', 'firebrick', 'firebrick', 'firebrick', 'firebrick', 'firebrick', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'forestgreen', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'blue', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange', 'orange', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple'])


bad_strains=['B18', 'ZS29' ,'ZW190', 'B17','T26A','ZW133', 'ZW149']

samples=['B04', 'B05', 'B10', 'B11', 'B12', 'B14', 'B17', 'B18', 'B23', 'B28', 'B38', 'B42', 'B43', 'B51', 'B52', 'B54', 'B59', 'I01', 'I02', 'I03', 'I04', 'I06', 'I07', 'I13', 'I16', 'I17', 'I22', 'I23', 'I24', 'I26', 'I29', 'I31', 'I33', 'I34', 'I35', 'I38', 'N01', 'N02', 'N03', 'N04', 'N07', 'N10', 'N11', 'N13', 'N14', 'N15', 'N16', 'N17', 'N18', 'N19', 'N22', 'N23', 'N25', 'N29', 'N30', 'T01', 'T04', 'T05', 'T07', 'T09', 'T10', 'T14A', 'T22A', 'T23', 'T24', 'T25A', 'T26A', 'T29A', 'T30A', 'T35', 'T36B', 'T39', 'T43A', 'T45B', 'ZH23', 'ZH26', 'ZH33', 'ZH42', 'ZS10', 'ZS29', 'ZW09', 'ZW133', 'ZW139', 'ZW140', 'ZW142', 'ZW144', 'ZW149', 'ZW155', 'ZW177', 'ZW184', 'ZW185', 'ZW190']
##
good_ind=numpy.array( [True, True, True, True, True, True, False, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, False, True, False, True, True, True, True, False, True, True, True, True, False])
refnames=['3r', '3l', '2r', '2l', 'x', '4']

ID_Header=['Subject', 'Target', 'Quad X', 'Quad Y', 'Feature Type']
Jxn_Header=['Jxn X', 'Jxn X (range)','Jxn Y', 'Jxn X (range)']
Position_Header=["Read Count","Reads 3' X", "Reads 3' Y"]

repeat_labels=['rDNA', 'Histone', 'NTS', 'DM359', 'R1', 'Bari']
#Style
##matplotlib.colors.cnames['firebrick']=u'#6f0000'
##matplotlib.colors.cnames['forestgreen']=matplotlib.colors.rgb2hex((0, 73./256,73/256.))
##matplotlib.colors.cnames['orange']=matplotlib.colors.rgb2hex((0, 209./256,219/256.))

# This sets reasonable defaults for font size for
# a figure that will go in a paper
seaborn.set_context("paper")

# Set the font to be serif, rather than sans
seaborn.set(font='serif')

# Make the background white, and specify the
# specific font family
seaborn.set_style("white", {
    "font.family": "serif",
    "font.serif": ["Times","Roman", "serif"]
})

matplotlib.rcParams['svg.fonttype'] = 'none'
##matplotlib.rcParams['eps.fonttype'] = 'none'


def set_style():
    # This sets reasonable defaults for font size for
    # a figure that will go in a paper
    seaborn.set_context("paper")

    # Set the font to be serif, rather than sans
    seaborn.set(font='serif')

    # Make the background white, and specify the
    # specific font family
    seaborn.set_style("white", {
        "font.family": "serif",
        "font.serif": ["Times","Roman", "serif"]
    })

    matplotlib.rcParams['svg.fonttype'] = 'none'


class SequenceSummary():
    def __init__(self, infile,ignore=[],sample_file='', pop_rule=1 , FDR=.01, error_rate=1e-3, depth_cutoff=1,mean_cutoff=0, nan_filter=True):
        """Computes two statistics from the SNP pileup.
        First, it defines the major allele at each positions as the allele present at
        the highest average allele proportion across the populations. It

        Second, it estimates the proportion of samples that contain a variant allele at a
        position in at least one repeat unit. Variant alleles are called by rejecting
        the null hypothesis that all observed variants are sequencing errors using a
        conservative error rate epsilon (the minimum quality score). This does not account for
        errors introduced by PCR, which may be important. For each position i in
        sample j total read count X_ij and variant read count k_ij are determined. We
        define the variant allele as the allele with the second greatest average allelw
        proportion. The null model is:

        k_ij ~ Binom(X_ij, epsilon)

        A p-value is computed for each sample and the null hypothesis is (or isn't)
        rejected using the Benjami-Hochsberg multiple testing correction for the
        number of samples at the given FDR (default .01). Note, this null model does
        not account for PCR-duplicates, and the assumption of a simple sampling
        distribution is usually simplistic.

        Parameters:
            snp_array: the I x J x K 3-dimensional snp pileup array

            FDR: Accepted false-discovery rate for a position.

            error_rate: The expected rate of sequencing errors. For base-calling errors,
            the default is 1e-3, which corresponds to lowest Phred quality score considered
            when ConTExt builds SNP pileups (30).

            depth_cutoff: A statistic will only be computed for a consensus position
            if the read depth at the position is higher than the cutoff in every sample.

        Returns:
            An object with the following attributes:
            major_allele_proportion: An I x K array indicating the major allele proportion
            at positon k in sample i

            allele_frequency: A length K array containing the fraction of samples containing a
            variant allele in at least one repeat unit for each position

            passed_positions: The indices of consensus allele positions which survive
            the read depth cutoff

            major_allele: The major allele at each position"""
        #Determine the parent directory
        if sample_file=='':
            indir='/'.join( infile.split('/')[:-1])
            sample_file=indir+'/samples.tsv'

        assert os.path.exists(sample_file)==True, "This function requires a comma separated list of samples labels \
in the same order as in the pile-up table. If the input's parent directory \
does not contain a such a file named samples.tsv, please specify the \
path to such a file in the sample_file parameter."

        self.samples=numpy.loadtxt(sample_file, str)

        print len( self.samples)
        snp_array=numpy.load(infile)
        if type( ignore)==str:
            assert os.path.exists(sample_file)==True, "A string was passed to the 'ignore' parameter,\
            but this was not a valid path."
            ignore=numpy.loadtxt(ignore, str)
        print self.samples
        if ignore!=[]:
            good_indices=~numpy.isin(samples, ignore)
            snp_array=snp_array[good_indices]
            self.samples=self.samples[good_indices]
    #Check for samples without the repeat
        read_counts=numpy.mean( snp_array.sum(1),1)
        good_indices=read_counts>=mean_cutoff
        snp_array=snp_array[good_indices]
        self.samples=self.samples[good_indices]
        if type( pop_rule)==str:
            #File path
            pass
        self.pop=numpy.array( [s[:pop_rule] for s in self.samples])
        print len(self.pop)
        num_samples=snp_array.shape[0]
        read_count=snp_array[:,:,:] .sum(1)
    ##    print read_count.shape


        major_allele_prop,major_allele, all_pos=GetMajorAlleleFreq(snp_array, cutoff=depth_cutoff)
        var_allele_prop=GetOtherAlleleProportion(snp_array, depth_cutoff, order=2)
        read_count=read_count[:,all_pos]
        var_counts=(read_count*var_allele_prop).astype(int)
        read_count=read_count.astype(int)
        min_cvg= numpy.min(read_count,0)
        allele_freq=[]
        var_pres=[]

        for i in range(read_count.shape[1]):
            p_vals=(1.- scipy.stats.binom.cdf(var_counts[:,i], read_count[:,i], error_rate))
    ##        p_vals=(1.- scipy.stats.binom.cdf(snp_array[:,:,i], read_count[:,i][:,None], error_rate))
    ##        print p_vals.shape
            rej,q_vals,sidak, bonf=statsmodels.stats.multitest.multipletests(p_vals,FDR, method='fdr_bh')
            allele_freq.append(float( rej.sum()) /num_samples)
            var_pres.append(rej)

        self.pileup=snp_array
        print self.pileup.shape
        self.major_allele_proportion = major_allele_prop
        self.nan_filter=~numpy.isnan(self.major_allele_proportion.sum(0))
##        print self.pop
        self.allele_frequency = numpy.array( allele_freq)
        self.variant_table = numpy.array( var_pres).transpose()
        self.passed_positions =  all_pos
        self.major_allele = major_allele
##        self.FST=ComputeFst(self.pileup, self.pop)

    def ComputeWeights(self, interpolate=True):
        N=self.pileup.sum(1)
        k=N*self.major_allele_proportion
        a=1+k
        b=1+N-k
        variance=(a*b)/((a+b+1)*(a+b)**2)
        p=self.major_allele_proportion
        q=1.+p
        se=(1/N)**.5

        self.weights=1./variance**.5
##        self.weights=1./se
        nan_ind=numpy.isnan(self.major_allele_proportion)
        self.weights[nan_ind]=0.
        nan_ind=numpy.isnan(self.weights)
        inf_ind=numpy.isinf(self.weights)
        self.weights[inf_ind]=0.
        if interpolate==True:
            self.weights[ self.weights!=0.]=1.

    def PerformPCA(self, method="PCA", interpolate=True, colors=colors, max_iter=100, random_state=100, nan_filter=False,selection_method=None):
    ##    data=data[numpy.array( good_ind),:]
        self.colors=COlorBlindColors( colors)
        data=self.major_allele_proportion
        pos=self.passed_positions
        if nan_filter==True:
            #We are worried about nans
            data=data[:,self.nan_filter]
        if samples!='':
            print 'b'
            keep_ind=numpy.array( [bad_strains.count(s)==0 for s in self.samples])
            colors=numpy.array(colors)
            data=data
            nan_ind=numpy.isnan(data.sum(0))

##            data=data[:,~nan_ind]
    ##        pca=sklearn.decomposition.IncrementalPCA(92-len(bad_strains),whiten=True, copy=False,batch_size= 100)
        if method=='PCA':

            pca=sklearn.decomposition.PCA(data.shape[0],whiten=True, copy=False)
            transformed=pca.fit_transform(AlleleRescale( data))
        else:
            #compute weights
            self.ComputeWeights(interpolate)
            if method=='EMPCA':
                pca=wpca.EMPCA(10,max_iter=max_iter,random_state=random_state)
            elif method=='WPCA':
                pca=wpca.WPCA(10)
            rescaled=AlleleRescale( data)
            sample_size=self.weights.shape[0]
            isnan=numpy.isnan(rescaled.sum(0))*(( self.weights==0).sum(0)>=.7*sample_size)
##            self.passed
            transformed=pca.fit_transform(rescaled[:,~isnan], weights=self.weights[:,~isnan])
    ##    if partial_fit==False:





        if selection_method!=None:
            lin_fit=ExplainPCs(transformed, AlleleRescale(data),3, method=selection_method)
            for i in range(3):
    ##            pyplot.plot(pos,lin_fit[i].coef_.transpose(), alpha=.9)

                if plot==True:
                    if selection_method!='RL':
                        pyplot.plot(lin_fit[i].coef_.transpose(), alpha=.9)
                    else: pyplot.plot(lin_fit[i].scores_.transpose(), alpha=.9)

        ##            pyplot.plot(pos, lin_fit[i].scores_)
                    pyplot.show()
        ##            pyplot.plot(numpy.cumsum( sorted(list(abs( lin_fit.coef_[i])) )))
                    pyplot.show()

        else:
            lin_fit=pca.components_


        self.components, self.PCA, self.contributions, self.rescaled= transformed, pca, lin_fit, AlleleRescale(data)

    def PlotScree(self):
        assert hasattr(self, 'PCA'), "The PerformPCA() method must be run before plotting PCA summaries."

        pyplot.plot(range(1 , len( self.PCA.explained_variance_)+1),self.PCA.explained_variance_ratio_, lw=2, c='forestgreen', alpha=.9)

        pyplot.ylabel('Variance Explained')
        pyplot.xlabel('Component')

##        pyplot.show()
    def PlotPCA(self, comp1=1, comp2=2):
        assert hasattr(self, 'PCA'), "The PerformPCA() method must be run before plotting PCA summaries."
        comp1-=1
        comp2-=1
        transformed=self.components
        pca=self.PCA
        pyplot.scatter(transformed[:,comp1], transformed[:,comp2], c=colors)
        pyplot.gca().set_aspect('equal')
##        pyplot.show()
##        trans_2=pca.fit_transform(AlleleRescale( data).transpose())
##        pyplot.scatter(pca.components_[0,:], pca.c//omponents_[1,:], c=colors, s=80, alpha=.7)
        pyplot.xlabel('PC{0} ({1}%)'.format(comp1+1, numpy.round( pca.explained_variance_ratio_[comp1]*100, 3)))
        pyplot.ylabel('PC{0} ({1}%)'.format(comp2+1, numpy.round( pca.explained_variance_ratio_[comp2]*100, 3)))
        pyplot.gca().set_aspect('equal')
##        pyplot.show()

    def PlotPCA3D(self, comp1=1, comp2=2, comp3=3):
        assert hasattr(self, 'PCA'), "The PerformPCA() method must be run before plotting PCA summaries."
        comp1-=1
        comp2-=1
        comp3-=1
        transformed=self.components
        pca=self.PCA
        fig=pyplot.figure()
        ax=mplot3d.axes3d.Axes3D(fig)
        ax.scatter(transformed[:,0], transformed[:,1],transformed[:,2], c=colors, s=60)

        ax.set_xlabel('PC{0} ({1}%)'.format(comp1+1, numpy.round( pca.explained_variance_ratio_[comp1]*100, 3)))
        ax.set_ylabel('PC{0} ({1}%)'.format(comp2+1, numpy.round( pca.explained_variance_ratio_[comp2]*100, 3)))
        ax.set_zlabel('PC{0} ({1}%)'.format(comp3+1, numpy.round( pca.explained_variance_ratio_[comp3]*100, 3)))
##        pyplot.show()

    def PlotContributions(self, component=1):
        component-=1
        pyplot.plot(self.contributions [component], alpha=.9)
##        pyplot.show()
        ##            pyplot.plot(numpy.cumsum( sorted(list(abs( lin_fit.coef_[i])) )))
    def SelectFeatures(self):
        """Not implemented."""

        pass
    def Biplot(self, pc1, pc2, head_mod=2, lw=1,alpha=.2, ax=None):
        pca=self.PCA
        lin_fit=self.contributions
        rel_ind=((lin_fit[0]!=0)+ (lin_fit[1]!=0) +(lin_fit[2]!=0))>0
        indices=numpy.where(rel_ind==True)[0]
        ax=pyplot.gca()
        seaborn.set_style('white')
        ax.set_aspect('equal')
    ##    ax.grid(False)

        ax.scatter(pca.components_[pc1,:],pca.components_[pc2,:], c=colors, s=20, alpha=.9, zorder=2 )
        scale_1=numpy.max(abs(pca.components_[pc1,:]))/(numpy.max(abs( lin_fit[pc1])))
        scale_2=numpy.max(abs(pca.components_[pc2,:]))/(numpy.max(abs( lin_fit[pc2])))
        scale=min(scale_1,scale_2)
        for ind in indices:
            ax.plot([0,-1*lin_fit[pc1][ind]*scale],[0,-1* lin_fit[pc2][ind]*scale], c='grey', alpha=.01, zorder=3)
            ax.arrow(0,0,-1*lin_fit[pc1][ind]*scale,-1* lin_fit[pc2][ind]*scale,lw=lw, head_width=0.005*head_mod, head_length=0.01*head_mod,color='black', alpha=alpha, zorder=3)


        pyplot.xticks(size=14)
        pyplot.yticks(size=14)
        pyplot.xlabel('PC{0} ({1:.1f}%)'.format (pc1+1, numpy.round( pca.explained_variance_ratio_[pc1]*100, 3)), size=16)
        pyplot.ylabel('PC{0} ({1:.1f}%)'.format (pc2+1, numpy.round( pca.explained_variance_ratio_[pc2]*100, 3)), size=16)


    def BuildSequenceVariantTable(self,outfile, seq):
        rel_ind=numpy.where( self.allele_frequency>0)[0]
        #local variables
    ##    ind=numpy.array(good_ind)
        nt_dict={'A':1, 'T':2, 'C':3,'G':4}
        nt_classes={'A':'PUR', 'T':'PYR', 'C':'PYR','G':'PUR'}
        mut_classes={True:'Transition', False:'Transversion'}
        inv_nt_dict={1:'A', 2:'T', 3:'C',4:'G'}
        array=numpy.copy( self.pileup[:,:,:])
        array/=array.sum(1)[:,None,:]

        workbook=xlsxwriter.Workbook(outfile)
        worksheet=workbook.add_worksheet()

        int_format=workbook.add_format()
        int_format.set_num_format('#,##0')

        float_format=workbook.add_format()
        float_format.set_num_format('0.00E+00')

    ##    rel_ind=((lin_fit[0].coef_!=0)+ (lin_fit[1].coef_!=0) +(lin_fit[2].coef_!=0))>0
        lin_indices=numpy.where(rel_ind==True)[0]
        relevant_positions=self.passed_positions[rel_ind]
        mean_major=numpy.mean(self.major_allele_proportion[:,rel_ind],0)
        std_major=numpy.std(self.major_allele_proportion[:,rel_ind],0)
        min_major=numpy.min(self.major_allele_proportion[:,rel_ind],0)
        var_freq=self.allele_frequency[rel_ind]
    ##    fst=ComputeFst(table)[relevant_positions]
        row=2
    ##    worksheet.write(0,0, 'Sequence Information')
    ##    worksheet.write(0,3, 'Regression Coefficients')
        header=['Position',	'Consensus',	'Variant','Type' ,	'Mean', 'Std', 'Min', 'Freq', 'Fst']
        for column in range(6):
            worksheet.write(1, column, header[column])
        for index in range(len( relevant_positions)):
            column=0
            pos=relevant_positions[index]
            cons_allele=seq[pos-1]      #Pile-up is 1-coordinate, consensus positions are 0-coordinate. CXOnvert to zero.
            cons_index=nt_dict[cons_allele]
            mean_freq=numpy.mean(array[:,:,pos],0)
            mean_freq[0]=0.
            mean_freq[cons_index]=0.
            try:
                major_var_ind=numpy.argmax(mean_freq)
                var_allele=inv_nt_dict[major_var_ind]
            except:
                print mean_freq
                print major_var_ind
                print array[:,:,pos]
                continue
            mut_type=mut_classes[nt_classes[cons_allele]==nt_classes[var_allele]]
            entries=[pos, cons_allele, var_allele, mut_type, mean_major[index], std_major[index],min_major[index],var_freq[index]]#, fst[index] ]#,  lin_fit[0].coef_[lin_indices[index]], lin_fit[1].coef_[lin_indices[index]], lin_fit[2].coef_[lin_indices[index]]]
            for item in entries:
                if column==0: cell_format=int_format
                else: cell_format=float_format
                worksheet.write(row, column,item, cell_format )
                column+=1
            row+=1
        merge_format = workbook.add_format({
        'bold': 1,
        'border': 1,
        'align': 'center',
        'valign': 'vcenter'})
        worksheet.merge_range('A1:D1', 'Sequence Information', merge_format)
        worksheet.merge_range('E1:G1', 'Regression Coefficients', merge_format)
        workbook.close()
##class ConTExtLine():
##    strand={'0':'+', '16':'-', '4':'*', 'strand': '*'}
##    def __init__(self, row):
##
##        self.Seq1=row[0]
##        self.Strand1=ConTExtLine.strand[row[1]]
##        self.Start1=int(row[2])
##        self.End1=int(row[3])
##        self.Seq2=row[6]
##        self.Strand2=ConTExtLine.strand[row[7]]
##        self.Start2=int(row[8])
##        self.End2=int(row[9])
##        self.Cigar1=row[12]
##        self.Cigar2=row[13]
##        self.MapQ1=int(row[14])
##        self.MapQ2=int(row[15])
##        #self.Mismatches1=row[16]
##        #self.Mismatches2=row[17]

class Datum(object):
    colorin = colorConverter.to_rgba('red')
    colorout = colorConverter.to_rgba('None')
    def __init__(self, x, y, include=False):
        self.x = x
        self.y = y
##        if include: self.color = self.colorin
##        else: self.color = self.colorout
        self.color=self.colorout

class LassoManager(object):
    def __init__(self, ax, data, colors, sample_list):
        self.axes = [ax]
        self.canvas = ax.figure.canvas
        self.data = [data]


        self.sample_list=numpy.array(sample_list)
        self.x_list=numpy.array([d.x for d in data])
        self.y_list=numpy.array([d.y for d in data])
##        self.sample_set=

        self.Nxy = [ len(data) ]

        facecolors = [d.color for d in data]
        self.xys = [ [(d.x, d.y) for d in data] ]
        fig = ax.figure
        self.collection = [ RegularPolyCollection(
            fig.dpi, 6, sizes=(60,),
            facecolors=colors,
            edgecolors=facecolors,
            offsets = self.xys[0],
            transOffset = ax.transData,linewidths=(1,), alpha=.7)]

        ax.add_collection(self.collection[0])

        self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)

    def callback(self, verts):
##        self.axes[1].cla()
        axind = self.axes.index(self.current_axes)
        facecolors = self.collection[axind].get_edgecolors()
        p = path.Path(verts)

        ind = p.contains_points(self.xys[axind])
##        print ind
##        print self.counts[ind]
##        self.axes[1].scatter(range(len(self.CN_map[ind])),sorted( self.CN_map[ind]))
##        print self.xys[axind]
##        print ind.shape, self.x_list.shape
        print numpy.array( self.sample_list) [ind]


        for i in range(len(self.xys[axind])):
            if ind[i]:
                facecolors[i] = Datum.colorin
            else:
                facecolors[i] = Datum.colorout
        self.canvas.draw_idle()
##        self.canvas.show()
##        self.canvas.draw()
        self.canvas.widgetlock.release(self.lasso)
##        self.axes[1].scatter(ind,ind, c='g')
        del self.lasso

    def onpress(self, event):
        if self.canvas.widgetlock.locked(): return
        if event.inaxes is None: return
        self.current_axes = event.inaxes

        self.lasso = Lasso(event.inaxes, (event.xdata, event.ydata), self.callback)
        # acquire a lock on the widget drawing
        self.canvas.widgetlock(self.lasso)

    def add_axis(self, ax,  data):
        self.axes.append(ax)
        self.data.append(data)

        self.Nxy.append( len(data) )

        facecolors = [d.color for d in data]
        self.xys.append( [(d.x, d.y) for d in data] )
        fig = ax.figure
##        self.collection.append( RegularPolyCollection(
##            fig.dpi, 6, sizes=(100,),
##            facecolors=facecolors,
##            offsets = self.xys[-1],
##            transOffset = ax.transData))
##
##        ax.add_collection(self.collection[-1])

def UpdateTable(table):
    return numpy.vstack((samples, table))

#Functions for a

def COlorBlindColors(colors):
    colors=numpy.copy(colors)
    colors[numpy.where(colors=='firebrick')]='#6f0000'
    colors[numpy.where(colors=='forestgreen')]=matplotlib.colors.rgb2hex((0, 146./256,146/256.))
    colors[numpy.where(colors=='blue')]=matplotlib.colors.rgb2hex((109/255., 182./255,255/255.))
    colors[numpy.where(colors=='purple')]=matplotlib.colors.rgb2hex((73./255., 0.,146/255.))
##    colors[numpy.where(colors=='orange')]=matplotlib.colors.rgb2hex((255./255., 85/255.,22/255.))
    return colors

def ComputeMajorAlleleStatistics(snp_array,FDR=.01, error_rate=1e-3, depth_cutoff=10):
    """
    Computes two statistics from the SNP pileup.
    First, it defines the major allele at each positions as the allele present at
    the highest average allele proportion across the populations. It

    Second, it estimates the proportion of samples that contain a variant allele at a
    position in at least one repeat unit. Variant alleles are called by rejecting
    the null hypothesis that all observed variants are sequencing errors using a
    conservative error rate epsilon (the minimum quality score). This does not account for
    errors introduced by PCR, which may be important. For each position i in
    sample j total read count X_ij and variant read count k_ij are determined. We
    define the variant allele as the allele with the second greatest average allelw
    proportion. The null model is:

    k_ij ~ Binom(X_ij, epsilon)

    A p-value is computed for each sample and the null hypothesis is (or isn't)
    rejected using the Benjami-Hochsberg multiple testing correction for the
    number of samples at the given FDR (default .01). Note, this null model does
    not account for PCR-duplicates, and the assumption of a simple sampling
    distribution is usually simplistic.

    Parameters:
        snp_array: the I x J x K 3-dimensional snp pileup array

        FDR: Accepted false-discovery rate for a position.

        error_rate: The expected rate of sequencing errors. For base-calling errors,
        the default is 1e-3, which corresponds to lowest Phred quality score consider
        when ConTExt builds SNP pileups (30).

        depth_cutoff: A statistic will only be computed for a consensus position
        if the read depth at the position is higher than the cutoff in every sample.

    Returns:
        major_allele_prop: An I x K array indicating the major allele proportion
        at positon k in sample i

        var_allele_freq: A length K array the fraction of samples containing a
        variant allele in at least one repeat unit

        pos: The indices of consensus allele positions which survive the cutoff


    """
    samples=snp_array.shape[0]
    read_count=snp_array[:,:,:] .sum(1)
##    print read_count.shape


    major_allele_prop,major_allele, all_pos=GetMajorAlleleFreq(snp_array, cutoff=depth_cutoff)
    var_allele_prop=GetOtherAlleleProportion(snp_array, depth_cutoff, order=2)
    read_count=read_count[:,all_pos]
    var_counts=(read_count*var_allele_prop).astype(int)
    read_count=read_count.astype(int)
    min_cvg= numpy.min(read_count,0)
    allele_freq=[]
    var_pres=[]

    for i in range(read_count.shape[1]):
        p_vals=(1.- scipy.stats.binom.cdf(var_counts[:,i], read_count[:,i], error_rate))
##        p_vals=(1.- scipy.stats.binom.cdf(snp_array[:,:,i], read_count[:,i][:,None], error_rate))
##        print p_vals.shape
        rej,q_vals,sidak, bonf=statsmodels.stats.multitest.multipletests(p_vals,FDR, method='fdr_bh')
        allele_freq.append(float( rej.sum()) /samples)
        var_pres.append(rej)
    return major_allele_prop, numpy.array( allele_freq) ,numpy.array( var_pres).transpose(),  all_pos, major_allele

def EmptyFunction(data, ax=0):
    return data
def DriveFxn( d):
    gen_dist= numpy.logspace(0, 2,100)
    print gen_dist
    gen_dist=(1-numpy.exp(-2*gen_dist/100))/2.
##    pyplot.plot(gen_dist)
##    pyplot.show()
    return .5 + .5*(d*(1-gen_dist)+(1-d)*gen_dist)

def PlotSummaries(samples, labels,fxn, log=False,scatter=True, mask=[], mask_value=None,plot='box', width=.1, ylim=None, figsize=None):
    fig,ax=pyplot.subplots()
    if figsize!=None:
        fig.set_size_inches(figsize)
    if ylim!=None:
        pyplot.ylim(ylim)
    pyplot.tick_params('y', length=6)


    variances=[]
    for i in range(len(samples)):
##        logvar=numpy.log10( numpy.mean (samples[i],0))
        logvar=fxn(samples [i],0)
        if mask!=[] and mask_value!=None:
            logvar[numpy.array(mask[i])]=mask_value
        if mask!=[] and mask_value==None:
            logvar= logvar[numpy.array(mask[i])]
        if log==True:
            logvar=numpy.log10 (logvar)
        nan_pos=numpy.isinf(logvar)
        variances.append(logvar[~nan_pos] )
##    pyplot.violinplot( variances)
    if plot=='violin':
        sns=seaborn.violinplot(data= variances,  orient='v')
    if plot=='box':
        sns=seaborn.boxplot(data= variances,  orient='v',fliersize=0)
##        pyplot.boxplot(variances)
    if scatter==True:
        for i in range(len(variances)):
            print numpy.nanmean(variances[i])
            jitter=numpy.random.uniform(-1*width,width,size=len(variances[i]))
            pyplot.scatter(jitter+i, variances[i],s=15, c='orange', alpha=.6)
    pyplot.xticks(numpy.arange(0,len(labels)), labels, size=16)
    pyplot.yticks( size=16)
    boxes=sns.artists
    for i,box in enumerate(boxes):
        box.set_facecolor('dodgerblue')

##    pyplot.ylabel('Standard Deviation', size=14)
    return variances

def PlotAllListedPositions( outDir, table,positions, condition_array,pop=colors):

    MakeDir(outDir)

    for pos in numpy.where( condition_array==True)[0] :
##        nt=seq[pos-1]
        plot=PlotMajorAlleleFreqByPopulation( table, positions, pos, 10.**-3.5, pop=pop)
        pyplot.savefig("{0}/{1}.jpeg".format(outDir, positions[ pos]))
        pyplot.close()

def BuildMajorConsensus(pos, major_ind, seq):
    consensus=numpy.array( ['']*len(seq), '|S1')
    nt_dict={'A':1, 'T':2, 'C':3,'G':4}
    inv_nt_dict={0:'Cons', 1:'A', 2:'T', 3:'C',4:'G'}
    test=[]
    for p in pos:
        seq_position=p-1
##        print p

        ind=major_ind[p]

        if ind==0:
##            print
            consensus[seq_position]=seq[seq_position]
        else:
            print p, seq[seq_position], inv_nt_dict[ind].upper()
            test.append( seq[seq_position]==inv_nt_dict[ind].upper())
            consensus[seq_position]=inv_nt_dict[ind].upper()
    print numpy.mean(test)

    return ''.join( consensus)

def PlotMajorAlleleFreqByPopulation(table, pos_list, pos1, offset, pop=colors):
    color_dict={'firebrick':.5, 'forestgreen':1.5, 'blue':2.5, "orange":3.5,'purple':4.5}
    array=numpy.copy( table[:, pos1])


##    freq1=table[ind, nt_dict[var1], pos1]/table[ind, :, pos1].sum(1)
    freq1=array
    print freq1.shape

##    print numpy.log10( freq1+offset)
    color=numpy.array( pop)
    pop_dict=[]
    color_set=sorted(list(set(color)))
    pop_freq=[]
    for c in color_set:
        ind=(color==c)
        pop_freq.append(numpy.log10( freq1+offset))
##    seaborn.violinplot(data=pop_freq)
    jitter=numpy.random.uniform(-.2, .2, size=len(ind))
    x_pos=numpy.array( [color_dict[c] for c in color])
##    print numpy.log10( freq1+offset)
##    print x_pos
##    pyplot.scatter(x_pos+jitter+.5, numpy.log10(1.- freq1+offset),c=color, s=80, alpha=.8)

    pyplot.scatter(x_pos+jitter+.5,  freq1,c=color, s=80, alpha=.8)

##    pyplot.ylabel('Cons Freq at Pos {0} + 1e-3.5 (Log10)'.format(pos1), size=20)
    pyplot.ylabel('Cons Freq at Pos {0}'.format(pos_list[ pos1]), size=20)
    pyplot.xticks(numpy.arange(1,6),['B','I','N','T','Z'], size=20)
    pyplot.yticks(size=20)
    pyplot.xlim(.5, 5.5)
##    pyplot.ylim(-3.6, .1)
    pyplot.ylim(-.05, 1.05)
##    pyplot.title(var1)
##    else: color='firebrick'
##    pyplot.show()



def PlotBiplots(pca,linfit,alpha=.05, colors=colors):
    ax=pyplot.subplot(121)
    Biplot(pca, linfit, 0,1,alpha=alpha ,colors=colors)
    pyplot.subplot(122)
    Biplot(pca, linfit, 1,2,alpha=alpha, colors=colors)
    pyplot.show()


def PlotSummaryPanel(table_list,summ_list, mask_list=None, mask_vals=None):
##    fig=pyplot.figure(figsize=(6.75,4.5 ))
    gs=gridspec.GridSpec(3, 1)
    ax_list=[]
    ax_list.append( pyplot.subplot(gs[0]))
    PlotSummaries([s[0] for s in summ_list], repeat_labels,mask=[s[1]==0. for s in summ_list],fxn=numpy.mean, mask_value=1 )
    return fig, gs, ax_list
def PlotKey(colors=COlorBlindColors(colors)):
    legend_key={'Bei':'#6f0000', 'Ith':'#009191', 'Net':'#6db6ff', 'Tas':'orange', 'Zim':'#490092'}
    fig,ax=pyplot.subplots()
    fig.set_size_inches((5,2))
    leg_obj=[]
    leg_name=[]
    for key in sorted( legend_key.keys()):
        leg_obj.append(pyplot.Line2D((0,1),(0,0), c=legend_key[key], marker='o', linestyle=''))
        leg_name.append(key)
    leg=pyplot.legend(leg_obj, leg_name, ncol=5,bbox_to_anchor=(1.1, 1.05), fontsize=10, frameon=False)
    leg.draggable()


def PlotBiplotPanel(pcas, shape, titles=[],pc=(0,1), colors=COlorBlindColors(colors)):
##    pyplot.gca().set_aspect('equal')
    legend_key={'Bei':'firebrick', 'Ith':'forestgreen', 'Net':'blue', 'Tas':'orange', 'Zim':'purple'}
    n_pcas=len(pcas)
##    pyplot.figure(figsize=(6,4))
    pyplot.close()
    fig=pyplot.figure( dpi=100)
    fig.set_size_inches((3.5,9))
    gs=gridspec.GridSpec(shape[0], shape[1], )#width_ratios=[1]*n_pcas, height_ratios=[1]*n_pcas)

    ax_list=[]
    i=0
    for x in range(shape[0]):
        for y in range(shape[1]):

            print i
            ax_list.append(fig.add_subplot(gs[x,y]))
            ax_list[-1].set_aspect('equal', adjustable='datalim')
##            ax_list[-1]. locator_params(axis='x',nticks=5)
            ax_list[-1].xaxis.major.locator.set_params(nticks=3)
            ax_list[-1].yaxis.major.locator.set_params(nticks=3)
            ax_list[-1].tick_params(axis=u'both', which=u'both',length=2)
##            ax_list.append(pyplot.subplot2grid(shape, (x,y)))
    ##        print pcas[i]

            Biplot(pcas[i][0],pcas[i][1], pc1=pc[0], pc2=pc[1],ax=ax_list[-1])
            if titles!=[]:
                pyplot.title(titles[i], size=20)
            i+=1
    leg_obj=[]
    leg_name=[]
    for key in sorted(legend_key.keys()):
        leg_obj.append(pyplot.Line2D((0,1),(0,0), c=legend_key[key], marker='o', linestyle=''))
        leg_name.append(key)
    leg=pyplot.legend(leg_obj, leg_name, ncol=3,bbox_to_anchor=(1.1, 1.05), fontsize=14, frameon=True)


    return  fig,gs, ax_list, leg

##    pyplot.tight_layout()


def ComputeFst(snp_array,ax =0 , pop=colors):
    """Computes Fst on allele proportion at each locus following Thomas Nagylaki "Fixation Indices
    in subdivided populations.

    Parameters:

    snp_array: the I x J x K 3-dimensional snp pileup array
    pop: An array indicating sample labels. For the GDL, this stored as a global variable which is used by default.
    """



    freq_array=snp_array/snp_array.sum(1)[:,None,:]
##    print freq_array
    #Eq 3 Total population mean
    avg_prop=numpy.nanmean( freq_array,0)

    #Eq 24 .a and .b
    div_array=avg_prop**2
    h_t=1-div_array.sum(0)
##    print h_t[1519]


    #Eq 2 : Subpopulation means
    pop=numpy.array(pop)
    pop_set=sorted( list(set( pop) ))
    pop_means=[]
    pop_count=[]
    for p in pop_set:
        pop_means.append( numpy.nanmean( freq_array[pop==p] ,0))
        pop_count.append(numpy.isnan(freq_array[pop==p].sum(1)).sum(0))

    pop_means=numpy.array(pop_means)
    pop_count=numpy.array(pop_count)
    #Eq 10.a
    var_i=numpy.var(pop_means,0 )
    print var_i
    exp_var_i=avg_prop*(1-avg_prop)
##    print exp_var_i[:,1519]
##    print var_i[:,1519]
    fst_i=var_i/(exp_var_i)
##    print fst_i[:,1519]
    #Eq 32.a
    fst=numpy.nansum( fst_i*exp_var_i,0)/h_t
    print var_i.shape
    print var_i[:,fst>1]
    print avg_prop[:,fst>1]
    print h_t[fst>1]
##    print freq_array[:,:,fst>1]
    print pop_count.shape

##    print fst[1519]
    fst[h_t==0]=0.
    fst[(pop_count>=5).sum(0)>0]=numpy.nan
    return fst


def AnovaOnPCA(PCA):
    pop_samples=[]
    for p in set(colors):
        pop_samples.append((PCA[colors==p, :3]))
    return scipy.stats.f_oneway(*pop_samples)

def PlotPCA(data,pos=[], selection_method=None, colors=COlorBlindColors( colors), plot=True, partial_fit=False):
##    data=data[numpy.array( good_ind),:]
    if samples!='':
        print 'b'
        keep_ind=numpy.array( [bad_strains.count(s)==0 for s in samples])
        colors=numpy.array(colors)
        data=data
        nan_ind=numpy.isnan(data.sum(1))
##        data=data[~nan_ind,:]
##        pca=sklearn.decomposition.IncrementalPCA(92-len(bad_strains),whiten=True, copy=False,batch_size= 100)

        pca=sklearn.decomposition.PCA(data.shape[0],whiten=True, copy=False)
    else:

        pca=sklearn.decomposition.PCA(92,whiten=True, copy=False)
##    if partial_fit==False:
    transformed=pca.fit_transform(AlleleRescale( data[keep_ind]))
    if plot==True:
        pyplot.scatter(transformed[:,0], transformed[:,1], c=colors)
        pyplot.gca().set_aspect('equal')
##        pyplot.show()
##        trans_2=pca.fit_transform(AlleleRescale( data).transpose())
##        pyplot.scatter(pca.components_[0,:], pca.c//omponents_[1,:], c=colors, s=80, alpha=.7)
        pyplot.xlabel('PC1 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[0]*100, 3)))
        pyplot.ylabel('PC2 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[1]*100, 3)))
        pyplot.gca().set_aspect('equal')
        pyplot.show()
##        pyplot.scatter(trans_1[:,0], pca.components_[0,:], c=colors)
##        pyplot.show()

    ##    return
        pyplot.scatter(transformed[:,1], transformed[:,2], c=colors)
##        pyplot.scatter(pca.components_[1,:], pca.components_[2,:], c=colors, s=80, alpha=.7)
        pyplot.xlabel('PC2 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[1]*100, 3)))
        pyplot.ylabel('PC3 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[2]*100, 3)))
        pyplot.gca().set_aspect('equal')
        pyplot.show()

    ##    pyplot.scatter(pca.components_[2,:], pca.components_[3,:], c=colors, s=80, alpha=.7)
    ##    pyplot.xlabel('PC1 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[0]*100, 3)))
    ##    pyplot.ylabel('PC2 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[1]*100, 3)))
    ##    pyplot.show()
    ##
    ##    pyplot.scatter(pca.components_[5,:], pca.components_[6,:], c=colors, s=80, alpha=.7)
    ##    pyplot.xlabel('PC2 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[1]*100, 3)))
    ##    pyplot.ylabel('PC3 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[2]*100, 3)))
    ##    pyplot.show()

        fig=pyplot.figure()
        ax=mplot3d.axes3d.Axes3D(fig)
        ax.scatter(transformed[:,0], transformed[:,1],transformed[:,2], c=colors, s=60)
        ax.set_xlabel('PC1 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[0]*100, 3)))
        ax.set_ylabel('PC2 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[1]*100, 3)))
        ax.set_zlabel('PC3 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[2]*100, 3)))
        pyplot.show()

        pyplot.plot(range(1 , len( pca.explained_variance_)+1),pca.explained_variance_ratio_, lw=2, c='forestgreen', alpha=.9)
        pyplot.ylabel('Variance Explained')
        pyplot.xlabel('Component')
        pyplot.show()

    if selection_method!=None:
        lin_fit=ExplainPCs(transformed, AlleleRescale(data),3, method=selection_method)
        for i in range(3):
##            pyplot.plot(pos,lin_fit[i].coef_.transpose(), alpha=.9)

            if plot==True:
                if selection_method!='RL':
                    pyplot.plot(lin_fit[i].coef_.transpose(), alpha=.9)
                else: pyplot.plot(lin_fit[i].scores_.transpose(), alpha=.9)

    ##            pyplot.plot(pos, lin_fit[i].scores_)
                pyplot.show()
    ##            pyplot.plot(numpy.cumsum( sorted(list(abs( lin_fit.coef_[i])) )))
                pyplot.show()
    else:
        lin_fit=pca.components_
        return transformed, pca, lin_fit, AlleleRescale(data)

    return transformed, pca, lin_fit, AlleleRescale(data)


def ExplainPCs(pca, data, n=3, method='EN'):
    """Methods for solving the linear model:

        EN: Elastic Net -- Accomplishes feature selection and captures correlated variables.
        This is the preferred method for selecting the positions driving a PCA, but can be very slow (may take a few minutes)

        RR: Ridge Regression -- Does not select features, but does identify how stronlgy each
        feature drives the PCA. Very fast.

        RL: Randomized Lasso -- Accomplishes feature select. Not particularly recommended. Elastic Net is  better.

        3-fold cross-validation determines the Elastic Net alpha and lambda parameters,
        and LOO CV determines Ridge Regression regularizartion parameter.
        """

    fits=[]
    for i in range(n):
        print data.shape
        print pca.transpose().shape
        if method=='EN':
            lin_fit=sklearn.linear_model.ElasticNetCV([.01,.1, .5, .7, .9, .95, .99], normalize=False,  fit_intercept=False, copy_X=False)
            lin_fit.fit(data, pca[:,i]-numpy.mean(pca[:,i]))
        if method=='RR':
            lin_fit=sklearn.linear_model.RidgeCV(10.**numpy.linspace(-5,3), normalize=False,  fit_intercept=False)
            lin_fit.fit(data, pca[:,i]-numpy.mean(pca[:,i]))
##        cv=sklearn.model_selection.GridSearchCV(sklearn.linear_model.ElasticNet( normalize=False, copy_X=False),param_grid={'alpha':10.**(numpy.arange(-20,0,1.)), 'l1_ratio':[.001, .01,.1, .5, .7, .9, .95, .99]})
##        cv.fit(data, pca.components_[i,:])
##        lin_fit=cv.best_estimator_
##        lin_fit=sklearn.linear_model.ElasticNetCV(.5,n_alphas=200, normalize=False, copy_X=False,cv=2)
##        lin_fit=sklearn.linear_model.LassoCV( normalize=False, copy_X=False)

##
        elif method=='RL':
            lin_fit=sklearn.linear_model.RandomizedLasso(fit_intercept=False)
            lin_fit.fit(data, pca[:,i]-numpy.mean(pca[:,i]))
##            print lin_fit.coef_.shape
##        cv_fit.fit(data, pca.components_[i,:])
        print data.shape
        print pca.transpose().shape
    ##    lin_fit.fit(data, pca.components_.transpose())
##        alphas=np.linspace(cv_fit.alphas_[0], .1 * cv_fit.alphas_[0], 6)
##        lin_fit=sklearn.linear_model.RandomizedLasso(alpha=alphas)

##        print lin_fit.coef_
##        print lin_fit.alpha_
##        print lin_fit.score()
        fits.append(lin_fit)
        if method!='RL':
            predict=lin_fit.score(data, pca[:,i]-numpy.mean(pca[:,i]))
    ##        predict=lin_fit.predict(data, pca.components_[i,:])
            print predict
            if method=='EN':           print lin_fit.l1_ratio_

            print lin_fit.alpha_
    ##        print lin_fit.alphas_
    ##        pyplot.scatter(pca.components_[0,:], predict[:,0])
##    print ((pca.components_[0,:]- predict[:,0])**2).sum()
##    pyplot.show()

    return fits


def PlotAlleleProportions(hist_freq, hist_pos, colors=colors, highlight=[] ):

    """Plot"""

    color_dict={'B':'firebrick', 'I':'forestgreen', "N":'blue', 'T':"orange",'Z':'purple'}

    offset=-.2

    rescale=hist_freq[:,:]
    for c in list(set(colors)):
##        f1=[ pyplot.scatter(hist_pos[rel_ind]+offset, hist_freq[i,rel_ind], c=c, alpha=.7) for i in range(len(colors)) if colors[i]==c]
        outlier=highlight!=0
        if outlier!=[] and highlight!=[]:
            f1=[ pyplot.scatter(hist_pos[:][outlier] +offset, rescale [i,:][outlier], c=c, edgecolors='dodgerblue',linewidths=1, s=80,  alpha=.7) for i in range(len(colors)) if colors[i]==c]
        f1=[ pyplot.scatter(hist_pos[:]+offset, rescale [i,:], c=c, alpha=.7) for i in range(len(colors)) if colors[i]==c]
        offset+=.1
    pyplot.ylabel('Consensus Allele Frequency')
    pyplot.xlabel('Consensus Position')
##    cutoff
##    for  i in range(len(lin_fit)):
    return


def TestPCA(pca, data,lin_fit, cutoffs):
    r=[]
    for c in cutoffs:
        ind=((lin_fit[0].coef_>=c)+(lin_fit[1].coef_>=c)+(lin_fit[2].coef_>=c)>0)
        pca_test=sklearn.decomposition.PCA(3)
        pca_test.fit(data[:, ind].transpose())
        rss=sum(sum((pca.components_[:3,:]-pca_test.components_[:3,:])**2))**.5
        print sum(ind)
        print rss
        r.append(rss)
    return(r)


def PlotSelectablePCA(pca, sample ):

    colors={'B':'firebrick', 'I':'forestgreen', "N":'blue', 'T':"orange",'Z':'purple'}
    heads=[s[0] for s in sample]
    color_list=[colors[h] for h in heads]
##    subplot_kw = dict( autoscale_on=True)
##    fig, ax = plt.subplots(subplot_kw=subplot_kw)
    scatter=zip(pca.components_[0,:], pca.components_[1,:],)
    ax = plt.subplot(1,1,1)
    x,y=zip(*scatter)
    data=[ Datum(x[i],y[i]) for i in range(len(x))]
    lman = LassoManager(ax, data,color_list, sample )
##    pyplot.ylabel('{0}, {1} strand'.format( seq2, quadrant[1]))
##    pyplot.xlabel('{0}, {1} strand'.format( name, quadrant[0]))
    ax.autoscale(True)
    pyplot.xlabel('PC1 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[0]*100, 3)))
    pyplot.ylabel('PC2 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[1]*100, 3)))
    pyplot.show()
##    lman.add_axis(ax2, data2)


##    pyplot.draw()
##    lman.outhandle.close()
    return lman, color_list

def Isomap(data, param1,param2):
    fitted=[]
    score=[]
    matrix=numpy.ndarray((max(param1)+1, max(param2)+1))
    matrix.fill(0.)
    for i in param1:
        for j in param2:
            iso=sklearn.manifold.Isomap(i,j)
            fitted.append( iso.fit(data))
            score.append(iso.reconstruction_error() )

            matrix[i,j]=score[-1]
    return fitted, score, matrix
def GetKeepInd(samples):
    return numpy.array( [bad_strains.count(s)==0 for s in samples])


def DummyData(colors, snps, group_1=['blue', 'orange']):
    freq=numpy.ndarray((len(colors), 500))
    snp_pos=numpy.random.choice(range(500), size=snps, replace=False)
    snp_freq=numpy.array([.95]*500)
    var_freq=numpy.array([.95]*500)
    var_freq[snp_pos]=numpy.random.uniform(.2,.8, size=snps)
##    return snp_freq, var_freq
    for i in range(len(colors)):
        if group_1.count( colors[i])==0:
            freq[i,:]=scipy.stats.binom.rvs(200, snp_freq)/200.
        else:
            freq[i,:]=scipy.stats.binom.rvs(200, var_freq)/200.

    return freq, snp_pos

def DummyDataByPop(colors, snps, group_1=['blue', 'orange']):
    freq=numpy.ndarray((len(colors), 500))
    pop_freq={}
    pop_pos={}
    pops= sorted(list(set(colors)))
    for p in range(len(pops)):
        snp_pos=numpy.random.choice(range(500), size=snps, replace=False)
        snp_freq=numpy.array([.95]*500)
        var_freq=numpy.array([.95]*500)
        var_freq[snp_pos]=numpy.random.uniform(.2,.8, size=snps)
        pop_pos[pops[p]]= snp_pos
        pop_freq[pops[p]]= var_freq
##    return snp_freq, var_freq
    for i in range(len(colors)):

        freq[i,:]=scipy.stats.binom.rvs(100, pop_freq[colors[i]] )/100.


    return freq, pop_pos

def DummyDataByPopWithMix(colors, snps, group_1=['blue', 'orange']):
    freq=numpy.ndarray((len(colors), 500))
    pop_freq={}
    pop_pos={}
    pops= sorted(list(set(colors)))
    for p in range(len(pops)-1):
        snp_pos=numpy.random.choice(range(500), size=snps, replace=False)
        snp_freq=numpy.random.uniform(.97,.99, size=500)
        var_freq=numpy.copy(snp_freq)
        var_freq[snp_pos]=numpy.random.uniform(.2,.8, size=snps)
        pop_pos[pops[p]]= snp_pos
        pop_freq[pops[p]]= var_freq

##    pop_freq[pops[ -1]]=(pop_freq[pops[0]]+pop_freq[pops[1]])
##    return snp_freq, var_freq
    for i in range(len(colors)):
        if colors[i]!=pops[-1] and colors[i]!=pops[-2]:
            if numpy.random.random()>.2:
                freq[i,:]=scipy.stats.binom.rvs(100, pop_freq[colors[i]] )/100.
            else:
                freq[i,:]=scipy.stats.binom.rvs(100, pop_freq[colors[numpy.random.randint(0,3)]] )/100.
        elif colors[i]==pops[-2]:
            if numpy.random.random()>.5:
                alpha=int( numpy.random.random()*100)
                freq[i,:]=( scipy.stats.binom.rvs(100-alpha, pop_freq[pops[1]] )+scipy.stats.binom.rvs(alpha, pop_freq[pops[0]] ) )/100.
            else:
                alpha=int( numpy.random.random()*100)
                freq[i,:]=( scipy.stats.binom.rvs(100-alpha, pop_freq[pops[2]] )+scipy.stats.binom.rvs(alpha, pop_freq[pops[0]] ) )/100.

        else:
            if numpy.random.random()>.5:
                alpha=int( numpy.random.random()*100)
                freq[i,:]=( scipy.stats.binom.rvs(100-alpha, pop_freq[pops[-2]] )+scipy.stats.binom.rvs(alpha, pop_freq[pops[1]] ) )/100.
            else:
                alpha=int( numpy.random.random()*100)
                freq[i,:]=( scipy.stats.binom.rvs(100-alpha, pop_freq[pops[-2]] )+scipy.stats.binom.rvs(alpha, pop_freq[pops[2]] ) )/100.


    pop_set=list( set([s for s in pop_pos.values]))
    return freq, pop_pos, pop_set

def SelectFeatures(regressors):
    """Identify the non-zeroe coeficients in a list of sparse regression."""
    relevant_positions=[]
    for r in regressors:
        relevant_positions.append(numpy.where(r.coef_!=0.)[0])
    return numpy.hstack(relevant_positions)
def PlotDummyPositions(pos_dict):
    for color in pos_dict.keys():
        pyplot.scatter(pos_dict[color], [0.]*len(pos_dict[color]), c=color,s=30, alpha=.8)


def PlotScatterAtPositions(table, pos1, pos2, var1, var2, offset=0, log=False, pop_colors=False):
    ind=numpy.array(good_ind)
    nt_dict={'A':1, 'T':2, 'C':3,'G':4}
    freq1=table[ind, nt_dict[var1], pos1]/table[ind, :, pos1].sum(1)
    freq2=table[ind, nt_dict[var2], pos2]/table[ind, :, pos2].sum(1)
    if pop_colors==True: color=colors
    else: color='firebrick'
    if log==False:
        pyplot.scatter(freq1, freq2, c=color, s=60, alpha=.8)
        pyplot.ylabel('Variant Freq at Position {0}'.format(pos1))
        pyplot.ylabel('Variant Freq at Position {0}'.format(pos2))
        pyplot.xlim(min(  freq1)-.01, max( freq1)+.01)
        pyplot.ylim(min(freq2)-.01, max( freq2)+.01)
    elif log==True and offset!=0:
        pyplot.scatter(numpy.log10( freq1+offset), numpy.log10( freq2+offset), c=color, s=60, alpha=.8)
        pyplot.xlabel('Variant Freq at Position {0} + 1e-3.5 (Log10)'.format(pos1))
        pyplot.ylabel('Variant Freq at Position {0} + 1e-3.5 (Log10)'.format(pos2))
        pyplot.xlim(min( numpy.log10( freq1+offset))-.1, max( numpy.log10( freq1+offset))+.1)
        pyplot.ylim(min(numpy.log10( freq2+offset))-.1, max(numpy.log10( freq2+offset))+.1)
    elif log==True and offset==0:
        pyplot.scatter(numpy.log10( freq1+offset), numpy.log10( freq2+offset), c=color, s=60, alpha=.8)
        pyplot.xlabel('Variant Freq at Position {0} (Log10)'.format(pos1))
        pyplot.ylabel('Variant Freq at Position {0} (Log10)'.format(pos2))
        pyplot.xlim(min( numpy.log10( freq1+offset))-.1, max( numpy.log10( freq1+offset))+.1)
        pyplot.ylim(min(numpy.log10( freq2+offset))-.1, max(numpy.log10( freq2+offset))+.1)
    pyplot.gca().set_aspect('equal')
    pyplot.show()
    return freq1, freq2

def PlotFreqVersusCoverage(table, pos1, var1, offset=0, highlight=[], log=False, pop_colors=False):
    """Diagnostic plot for thinking about allele frequency."""
    ind=numpy.array(good_ind)
    nt_dict={'A':1, 'T':2, 'C':3,'G':4}
    freq1=table[ind, nt_dict[var1], pos1]/table[ind, :, pos1].sum(1)
    cvg=table[ind, :, pos1].sum(1)
    if pop_colors==True: color=numpy.array( colors)
    else: color='firebrick'
    if log==False:
        pyplot.scatter(cvg, freq1, c=color, s=60, alpha=.8)
        pyplot.ylabel('Variant Freq at Position {0}'.format(pos1))
        pyplot.xlabel('Coverage at Position {0}'.format(pos1))
        pyplot.xlim(min(  cvg)-.01, max( cvg)+.01)
        pyplot.ylim(min(freq1)-.01, max( freq1)+.01)
    elif log==True and offset!=0:
        pyplot.scatter(cvg, numpy.log10( freq1+offset), c=color, s=60, alpha=.8)
        pyplot.scatter(cvg[highlight], numpy.log10( freq1[highlight]+offset), c=color[highlight] , edgecolor='black', s=150, alpha=1, zorder=5)
        pyplot.xlabel('Coverage'.format(pos1))
        pyplot.ylabel('Variant Freq at Position {0} + 1e-3.5 (Log10)'.format(pos1))
        pyplot.xlim(min(cvg)-100, max( cvg)+100)
        pyplot.ylim(min(numpy.log10( freq1+offset))-.1, max(numpy.log10( freq1+offset))+.1)
    elif log==True and offset==0:
        pyplot.scatter(numpy.log10( freq1+offset), numpy.log10( freq2+offset), c=color, s=60, alpha=.8)
        pyplot.xlabel('Variant Freq at Position {0} (Log10)'.format(pos1))
        pyplot.ylabel('Variant Freq at Position {0} (Log10)'.format(pos2))
        pyplot.xlim(min( numpy.log10( freq1+offset))-.1, max( numpy.log10( freq1+offset))+.1)
        pyplot.ylim(min(numpy.log10( freq2+offset))-.1, max(numpy.log10( freq2+offset))+.1)
##    pyplot.gca().set_aspect('equal')
    pyplot.show()

def ChooseMostCorrelatedSNP(freq, pos, pca, lin_fit):
    rel_ind=numpy.where(( ((lin_fit[0].coef_!=0)+ (lin_fit[1].coef_!=0) +(lin_fit[2].coef_!=0))>0)==True)[0]
    print rel_ind
##    ind=numpy.array(good_ind)
    rescaled=AlleleRescale(freq[:,:])

    for j in range(3):
        corr_list=[]
        for i in rel_ind:
            corr_list.append(scipy.stats.pearsonr(rescaled[:,i], pca.components_[j,:])[0])
##        print len(corr_list)
        best_ind=numpy.argmax(abs( numpy.array( corr_list)))
        print numpy.sign(corr_list[best_ind])

##        print best_ind
        print pos[rel_ind[best_ind]]
        pyplot.scatter(rescaled[:,rel_ind[best_ind]], pca.components_[j,:],c=colors)
        pyplot.ylabel('PC{0} Loading'.format(j+1))
        pyplot.xlabel('Centered Cons Allele Proportion at {0}'.format(pos[rel_ind[best_ind]]))
##        pyplot.gca().set_aspect('equal')
        pyplot.show()

def PlotFreqByPopulation(table,var1, pos1, offset):
    color_dict={'firebrick':.5, 'forestgreen':1.5, 'blue':2.5, "orange":3.5,'purple':4.5}
##    ind=numpy.array(good_ind)
    nt_dict={'A':1, 'T':2, 'C':3,'G':4}
    array=numpy.copy( table[:,:, pos1])
    array[:, nt_dict[var1]]+=array[:,0]
    array=array[:,1:]
    array/=array.sum(1)[:,None]

##    freq1=table[ind, nt_dict[var1], pos1]/table[ind, :, pos1].sum(1)
    freq1=array[:,nt_dict[var1]-1]
    print freq1.shape

##    print numpy.log10( freq1+offset)
    color=numpy.array( colors)
    pop_dict=[]
    color_set=sorted(list(set(color)))
    pop_freq=[]
    for c in color_set:
        ind=(color==c)
        pop_freq.append(numpy.log10( freq1+offset))
##    seaborn.violinplot(data=pop_freq)
    jitter=numpy.random.uniform(-.2, .2, size=len(ind))
    x_pos=numpy.array( [color_dict[c] for c in color])
##    print numpy.log10( freq1+offset)
##    print x_pos
##    pyplot.scatter(x_pos+jitter+.5, numpy.log10( freq1+offset),c=color, s=80, alpha=.8)
    pyplot.scatter(x_pos+jitter+.5,  freq1,c=color, s=80, alpha=.8)
##    pyplot.ylabel('Cons Freq at Pos {0} + 1e-3.5 (Log10)'.format(pos1), size=20)
    pyplot.ylabel('Cons Freq at Pos {0}'.format(pos1), size=20)
    pyplot.xticks(numpy.arange(1,6),['B','I','N','T','Z'], size=20)
    pyplot.yticks(size=20)
    pyplot.xlim(.5, 5.5)
##    pyplot.ylim(-3.6, .1)
##    pyplot.ylim(-.05, 1.05)
    pyplot.title(var1)
##    else: color='firebrick'
##    pyplot.show()

def DetermineConsensusAndVariantAllele(array, cons_allele):
    ind=numpy.array(good_ind)
    nt_dict={'A':1, 'T':2, 'C':3,'G':4}
    array=numpy.copy( table[ind,:, pos1])
    array[:, nt_dict[var1]]+=array[:,0]
    array=array[:,1:]
    array/=array.sum(1)[:,None]

##    freq1=table[ind, nt_dict[var1], pos1]/table[ind, :, pos1].sum(1)
    freq1=array[:,nt_dict[var1]-1]

def SelectableBiplot():
    """A biplot that can be interacted with by clicking arrows to plot the
    consensus allele frequency at a position."""
    pass

##    pyplot.show()



def PlotAlleleProfiles(freq, pos, lin_fit):
    ind=numpy.array(good_ind)
    rescaled=AlleleRescale(freq[ind,:])

    for i in range(3):
        positive=lin_fit[i].coef_!=0
        std=numpy.std( rescaled[:, positive],0)
        print std.shape
        f=pyplot.plot(rescaled[:, positive]/std[None,:], c='red', alpha=.1)
        pyplot.show()
        corr=numpy.corrcoef(rescaled[:, positive].transpose())
        seaborn.clustermap(corr)
        pyplot.show()
def PlotHistone(hist_freq, hist_pos, lin_fit, PCA, colors=colors, pc=0, hist=False, highlight=[], outlier=[] ):
##    ind=numpy.array(good_ind)
    color_dict={'B':'firebrick', 'I':'forestgreen', "N":'blue', 'T':"orange",'Z':'purple'}
    sample_dict={v:k for v,k in color_dict.iteritems()}
    Samples=numpy.array(samples)[numpy.array(good_ind)]

##    rel_ind=(lin_fit[1].coef_!=0)
##    print sum(rel_ind)
##/    stacked_coef=numpy.vstack(( lin_fit[0].coef_[rel_ind],lin_fit[1].coef_[rel_ind],lin_fit[2].coef_[rel_ind]))
    if hist==True:
        pyplot.plot([129,930],[-.5,-.5], c='r', lw=2)
        pyplot.plot([1932,2309],[-.5,-.5], c='r', lw=2)
        pyplot.plot([1336,1708],[-.5,-.5], c='r', lw=2)
        pyplot.plot([2783,3095],[-.5,-.5], c='r', lw=2)
        pyplot.plot([3391,3802],[-.5,-.5], c='r', lw=2)
####    return ()
##    pyplot.scatter(hist_pos[lin_fit[pc].coef_<.0], [-.1]*len(hist_pos[lin_fit[pc].coef_<.0]),s=60, c='orange', marker='v')
##    pyplot.scatter(hist_pos[lin_fit[1].coef_<0], [-.13]*len(hist_pos[lin_fit[1].coef_<0]),s=60, c='green', marker='v')
##    pyplot.scatter(hist_pos[lin_fit[2].coef_<.0], [-.16]*len(hist_pos[lin_fit[2].coef_<.0]),s=60, c='blue',marker='v')
##    pyplot.scatter(hist_pos[lin_fit[pc].coef_>.0], [-.1]*len(hist_pos[lin_fit[pc].coef_>.0]),s=60, c='orange', marker='^')
##    pyplot.scatter(hist_pos[lin_fit[1].coef_>0], [-.13]*len(hist_pos[lin_fit[1].coef_>0]),s=60, c='green', marker='^')
##    pyplot.scatter(hist_pos[lin_fit[2].coef_>.0], [-.16]*len(hist_pos[lin_fit[2].coef_>.0]),s=60, c='blue',marker='^')
##    pyplot.scatter(hist_pos[lin_fit[pc].scores_!=0], [-.1]*len(hist_pos[lin_fit[pc].scores_!=0]),s=60, c='orange')
##    pyplot.scatter(hist_pos[lin_fit[1].scores_!=0], [-.13]*len(hist_pos[lin_fit[1].scores_!=0]),s=60, c='green')
##    pyplot.scatter(hist_pos[lin_fit[2].scores_!=0], [-.16]*len(hist_pos[lin_fit[2].scores_!=0]),s=60, c='blue')
    offset=-.2
##    rescale=AlleleRescale(hist_freq[ind,:])
    rescale=hist_freq[:,:]
    for c in list(set(colors)):
##        f1=[ pyplot.scatter(hist_pos[rel_ind]+offset, hist_freq[i,rel_ind], c=c, alpha=.7) for i in range(len(colors)) if colors[i]==c]
        outlier=highlight!=0
        if outlier!=[]:
            f1=[ pyplot.scatter(hist_pos[:][outlier] +offset, rescale [i,:][outlier], c=c, edgecolor='dodgerblue', s=80,  alpha=.7) for i in range(len(colors)) if colors[i]==c]
        f1=[ pyplot.scatter(hist_pos[:]+offset, rescale [i,:], c=c, alpha=.7) for i in range(len(colors)) if colors[i]==c]
        offset+=.1
    pyplot.ylabel('Consensus Allele Frequency')
    pyplot.xlabel('Consensus Position')
##    cutoff
##    for  i in range(len(lin_fit)):
    return
    rel_ind=((lin_fit[0].coef_!=0)+ (lin_fit[1].coef_!=0) +(lin_fit[2].coef_!=0))>0
##    rel_ind=((lin_fit[0].coef_!=0)+ (lin_fit[1].coef_!=0) +(lin_fit[2].coef_!=0)+(lin_fit[3].coef_!=0)+ (lin_fit[4].coef_!=0) +(lin_fit[5].coef_!=0))>0
##    rel_ind=((lin_fit[0].coef_!=0)+ (lin_fit[1].coef_!=0) )>0
    print sum(rel_ind)
    stacked_coef=numpy.vstack(( lin_fit[0].coef_[rel_ind],lin_fit[1].coef_[rel_ind],lin_fit[2].coef_[rel_ind]))
##    stacked_coef=numpy.vstack(( lin_fit[0].coef_[rel_ind],lin_fit[1].coef_[rel_ind],lin_fit[2].coef_[rel_ind], lin_fit[3].coef_[rel_ind],lin_fit[4].coef_[rel_ind],lin_fit[5].coef_[rel_ind]))
    print stacked_coef.shape
##    rel_ind=((lin_fit[0].scores_!=0)+ (lin_fit[1].scores_!=0) +(lin_fit[2].scores_!=0))>0
    pyplot.show()
##    return ()
    ax=pyplot.subplot(411)
    seaborn.heatmap(hist_freq[ind,:] [:, rel_ind], xticklabels=hist_pos[rel_ind], cbar_kws={'label':'Consensus Freq'})
##    return rel_ind
    yticks=pyplot.yticks( numpy. arange(len(Samples)-1,-1,-1)+.5 , ['-']*len(Samples), size=20)
    f=[yticks[1][i].set_color(color_dict[Samples[i][0]]) for i in range(len (Samples))]
    f=[yticks[1][numpy.where(Samples== i)[0]].set_color('black') for i in highlight]

##    pyplot.subplot(212, sharex=ax)
##    seaborn.heatmap(stacked_coef/stacked_coef.sum(1)[:,None], xticklabels=hist_pos[rel_ind], yticklabels=['PC1', 'PC2','PC3'], cbar_kws={'label':'Coef'})
    for i in range(0,3):
        pyplot.subplot(412+i, sharex=ax)
        print stacked_coef[i,:]
        seaborn.heatmap(stacked_coef[i,:][None,:]*0.5, annot_kws={'size':14}, yticklabels=['PC1', 'PC2','PC3'][i],cmap='Greys', cbar=False, mask=abs( (stacked_coef[i,:][None,:])!=0), center=0.)
        seaborn.heatmap(stacked_coef[i,:][None,:],xticklabels=hist_pos[rel_ind], annot_kws={'size':14}, yticklabels=['PC1', 'PC2','PC3'][i], cmap='RdBu_r', cbar_kws={'label':'Coef'}, mask=(stacked_coef[i,:][None,:]==0), center=0.)
##        seaborn.heatmap(stacked_coef[i,:][None,:], xticklabels=hist_pos[rel_ind], yticklabels=['PC1', 'PC2','PC3'][i], cbar_kws={'label':'Coef'}, mask=(stacked_coef[i,:][None,:]!=0))
##    pyplot.xticks(numpy.arange(len(rel_ind))+.5, hist_pos[rel_ind], size=14)
    pyplot.show()
##    pca_part=PlotPCA(hist_freq[:, rel_ind])
    pca_rem=PlotPCA(hist_freq[:, ~rel_ind])
##    for i in range(3):
##        print scipy.stats.pearsonr(pca_part.components_[i,:], PCA.components_[i,:])
##    for i in range(3):
##        print scipy.stats.pearsonr(pca_rem.components_[i,:], PCA.components_[i,:])
    return stacked_coef
def TestPCA(pca, data,lin_fit, cutoffs):
    r=[]
    for c in cutoffs:
        ind=((lin_fit[0].coef_>=c)+(lin_fit[1].coef_>=c)+(lin_fit[2].coef_>=c)>0)
        pca_test=sklearn.decomposition.PCA(3)
        pca_test.fit(data[:, ind].transpose())
        rss=sum(sum((pca.components_[:3,:]-pca_test.components_[:3,:])**2))**.5
        print sum(ind)
        print rss
        r.append(rss)
    return(r)

def SemiPartialRS(pca, data, lin_fit_full, N=1):
    predict_full=lin_fit_full.predict(data)
    RSS_full= ((pca.components_[0,:]- predict_full)**2).sum()**.5
    ordering=[]
    R_s=[]
    lin_fit=sklearn.linear_model.LinearRegression()
    ind=numpy.arange(data.shape[1])
    reg=[]
    for i in range(len(ind)):
        lin_fit.fit(data[:, ind!=i], pca.components_[0,:])
        predict=lin_fit.predict(data[:, ind!=i])

##    pyplot.scatter(pca.components_[0,:], predict[:,0])
        RSS= ((pca.components_[0,:]- predict)**2).sum()**.5
        reg.append(RSS_full-RSS)

    return reg



def GetCorrelation(data,pca):
    corr1=[]
    corr2=[]
    samples, pos=data.shape
    for i in range(pos):
        corr1.append(numpy.corrcoef(data[:,i], pca.components_[0,:])[0,1])
        corr2.append(numpy.corrcoef(data[:,i], pca.components_[1,:])[0,1])
    return corr1, corr2

def PlotSparsePCA(data, colors,samples=''):
    if samples!='':
        print 'b'
        keep_ind=numpy.array( [bad_strains.count(s)==0 for s in samples])
        colors=numpy.array(colors)[keep_ind]
        data=data[ keep_ind,:]
        nan_ind=numpy.isnan(data.sum(1))
##        data=data[~nan_ind,:]
        pca=sklearn.decomposition.SparsePCA(92-len(bad_strains), alpha=.1)
    else:

        pca=sklearn.decomposition.PCA(92,whiten=True)
    trans= pca.fit_transform(AlleleRescale( data))
    pyplot.scatter(trans[:,0], trans[:,1], c=colors, s=80, alpha=.7)
##    pyplot.xlabel('PC1 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[0]*100, 3)))
##    pyplot.ylabel('PC2 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[1]*100, 3)))
    pyplot.show()

##    pyplot.scatter(pca.components_[1,:], pca.components_[2,:], c=colors, s=80, alpha=.7)
    pyplot.scatter(trans[:,1], trans[:,2], c=colors, s=80, alpha=.7)
##    pyplot.xlabel('PC2 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[1]*100, 3)))
##    pyplot.ylabel('PC3 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[2]*100, 3)))
    pyplot.show()


    fig=pyplot.figure()
    ax=mplot3d.axes3d.Axes3D(fig)
    ax.scatter(trans[:,0],trans[:,1], trans[:,2], c=colors, s=60)
##    ax.set_xlabel('PC1 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[0]*100, 3)))
##    ax.set_ylabel('PC2 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[1]*100, 3)))
##    ax.set_zlabel('PC3 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[2]*100, 3)))
    pyplot.show()

##    pyplot.plot(range(1 , len( pca.explained_variance_)+1),pca.explained_variance_ratio_, lw=2, c='forestgreen', alpha=.9)
##    pyplot.ylabel('Variance Explained')
##    pyplot.xlabel('Component')
    pyplot.plot(pca.components_[0,:],c='r', alpha=.8, lw=2, label='PC1')
    pyplot.plot(pca.components_[1,:],c='g', alpha=.8, lw=2, label='PC2')
    pyplot.plot(pca.components_[2,:],c='orange', alpha=.8, lw=2, label='PC3')
##    pyplot.plot(pca.components_[3,:],c='g', alpha=.8, lw=2, label='PC4')
##    pyplot.plot(pca.components_[4,:],c='orange', alpha=.8, lw=2, label='PC5')
    pyplot.legend()
    pyplot.show()

    pyplot.show()

    return pca


import sklearn.preprocessing
def PlotPCATxn(data, colors,samples=''):
    scaled=sklearn.preprocessing.robust_scale(data,0)
    pca=sklearn.decomposition.PCA(92,whiten=True)
    pca.fit(AlleleRescale( scaled).transpose())
    pyplot.scatter(pca.components_[0,:], pca.components_[1,:], c=colors, s=80, alpha=.7)
    pyplot.xlabel('PC1 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[0]*100, 3)))
    pyplot.ylabel('PC2 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[1]*100, 3)))
    pyplot.show()

    pyplot.scatter(pca.components_[1,:], pca.components_[2,:], c=colors, s=80, alpha=.7)
    pyplot.xlabel('PC2 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[1]*100, 3)))
    pyplot.ylabel('PC3 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[2]*100, 3)))
    pyplot.show()


    fig=pyplot.figure()
    ax=mplot3d.axes3d.Axes3D(fig)
    ax.scatter(pca.components_[0,:], pca.components_[1,:],pca.components_[2,:], c=colors, s=60)
    ax.set_xlabel('PC1 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[0]*100, 3)))
    ax.set_ylabel('PC2 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[1]*100, 3)))
    ax.set_zlabel('PC3 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[2]*100, 3)))
    pyplot.show()

    pyplot.plot(range(1 , len( pca.explained_variance_)+1),pca.explained_variance_ratio_, lw=2, c='forestgreen', alpha=.9)
    pyplot.ylabel('Variance Explained')
    pyplot.xlabel('Component')
    pyplot.show()

    return pca

def ClusterMap(data, colors,samples):
    samples=numpy.array(samples)
    keep_ind=numpy.array( [bad_strains.count(s)==0 for s in samples])
    colors=numpy.array(colors)[keep_ind]
    data=data[keep_ind,:]
    nan_ind=numpy.isnan(data.sum(1))
##    data=data[~nan_ind,:]
##    corr=numpy.corrcoef(AlleleRescale(data))
    corr=AlleleRescale(data)
    seaborn.clustermap(corr, xticklabels=samples[keep_ind], yticklabels=samples[keep_ind])
    pyplot.show()

def Replace_Marginal(vertical,real_x=[]):
    pyplot.plot(real_x)


def ManualJointPlot():
    gs=pyplot.GridSpec(3,3)
    jnt=pyplot.subplot(gs[ 1:,:-1])
    marg_y=pyplot.subplot(gs[ 1:,-1], sharey=jnt)
    marg_x=pyplot.subplot(gs[0,:-1], sharex=jnt)
    return jnt,marg_x, marg_y

def PlotPCAKernelDensities(data, colors, samples):
    keep_ind=numpy.array( [bad_strains.count(s)==0 for s in samples])
    colors=numpy.array(colors)
    data=data[:, keep_ind]
    nan_ind=numpy.isnan(data.sum(1))
    data=data[~nan_ind,:]
    jnt,marg_x, marg_y=ManualJointPlot()
    pca=sklearn.decomposition.PCA(92-len(bad_strains),whiten=True)
    pca.fit(data)
##    jnt.set_aspect('equal')
    jnt.scatter(pca.components_[0,:], pca.components_[1,:], c=colors, s=80, alpha=.7)
    jnt.set_xlabel('PC1 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[0]*100, 3)))
    jnt.set_ylabel('PC2 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[1]*100, 3)))
    color_set=set(colors)
##    ax1=pyplot.subplot(211)
##    ax2=pyplot.subplot(212)
    for c in color_set:

        pop_data=pca.components_[:2,colors==c]
        print pca.components_.shape
        print pop_data
        print colors
        print colors==c
        std=numpy.std(pop_data,1)
        print std
        print numpy.max(pop_data,0)
        pop_max,pop_min=numpy.max(pop_data,1)+std*2, numpy.min(pop_data,1)-std*2
        scaled_data=pop_data/std[:,None]
        cv1=sklearn.model_selection.GridSearchCV(sklearn.neighbors.KernelDensity(),param_grid={ 'bandwidth':numpy.logspace(-2,2,500)})
        cv1.fit(pop_data[0,:][:,None])
##        bw_x=cv.best_estimator_.bandwidth*std[0]
##        bw_y=cv.best_estimator_.bandwidth*std[1]
##        kde_1=sklearn.neighbors.KernelDensity(bw_x)
        print pop_data[0,:], pop_data[1,:]
##        kde_1.fit(pop_data[0,:][:,None])
        kde_1=cv1.best_estimator_

        cv2=sklearn.model_selection.GridSearchCV(sklearn.neighbors.KernelDensity(),param_grid={ 'bandwidth':numpy.logspace(-2,2,500)})
        cv2.fit(pop_data[1,:][:,None])

##        kde_2=sklearn.neighbors.KernelDensity(bw_y)
##
##        kde_2.fit(pop_data[1,:][:,None])
        kde_2=cv2.best_estimator_
        mgrid=numpy.mgrid
##        return kde_1
        print numpy.arange(pop_min[0],pop_max[0],(pop_max[0]-pop_min[0])/ 1000.)
        pop_x=numpy.exp( kde_1.score_samples(numpy.arange(pop_min[0],pop_max[0],(pop_max[0]-pop_min[0])/ 1000.)[:,None]))
        pop_y=numpy.exp( kde_2.score_samples(numpy.arange(pop_min[1],pop_max[1],(pop_max[1]-pop_min[1])/ 1000.)[:,None]))
        marg_x.plot(numpy.arange(pop_min[0],pop_max[0],(pop_max[0]-pop_min[0])/ 1000.), pop_x,lw=2,alpha=.8, c=c)
        marg_y.plot( pop_y,numpy.arange(pop_min[1],pop_max[1],(pop_max[1]-pop_min[1])/ 1000.), c=c,lw=2,alpha=.8)
    pyplot.show()

def JointFreqSpectrum(data, pop_assignments):
    pops=sorted(list(set(pop_assignments)))
    count=2
    ax=pyplot.subplot('611')
    for p in pops:
        pyplot.subplot('61{0}'.format(count), sharex=ax)
        sliced=data[pop_assignments==p,:].flatten()
        seaborn.distplot(sliced, kde=False, bins=200)
        count+=1
    pyplot.show()


def AlleleRescale(data,axis=0, std=False):
    mean_f=numpy.nanmean(data,axis)
##    mean_f=data
    var_f=(mean_f*(1-mean_f))**.5
    if std==True: var_f=numpy.std(data,axis)
    if axis==0:
        return (data-mean_f)
    else:
        return (data-mean_f[:, None])

def PlotIsomap(data, colors, param_range=range(2, 50)):
    fits, scores, mats=Isomap(data, param_range, [2,3])
    best_ind=numpy.argmin(scores)
    best_fit=fits[best_ind]
    pyplot.scatter(best_fit.embedding_[:,0], best_fit.embedding_[:,1], c=colors, s=80, alpha=.7)
##    pyplot.xlabel('PC1 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[0]*100, 3)))
##    pyplot.ylabel('PC2 ({0}%)'.format(numpy.round( pca.explained_variance_ratio_[1]*100, 3)))
    pyplot.show()
    fig=pyplot.figure()
    ax=mplot3d.axes3d.Axes3D(fig)
    ax.scatter(best_fit.embedding_[:,0], best_fit.embedding_[:,1], best_fit.embedding_[:,2], c=colors, s=60)
    pyplot.show()

def PlotNtComposition(table, position, seq):
    nt_dict={'A':1, 'T':2, 'C':3,'G':4}
    color_dict={'B':'firebrick', 'I':'forestgreen', "N":'blue', 'T':"orange",'Z':'purple'}
    nt=seq[position-1]
    print nt
    array=numpy.copy( table[:,:, position])
    array[:, nt_dict[nt]]+=array[:,0]
    array=array[:,1:]
    array/=array.sum(1)[:,None]
    seaborn.heatmap(array )
    xticks=pyplot.xticks(numpy.arange(4)+.5, ['A','T','C','G'], size=40)
    xticks[1][nt_dict[nt]-1].set_color('r')
    yticks=pyplot.yticks( numpy. arange(len(samples)-1,-1,-1)+.5 , ['-']*len(samples), size=20)
    f=[yticks[1][i].set_color(color_dict[samples[i][0]]) for i in range(len (samples))]
##    pyplot.show()
##    pyplot.savefig("c:/barbashlab/rdna_pos_{0}.jpeg".format(position))
##    pyplot.close()

def MakeDir(newdir):
    if os.path.exists(newdir)==False:
        os.mkdir(newdir)


def PlotOutliers(freq):
    rescaled=AlleleRescale( AlleleRescale(freq),1)
    pca=sklearn.decomposition.PCA()
    pca.fit(rescaled)
##    cv=sklearn.model_selection.GridSearchCV
    iso=sklearn.ensemble.IsolationForest()
##    cv=sklearn.model_selection.GridSearchCV(iso, {'contamination':10.**numpy.arange(-3,-1,.5)},scoring='adjusted_rand_score')
    iso.fit(rescaled.transpose())

    pred=iso.predict(rescaled.transpose())
##    print cv.best_params_
    print pca.components_.shape
    print pred.shape
    pyplot.scatter(pca.components_[0,:], pca.components_[1,:], c=pred, alpha=.5)
    pyplot.show()
    return iso.decision_function(rescaled.transpose())
def FindRelevantPositions(lin_fit):
    rel_ind=((lin_fit[0].coef_!=0)+ (lin_fit[1].coef_!=0))>0 #+(lin_fit[2].coef_!=0))>0
    return rel_ind

def PlotAllRelevantPositions( outfile, table, positions, lin_fit, seq):
    rel_ind=((lin_fit[0].coef_!=0)+ (lin_fit[1].coef_!=0) +(lin_fit[2].coef_!=0))>0
##    print sum(rel_ind)
    stacked_coef=numpy.vstack(( lin_fit[0].coef_[rel_ind],lin_fit[1].coef_[rel_ind],lin_fit[2].coef_[rel_ind]))
    relevant_positions=positions[rel_ind]
    MakeDir(outfile)
    for pos in relevant_positions:
        plot=PlotNtComposition(table, pos, seq)
        pyplot.savefig("{0}/{1}.jpeg".format(outfile, pos))
        pyplot.close()

def PlotAllRelevantPositionsScatter( outfile, table, positions, lin_fit, seq):
    rel_ind=((lin_fit[0].coef_!=0)+ (lin_fit[1].coef_!=0) +(lin_fit[2].coef_!=0))>0
##    print sum(rel_ind)
    stacked_coef=numpy.vstack(( lin_fit[0].coef_[rel_ind],lin_fit[1].coef_[rel_ind],lin_fit[2].coef_[rel_ind]))
    relevant_positions=positions[rel_ind]
    MakeDir(outfile)
    for pos in relevant_positions:
        nt=seq[pos-1]
        plot=PlotFreqByPopulation( table,nt, pos, 10.**-3.5)
        pyplot.savefig("{0}/{1}.jpeg".format(outfile, pos))
        pyplot.close()



def ClusterIsomap(data, colors, param_range=range(2, 50), k_range=range(2, 8)):
    fits, scores, mats=Isomap(data, param_range, [2,3])
    best_ind=numpy.argmin(scores)
    best_fit=fits[best_ind]
    label_list=[]
    score_list=[]
    for k in k_range:
        GMM=sklearn.mixture.GaussianMixture(k)
        GMM.fit(best_fit.embedding_)
        labels= GMM.predict(best_fit.embedding_)
        label_list.append(labels)
        score_list.append(metrics.silhouette_score(best_fit.embedding_,labels))

    clust_ind=numpy.argmax(score_list)
    colors=label_list[clust_ind]
    pyplot.scatter(best_fit.embedding_[:,0], best_fit.embedding_[:,1], c=colors, s=80, alpha=.7, cmap='gist_rainbow')
    pyplot.show()
    fig=pyplot.figure()
    ax=mplot3d.axes3d.Axes3D(fig)
    ax.scatter(best_fit.embedding_[:,0], best_fit.embedding_[:,1], best_fit.embedding_[:,2], c=colors, cmap='gist_rainbow', s=60)
    pyplot.show()
    return colors

def GetFreqArray(snp_array, cutoff=10):
    cvg=snp_array.sum(1)
    min_cvg=numpy.min(cvg,0)
##    pyplot.plot(min_cvg)
##    pyplot.show()
    good_ind=min_cvg>=cutoff
    freq_array=snp_array[:,0,good_ind]/snp_array[:,:,good_ind].sum(1)
    positions=numpy.arange(snp_array.shape[2])[good_ind]

    return freq_array, positions

def GetMajorAlleleFreq(snp_array, cutoff=10):
    cvg=snp_array.sum(1)
    min_cvg=numpy.min(cvg,0)
##    pyplot.plot(min_cvg)
##    pyplot.show()
    good_ind=min_cvg>=cutoff
    print good_ind
    freq_table=snp_array/snp_array.sum(1)[:,None,:]
    pop_avg=numpy.mean(freq_table,0)
    ind=numpy.argmax(pop_avg,0)

    pop_avg=numpy.mean(freq_table,0)
    ind=numpy.argmax(pop_avg,0)

    positions=numpy.arange(snp_array.shape[2])[good_ind]


    return numpy.array([ freq_table[:,ind[ i],i] for i in numpy.where(good_ind==True)[0]]).transpose(),ind, positions

def GetOtherAlleleProportion(snp_array, cutoff=10, order=2):
    cvg=snp_array.sum(1)
    min_cvg=numpy.min(cvg,0)
##    pyplot.plot(min_cvg)
##    pyplot.show()
    good_ind=min_cvg>=cutoff
    freq_table=snp_array/snp_array.sum(1)[:,None,:]
    pop_avg=numpy.mean(freq_table,0)
    ind=numpy.argsort(pop_avg,0)[-1*order, :]

##    pop_avg=numpy.mean(freq_table,0)
####    print pop_avg.shape
####    print jabber
##    ind=numpy.argmax(pop_avg,0)

    positions=numpy.arange(snp_array.shape[2])[good_ind]


    return numpy.array([ freq_table[:,ind[ i],i] for i in numpy.where(good_ind==True)[0]]).transpose()


def GetSecondaryAlleleFreq(snp_array, cutoff=10):

    table=numpy.copy(snp_array)
    sample_num,nt_num,pos_num=table.shape
    ind= numpy.argmax(snp_array,1)[:,-2,:]
    read_counts-table.sum(1)
    freq_array=numpy.ndarray((sample_num,pos_num))
    freq_array.fill(0.)
    for i in pos_num:
        freq_array[:,i]=table[s, ind[s,:],:]
    return freq_array

import sys
import pandas
import sklearn.preprocessing

#Simulation to model gene conversion

def GeneConvOverTime(var_num , CN, rate, generations):
    prop=[]
    for t in range(generations):
        number_of_events=scipy.stats.poisson.rvs(rate)
        for event in range(number_of_events):
            p=float(var_num)/CN
            q=1.-p
            prob=2*p*q
            bern_trial=scipy.stats.bernoulli.rvs(prob)
            if bern_trial==1:
                step=(scipy.stats.bernoulli.rvs(.5)-.5)*2
                var_num+=step
        prop.append(float(var_num)/CN)
    return prop

def UEOverTime(var_num , CN, rate, generations):
    prop=[]
    array=['A']*CN
##    indices=numpy.random.choice(numpy.arange(CN), replace=False, size=var_num)
    indices=numpy.arange(int( CN/10),int( CN/10)+var_num)
    for i in indices:
        array[i]='B'
    CN_list=[]
    for t in range(generations):
        #Gene Conversion

        number_of_events=scipy.stats.poisson.rvs(rate)
        for n in range(number_of_events):
            i=numpy.random.randint(0,len(array))
            j=numpy.random.randint(0,len(array))
            array[i]=array[j]

        i=numpy.random.randint(1,len(array)+1)
        j=numpy.random.randint(1,len(array)+1)
        if len( array[:i]+array[j:]) >100 and len( array[:i]+array[j:]) <2000:array=array[:i]+array[j:]
        prop.append (float( array.count('A'))/len(array))
        CN_list.append(len(array))
##        if CN_list[-1]<100:
##            CN_list[-1]=CN_list[-2]

    return prop, CN_list



def UEOverTime(var_num , CN, rate, generations):
    prop=[]
    array=['A']*CN
    indices=numpy.random.choice(numpy.arange(CN), replace=False, size=var_num)
##    indices=numpy.arange(int( CN/10),int( CN/10)+var_num)
    for i in indices:
        array[i]='B'
    CN_list=[]
    for t in range(generations):
        #Gene Conversion

        number_of_events=scipy.stats.poisson.rvs(rate)
        for n in range(number_of_events):
            i=numpy.random.randint(0,len(array))
            j=numpy.random.randint(0,len(array))
            array[i]=array[j]

        i=numpy.random.randint(1,len(array)+1)
        j=numpy.random.randint(1,len(array)+1)
        if len( array[:i]+array[j:]) >100 and len( array[:i]+array[j:]) <2000:array=array[:i]+array[j:]
        prop.append (float( array.count('A'))/len(array))
        CN_list.append(len(array))
##        if CN_list[-1]<100:
##            CN_list[-1]=CN_list[-2]

    return prop, CN_list


def CLTNorm(x, obs, mean,std):
    return scipy.stats.norm.pdf(int( obs), x*mean, std*x**.5)

def VarianceAtRelevantFeatures(freq, pos, lin_fit):
    rel_ind=numpy.where(( ((lin_fit[0].coef_!=0)+ (lin_fit[1].coef_!=0) +(lin_fit[2].coef_!=0))>0)==True)[0]
    ind=numpy.array(good_ind)
    rescaled=AlleleRescale(freq[ind,:])
    return numpy.std(rescaled[:,rel_ind],0)

def FstAtRelevantFeatures(freq, pos, lin_fit):
    rel_ind=numpy.where(( ((lin_fit[0].coef_!=0)+ (lin_fit[1].coef_!=0) +(lin_fit[2].coef_!=0))>0)==True)[0]
    ind=numpy.array(good_ind)
    rescaled=AlleleRescale(freq[ind,:])
    total_var= numpy.var(rescaled[:,rel_ind],0)
    color_set=sorted(list(colors))
    color_array=numpy.array(colors)
    sub_var=[]
    for c in color_set:
        color_ind= (color_array==c)
        pop_var=numpy.var(rescaled[:,rel_ind][color_ind,:],0 )
        sub_var.append(pop_var)
    sub_var=numpy.array(sub_var)
    return numpy.mean(( total_var- sub_var)/total_var,0)

def FstAtIrrelevantFeatures(freq, pos, lin_fit):
    rel_ind=numpy.where(( ((lin_fit[0].coef_!=0)+ (lin_fit[1].coef_!=0) +(lin_fit[2].coef_!=0))>0)==True)[0]
    ind=numpy.array(good_ind)
    rescaled=AlleleRescale(freq[ind,:])
    total_var= numpy.var(rescaled[:,~rel_ind],0)
    color_set=sorted(list(colors))
    color_array=numpy.array(colors)
    sub_var=[]
    for c in color_set:
        color_ind= (color_array==c)
        pop_var=numpy.var(rescaled[:,~rel_ind][color_ind,:],0 )
        sub_var.append(pop_var)
    sub_var=numpy.array(sub_var)
##    print sub_var
    return numpy.mean(( total_var- sub_var)/total_var,0)

def SequenceDiversity(ind1, ind2):
    total_count_1=ind1.sum(0)
    total_count_2=ind2.sum(0)
    allele_freq_1=ind1/total_count_1
    allele_freq_2=ind2/total_count_2
    nt, length=allele_freq_1.shape
    div=numpy.array([0.]*length)
    for n in range(nt):
        div+=allele_freq_1[n,:]*allele_freq_2[n,:]

    return 1.-div

##def ComputeFst(snp_array,positions, pop  ):
##    pi_within=[]
##    pi_between=[]
##    for j in range(len(pop)):
##        for i in range(j):
##            if pop[i]==pop[j]:
##                pi_within.append(SequenceDiversity(snp_array[i,:,:], snp_array[j,:,:])[positions])
##            elif pop[i]!=pop[j]:
##                pi_between.append(SequenceDiversity(snp_array[i,:,:], snp_array[j,:,:])[positions])
##    pi_between=numpy.array(pi_between)
##    pi_within=numpy.array(pi_within)
##    mean_between=numpy.mean(pi_between,0)
##    mean_within=numpy.mean(pi_within,0)
##    fst=(mean_between-mean_within)/mean_between
##    return fst
##
##def ComputeFst_efficient(snp_array,positions, pop  ):
##    pi_within=numpy.array([0.]* len(positions))
##    pi_between=numpy.array([0.]* len(positions))
##    n_within=0.
##    n_between=0.
##    for j in range(len(pop)):
##        for i in range(j):
##            if pop[i]==pop[j]:
##                pi_within+=SequenceDiversity(snp_array[i,:,:], snp_array[j,:,:])[positions]
####                pi_between+=SequenceDiversity(snp_array[i,:,:], snp_array[j,:,:])[positions]
####                n_between+=1
##                n_within+=1
##            elif pop[i]!=pop[j]:
##                pi_between+=SequenceDiversity(snp_array[i,:,:], snp_array[j,:,:])[positions]
##                n_between+=1
##    pi_between=numpy.array(pi_between)
##    pi_within=numpy.array(pi_within)
##    mean_between=pi_between/n_between
##    mean_within=pi_within/n_within
##    fst=(mean_between-mean_within)/(mean_between)
##    return fst
##
##
##
##def Fst(freq,ax=0):
##
##    rescaled=AlleleRescale(freq[:,:])
##    total_var= numpy.var(rescaled[:,:],0)
##    print total_var.shape
##    color_set=sorted(list(set( colors)))
##    color_array=numpy.array(colors)
##    sub_var=[]
##    for c in color_set:
##        color_ind= (color_array==c)
##        pop_var=numpy.var(rescaled[:,:][color_ind,:],0 )
##        sub_var.append(pop_var)
##    sub_var=numpy.array(sub_var)
##    print sub_var.shape
##
####    print pw_list.shape
##    return numpy.max(sub_var,ax)

def VarianceAtIrrelevantFeatures(freq, pos, lin_fit):
    rel_ind=numpy.where(( ((lin_fit[0].coef_!=0)+ (lin_fit[1].coef_!=0) +(lin_fit[2].coef_!=0))>0)==True)[0]
    ind=numpy.array(good_ind)
    rescaled=AlleleRescale(freq[ind,:])
    return numpy.std(rescaled[:,~rel_ind],0)
def MeanAtRelevantFeatures(freq, pos, lin_fit):
    rel_ind=numpy.where(( ((lin_fit[0].coef_!=0)+ (lin_fit[1].coef_!=0) +(lin_fit[2].coef_!=0))>0)==True)[0]
    ind=numpy.array(good_ind)
    return numpy.mean(freq[ind,:] [:,rel_ind],0)


def MeanAlleleProportion(table, ax=0):
    #Find major allele
    freq_table=table/table.sum(1)[:,None,:]
    pop_avg=numpy.mean(freq_table,0)

    return numpy.max(pop_avg,0)

def VarAlleleProportion(table, ax=0):
    #Find major allele
    freq_table=table/table.sum(1)[:,None,:]
    pop_avg=numpy.mean(freq_table,0)
    ind=numpy.argmax(pop_avg,0)
##    print freq_table[:,numpy.array( ind),:].shape
    return numpy.array([ numpy.std(freq_table[:,ind[ i],i],0) for i in range(len(ind))])

def PlotVarianceBoxplot(samples,lin_fits, labels):

    variances=[]
    for i in range(len(samples)):
        variances.append(VarianceAtIrrelevantFeatures(samples[i],0, lin_fits[i]))
    pyplot.boxplot(variances, widths=.8)
    for i in range(len(variances)):
        jitter=numpy.random.uniform(-.4,.4,size=len(variances[i]))
        pyplot.scatter(jitter+i+1, variances[i], alpha=.3)
    pyplot.xticks(numpy.arange(1,len(labels)+1), labels, size=14)
    pyplot.ylabel('Standard Deviation', size=14)
    pyplot.show()

def PlotrDNA(y, space_size=4, SU_size=4):
    ETS='grey'
    ITS='grey'
    SU='black'
    pyplot.plot([0,864],[y,y], c=ETS, lw=space_size)
    pyplot.plot([864,2859],[y,y], c=SU, lw=SU_size)
    pyplot.plot([2859,3585],[y,y], c=ITS, lw=space_size)
    pyplot.plot([3585,3708],[y,y], c=SU, lw=SU_size)
    pyplot.plot([3708,3736],[y,y], c=ITS, lw=space_size)
    pyplot.plot([3736,3766],[y,y], c=SU, lw=SU_size)
    pyplot.plot([3766,4151],[y,y], c=ITS, lw=space_size)
    pyplot.plot([4151,8120],[y,y], c=SU, lw=SU_size)

def PlotHistoneGene(y, color='black', size=4):
    pyplot.plot([0,5024],[y,y], c='grey', lw=size)
    pyplot.plot([129,930],[y,y], c='black', lw=size)
    pyplot.plot([1932,2309],[y,y], c='black', lw=size)
    pyplot.plot([1336,1708],[y,y], c='black', lw=size)
    pyplot.plot([2783,3095],[y,y], c='black', lw=size)
    pyplot.plot([3391,3802],[y,y], c='black', lw=size)


def PlotFst(fsts, labels, width):
    fst_list=[]
    for f in fsts:
        nan_ind=numpy.isnan(f)
        fst_list.append(f[~nan_ind])
    pyplot.boxplot(fst_list)
    pyplot.xticks(range(1, len(labels)+1), labels, size=16)
    pyplot.ylabel('Fst', size=20)
    pyplot.yticks(size=16)

    for i in range(len(fst_list)):

        jitter=numpy.random.uniform(-1*width,width,size=len(fst_list[i]))
        pyplot.scatter(jitter+i+1, fst_list[i],s=10, c='red',edgecolor=None, alpha=.6)
    pyplot.show()





def PlotCDF(data, steps, label='', c='red', log=True):
    x=numpy.copy(data)
    if log==True: x=numpy.log10(x)
    is_inf=numpy.isinf(x)
    x=x[~is_inf]
    cdf=[]
    step_size=1./steps
    if log==False: domain= numpy.linspace(min(x), max(x), steps)
    if log==True:domain= numpy.log10( numpy.logspace(min(x), max(x), steps))
    for s in domain:

        cdf.append( scipy.stats.percentileofscore(x,  s))
    pyplot.plot(domain,numpy.array( cdf) /100., c=c, label=label)
##    pyplot.show()


def PlotVarianceBoxplot(samples,lin_fits, labels):

    variances=[]

    for i in range(len(samples)):
        variances.append(VarianceAtRelevantFeatures(samples[i],0, lin_fits[i]))
        variances.append(VarianceAtIrrelevantFeatures(samples[i],0, lin_fits[i]))
    bx_plt=pyplot.boxplot(variances, widths=.8)
    color_dict={0:'red', 1:'blue'}
    for i in range(len(variances)):
        jitter=numpy.random.uniform(-.1,.1,size=len(variances[i]))
        pyplot.scatter(jitter+i+1, variances[i], alpha=.3, c=color_dict[i%2])
    set_med_colors=[med.set_c('firebrick') for med in    bx_plt['medians']]
    set_med_colors=[med.set_linewidth(3) for med in    bx_plt['medians']]
    set_med_colors=[med.set_zorder(7) for med in    bx_plt['medians']]
    set_med_colors=[med.set_zorder(8) for med in    bx_plt['means']]
    set_med_colors=[med.set_c ('black') for med in    bx_plt['whiskers']]
    set_med_colors=[med.set_c ('black') for med in    bx_plt['caps']]
    set_med_colors=[med.set_linewidth(3) for med in    bx_plt['boxes']]
    set_med_colors= [bx_plt['boxes'][ med] .set_c('orange') for med in range(0,len(  bx_plt['boxes']),2)]
    set_med_colors= [bx_plt['boxes'][ med] .set_alpha(.7) for med in range(0,len(  bx_plt['boxes']),2)]
    set_med_colors= [bx_plt['medians'][ med] .set_c('orange') for med in range(0,len(  bx_plt['boxes']),2)]
    set_med_colors= [bx_plt['boxes'][ med] .set_c('dodgerblue') for med in range(1,len(  bx_plt['boxes']),2)]
    set_med_colors= [bx_plt['boxes'][ med] .set_alpha(.7) for med in range(1,len(  bx_plt['boxes']),2)]
    set_med_colors= [bx_plt['medians'][ med] .set_c('dodgerblue') for med in range(1,len(  bx_plt['boxes']),2)]
    pyplot.xticks(numpy.arange(0,len(labels))*2+1.5, labels, size=14)
    pyplot.ylabel('Standard Deviation', size=14)
    max_var=max([max(v) for v in variances])
    pyplot.ylim(-.01, max_var+.01)
    pyplot.yticks(size=14)
    pyplot.show()


def PlotFstBoxplot(samples,lin_fits, labels):

    variances=[]

    for i in range(len(samples)):
        variances.append(FstAtRelevantFeatures(samples[i],0, lin_fits[i]))
        variances.append(FstAtIrrelevantFeatures(samples[i],0, lin_fits[i]))
    bx_plt=pyplot.boxplot(variances, widths=.8)
    color_dict={0:'red', 1:'blue'}
    for i in range(len(variances)):
        jitter=numpy.random.uniform(-.4,.4,size=len(variances[i]))
        pyplot.scatter(jitter+i+1, variances[i], alpha=.3, c=color_dict[i%2])
    set_med_colors=[med.set_c('firebrick') for med in    bx_plt['medians']]
    set_med_colors=[med.set_linewidth(3) for med in    bx_plt['medians']]
    set_med_colors=[med.set_zorder(7) for med in    bx_plt['medians']]
    set_med_colors=[med.set_zorder(8) for med in    bx_plt['means']]
    set_med_colors=[med.set_c ('black') for med in    bx_plt['whiskers']]
    set_med_colors=[med.set_c ('black') for med in    bx_plt['caps']]
    set_med_colors=[med.set_linewidth(3) for med in    bx_plt['boxes']]
    set_med_colors= [bx_plt['boxes'][ med] .set_c('orange') for med in range(0,len(  bx_plt['boxes']),2)]
    set_med_colors= [bx_plt['boxes'][ med] .set_alpha(.7) for med in range(0,len(  bx_plt['boxes']),2)]
    set_med_colors= [bx_plt['medians'][ med] .set_c('orange') for med in range(0,len(  bx_plt['boxes']),2)]
    set_med_colors= [bx_plt['boxes'][ med] .set_c('dodgerblue') for med in range(1,len(  bx_plt['boxes']),2)]
    set_med_colors= [bx_plt['boxes'][ med] .set_alpha(.7) for med in range(1,len(  bx_plt['boxes']),2)]
    set_med_colors= [bx_plt['medians'][ med] .set_c('dodgerblue') for med in range(1,len(  bx_plt['boxes']),2)]
    pyplot.xticks(numpy.arange(0,len(labels))*2+1.5, labels, size=14)
    pyplot.ylabel('Standard Deviation', size=14)
    max_var=max([max(v) for v in variances])
    pyplot.ylim(-.01, max_var+.01)
    pyplot.yticks(size=14)
    pyplot.show()



def ComparePCAtoRegression(freq, pca, lin_fit):
    ind=numpy.array(good_ind)
    rescaled=AlleleRescale(freq[ind,:])
    for i in range(3):
        ax=pyplot.gca()
        ax.yaxis. grid(True)
        ax.xaxis. grid(True)
        predicted=lin_fit[i].predict(rescaled)
        score=lin_fit[i].score(rescaled,pca.components_[i,:]+numpy.mean( pca.components_[i,:]))
        pyplot.scatter (pca.components_[i,:] , predicted+numpy.mean( pca.components_[i,:]) , s=60, alpha=.7,c=colors)
        pyplot.scatter([],[], s=0,label='$r^2$={:.4f}'.format(score))
        pyplot.legend(loc=9, fontsize=20)
        pyplot.yticks(size=14)
        pyplot.xticks(size=14)
        pyplot.ylabel('Predicted Loadings', size=20 )
        pyplot.xlabel('PC{0} Loadings'.format(i+1), size=20)
        ax.set_aspect('equal')
        pyplot.tight_layout()
        pyplot.show()


#-----------------------------------------------------------------------------
#Identify real variant positions

def SimulateVariants(readcounts, prob=1e-3):
    return scipy.stats.binom.rvs(readcounts, prob)

def ChiTestRepeat(table, error_rate=1e-3):
    read_count=table.sum(1)
    var_counts=read_count-numpy.max(table,1)
    min_cvg= numpy.min(read_count,0)
    ind=numpy.where(min_cvg>=10)[0]

    p_list=[]
    for i in ind:

        chi2, p=scipy.stats.chisquare( var_counts[:,i], read_count[:,i].astype(float)*error_rate)

        p_list.append(p)
    return numpy.array( p_list)

def MonteCarloTestRepeat(table,pos, reps=100, error_rate=1e-3, verbose=False):
    samples=table.shape[0]
    read_count=table[:,:,:] .sum(1)
    print read_count.shape
    read_count=read_count[:,pos]

    major_allele_prop,major_allele, all_pos=GetMajorAlleleFreq(table)

    var_counts=(read_count*(1-major_allele_prop)).astype(int)
    read_count=read_count.astype(int)
    min_cvg= numpy.min(read_count,0)


    p_list=[]
    for i in range(read_count.shape[1]):
        P_obs=numpy.log( (1.- scipy.stats.binom.cdf(var_counts[:,i], read_count[:,i], error_rate))).sum()

        more_extreme=0.
        monte_carlo_sample=scipy.stats.binom.rvs(read_count[:,i], error_rate, size=( reps,samples))

        logcdf=numpy.log (1.-scipy.stats.binom.cdf(monte_carlo_sample, read_count[:,i], error_rate)).sum(1)
        if verbose==True:
            print monte_carlo_sample.shape
            print numpy.log( (1.- scipy.stats.binom.cdf(var_counts[:,i], read_count[:,i], error_rate)))
            print numpy.log (1.-scipy.stats.binom.cdf(monte_carlo_sample, read_count[:,i], error_rate))[0,:]
            print logcdf.shape
            print i, P_obs
            print logcdf
        p_list.append(float((P_obs<logcdf).sum())/reps)
##        print jabber
##        p_list.append(p)
    return numpy.array( p_list)


def MonteCarloOutlierTestRepeat(table,pos, reps=100, error_rate=1e-3, verbose=False):
    samples=table.shape[0]
    read_count=table[:,:,:] .sum(1)
    print read_count.shape
    read_count=read_count[:,pos]

    major_allele_prop,major_allele, all_pos=GetMajorAlleleFreq(table)

    var_counts=(read_count*(1-major_allele_prop)).astype(int)
    read_count=read_count.astype(int)
    min_cvg= numpy.min(read_count,0)


    p_list=[]
    for i in range(read_count.shape[1]):
        P_obs=numpy.max( abs(numpy.mean( var_counts[:,i])-var_counts[:,i] ))

        more_extreme=0.
        monte_carlo_sample=scipy.stats.binom.rvs(read_count[:,i], error_rate, size=( reps,samples))

        logcdf=numpy.max( abs( numpy.mean(monte_carlo_sample,1)[:,None]-monte_carlo_sample),1)
        if verbose==True:
            print monte_carlo_sample.shape
            print numpy.log( (1.- scipy.stats.binom.cdf(var_counts[:,i], read_count[:,i], error_rate)))
            print numpy.log (1.-scipy.stats.binom.cdf(monte_carlo_sample, read_count[:,i], error_rate))[0,:]
            print logcdf.shape
            print i, P_obs
            print logcdf
        p_list.append(float((P_obs >logcdf).sum())/reps)
##        print jabber
##        p_list.append(p)
    return numpy.array( p_list)


def PagesTest(obs, exp):
    pass


def main(argv):
    param={}
    print argv
    for i in range(1, len(argv), 2):
        param[argv[i]]= argv[i+1]
    print param

    if param.has_key('-i')==True:
##        snps=numpy.load(param['-i'])
        snps=pandas.read_csv(param['-i'], delimiter=',')
        print snps.shape
        snps[snps=='A/A']=1.
        snps[snps=='A/B']=.5
        snps[snps=='B/B']=0.
        snps[snps=='']=numpy.nan
        snps=snps[1:][1:].as_matrix()
        imputer=sklearn.preprocessing.Imputer()
        snps=imputer.fit_transform(snps[1:,1:])
##        snp_freq, snp_pos=GetFreqArray(snps, 10)
        pca, lin_fit, data=PlotPCA(snps,[], True)
        PlotHistone(snps, numpy.arange(0,28501), lin_fit, colors, samples)

if __name__ == '__main__':
    main(sys.argv)