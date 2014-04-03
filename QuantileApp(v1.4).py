import time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from scipy.stats.stats import pearsonr, spearmanr, gmean, ttest_1samp

def clear_terminal():
	import os
	os.system('cls' if os.name=='nt' else 'clear')	

def start_screen():
	print '''
===========================================
	    Welcome to the QuantileApp
    	       Fatty Acids

        	   Version 1.4
         Last updated 01.20.2014

===========================================
Written By: Pho Diep (phodiep@gmail.com)

          Written in Python 2.7.3

-------------------------------------------'''

def getPeakList(data):
	tempList = list()
	tempList += {'p':data['p1'],'q':data['q1'],'peakName':'14:0','FAtype':'Saturated','common':'Myristic'},
	tempList += {'p':data['p2'],'q':data['q2'],'peakName':'14:1n5','FAtype':'Monounsaturated','common':'<Fix this - unknown>'},
	tempList += {'p':data['p3'],'q':data['q3'],'peakName':'15:0','FAtype':'Saturated','common':'Pentadecylic'},
	tempList += {'p':data['p4'],'q':data['q4'],'peakName':'16:0','FAtype':'Saturated','common':'Palmitic'},
	tempList += {'p':data['p5'],'q':data['q5'],'peakName':'16:1n9t','FAtype':'Trans','common':'7 trans heyadecenoic'},
	tempList += {'p':data['p6'],'q':data['q6'],'peakName':'16:1n7t','FAtype':'Trans','common':'Palmitelaidic'},
	tempList += {'p':data['p7'],'q':data['q7'],'peakName':'16:1n9c','FAtype':'Monounsaturated','common':'7-hexadecenoic'},
	tempList += {'p':data['p8'],'q':data['q8'],'peakName':'16:1n7c','FAtype':'Monounsaturated','common':'Palmitoleic'},
	tempList += {'p':data['p9'],'q':data['q9'],'peakName':'17:0','FAtype':'Saturated','common':'Margaric'},
	tempList += {'p':data['p10'],'q':data['q10'],'peakName':'U1','FAtype':'unknown','common':'unknown'},
	tempList += {'p':data['p11'],'q':data['q11'],'peakName':'17:1n9c','FAtype':'Monounsaturated','common':'heptadecenoic'},
	tempList += {'p':data['p12'],'q':data['q12'],'peakName':'18:0', 'FAtype':'Saturated','common':'Stearic'},
	tempList += {'p':data['p13'],'q':data['q13'],'peakName':'18:1n10-12t','FAtype':'Trans','common':'transoctadecenoic'},
	tempList += {'p':data['p14'],'q':data['q14'],'peakName':'18:1n9t','FAtype':'Trans','common':'Elaidic'},
	tempList += {'p':data['p15'],'q':data['q15'],'peakName':'18:1n8t','FAtype':'Trans','common':'transoctadecenoic'},
	tempList += {'p':data['p16'],'q':data['q16'],'peakName':'18:1n7t','FAtype':'Trans','common':'transvaccenic'},
	tempList += {'p':data['p17'],'q':data['q17'],'peakName':'18:1n6t','FAtype':'Trans','common':'transoctadecenoic'},
	tempList += {'p':data['p18'],'q':data['q18'],'peakName':'18:1n8c','FAtype':'Monounsaturated','common':'10-octadecenoic'},
	tempList += {'p':data['p19'],'q':data['q19'],'peakName':'18:1n9c','FAtype':'Monounsaturated','common':'Oleic'},
	tempList += {'p':data['p20'],'q':data['q20'],'peakName':'18:1n7c','FAtype':'Monounsaturated','common':'cis-vaccenic'},
	tempList += {'p':data['p21'],'q':data['q21'],'peakName':'18:1n5c','FAtype':'Monounsaturated','common':'13-octadecenoic'},
	tempList += {'p':data['p22'],'q':data['q22'],'peakName':'18:2n6tt','FAtype':'Trans','common':'6-neolaiolic'},
	tempList += {'p':data['p23'],'q':data['q23'],'peakName':'U2','FAtype':'unknown','common':'unknown'},
	tempList += {'p':data['p24'],'q':data['q24'],'peakName':'18:2n6ct','FAtype':'Trans','common':'cistrans linoelaiolic'},
	tempList += {'p':data['p25'],'q':data['q25'],'peakName':'18:2n6tc','FAtype':'Trans','common':'transcis linoelaiolic'},
	tempList += {'p':data['p26'],'q':data['q26'],'peakName':'18:2n6','FAtype':'Omega-6','common':'Linoleic'},
	tempList += {'p':data['p27'],'q':data['q27'],'peakName':'20:0','FAtype':'Saturated','common':'Arachidic'},
	tempList += {'p':data['p28'],'q':data['q28'],'peakName':'18:3n6','FAtype':'Omega-6','common':'Gamma-linolenic'},
	tempList += {'p':data['p29'],'q':data['q29'],'peakName':'20:1n9','FAtype':'Monounsaturated','common':'Gondoic'},
	tempList += {'p':data['p30'],'q':data['q30'],'peakName':'18:3n3','FAtype':'Omega-3','common':'alpha-Linolenic'},
	tempList += {'p':data['p31'],'q':data['q31'],'peakName':'20:2n6','FAtype':'Omega-6','common':'Eicosadienoic'},
	tempList += {'p':data['p32'],'q':data['q32'],'peakName':'22:0','FAtype':'Saturated','common':'Behenic'},
	tempList += {'p':data['p33'],'q':data['q33'],'peakName':'20:3n6','FAtype':'Omega-6','common':'Dihomo-gamma-linolenic'},
	tempList += {'p':data['p34'],'q':data['q34'],'peakName':'22:1n9','FAtype':'Monounsaturated','common':'Erucic'},
	tempList += {'p':data['p35'],'q':data['q35'],'peakName':'20:3n3','FAtype':'Omega-3','common':'EicoSaturatedrienoic'},
	tempList += {'p':data['p36'],'q':data['q36'],'peakName':'20:4n6','FAtype':'Omega-6','common':'Arachidonic'},
	tempList += {'p':data['p37'],'q':data['q37'],'peakName':'23:0','FAtype':'Saturated','common':'Tricosylic'},
	tempList += {'p':data['p38'],'q':data['q38'],'peakName':'22:2n6','FAtype':'Omega-6','common':'Docosadienoic'},
	tempList += {'p':data['p39'],'q':data['q39'],'peakName':'24:0','FAtype':'Saturated','common':'Lignoceric'},
	tempList += {'p':data['p40'],'q':data['q40'],'peakName':'20:5n3','FAtype':'Omega-3','common':'Eicosapentaenoic'},
	tempList += {'p':data['p41'],'q':data['q41'],'peakName':'24:1n9','FAtype':'Monounsaturated','common':'Nervonic'},
	tempList += {'p':data['p42'],'q':data['q42'],'peakName':'22:4n6','FAtype':'Omega-6','common':'Adrenic'},
	tempList += {'p':data['p43'],'q':data['q43'],'peakName':'22:5n6','FAtype':'Omega-6','common':'Docosapentaenoic'},
	tempList += {'p':data['p44'],'q':data['q44'],'peakName':'U5','FAtype':'unknown','common':'unknown'},
	tempList += {'p':data['p45'],'q':data['q45'],'peakName':'22:5n3','FAtype':'Omega-3','common':'Docosapentaenoic'},
	tempList += {'p':data['p46'],'q':data['q46'],'peakName':'22:6n3','FAtype':'Omega-3','common':'Docosahexaenoic'},
	return tempList

def add_ScatterPlot(dataX,dataY,fig,subR,subC,subN,title):
	
	fit = np.polyfit(dataX,dataY,1)	#calculate trendline
	fit_fn = np.poly1d(fit)

	ax = fig.add_subplot(subR,subC,subN)

	ax.scatter(dataX,dataY,color='b', marker='.') #add scatter plot
	ax.plot(dataX,fit_fn(dataX),color='r',linewidth=1.0) #add trendline in red
	
	plt.title(title, fontsize = 12) #add plot title

	ax.locator_params(nbins=4)
	plt.setp(ax.get_xticklabels(), fontsize=6)
	plt.setp(ax.get_yticklabels(), fontsize=6)

	ax.set_xlabel('%', fontsize=10)
	ax.set_ylabel('Abs', fontsize=10)

	return fig

def scatterPlot(data,tempTitle):
	fig = plt.figure(figsize=(12,9), dpi=100)
	fig.suptitle(tempTitle + ' (n = '+str(len(data))+')')
	peakList = getPeakList(data)
	countLocation = 0
	for entry in peakList:
		countLocation += 1
		try:
			add_ScatterPlot(entry['p'],entry['q'],fig,7,7,countLocation,entry['peakName'])
		except: pass

	plt.tight_layout()
	plt.subplots_adjust(top=0.92)
	return plt

def add_HistPlot(dataX,fig,subR,subC,subN,title):
	ax = fig.add_subplot(subR,subC,subN)

	ax.hist(dataX, bins =10) #add hist plot

	plt.title(title, fontsize = 12) #add plot title
	
	plt.setp(ax.get_xticklabels(), fontsize=6, rotation=90)
	plt.setp(ax.get_yticklabels(), fontsize=6)

	ax.set_xlabel('Abs', fontsize=10)
	ax.set_ylabel('Count', fontsize=10)

	return fig


def histPlot(data,tempTitle):
	fig = plt.figure(figsize=(12,9), dpi=100)
	fig.suptitle(tempTitle + ' (n = '+str(len(data))+')')
	peakList = getPeakList(data)
	countLocation = 0
	for entry in peakList:
		countLocation += 1
		try:
			add_HistPlot(list(entry['q']),fig,7,7,countLocation,entry['peakName'])
		except: pass

	plt.tight_layout()
	plt.subplots_adjust(top=0.92)
	return plt


def quantile(column,quantile=5):
	# categorizes each entry into quantile bin
	try:
		q = pd.qcut(column, quantile)
		return q.labels + 1
	except:
		return 'NaN'

def apply_quantile(data,bins):
	# reads csv raw data, applies quantile to data
	return data.apply(quantile,quantile=bins)

def get_buckets(bins):
	tempDict = dict()
	
	for row in range(1,bins+1,1):
		for col in range(1,bins+1,1):
			tempDict[str(row)+str(col)] = ''

	return tempDict

def get_labels(bins): 
	#columns 1...x    rows x...1
	#   1 2 3 4
	# 4
	# 3
	# 2
	# 1
	# return list(range(1,bins+1,1)), list(range(bins,0,-1)) 
	
	#columns 1...x    rows 1...x
	#   1 2 3 4
	# 1
	# 2
	# 3
	# 4
	return list(range(1,bins+1,1)), list(range(1,bins+1,1))

def make_QuantileSummary(dataX,dataY,tempDict):
	# tempDict = get_buckets(bins)

	for rowX, rowY in zip(dataX, dataY):
		try:
			if tempDict[str(rowY)+str(rowX)] == '':
				tempDict[str(rowY)+str(rowX)] = 1
			else:
				tempDict[str(rowY)+str(rowX)] += 1
		except: pass

	return tempDict

def make_table(row_labels,col_labels,tempDict):
	table_vals = list()
	for row in row_labels:
		temp_vals = list()
		for col in col_labels:
			temp_vals += tempDict[str(row)+str(col)], #pulls values from '1x...11' 
		table_vals += [temp_vals]					#add row to table
	return table_vals

def calc_QuantileSummary(dataX,dataY,fig,subR,subC,subN,title,bins):
	# creates a summary of each category

	tempDict = make_QuantileSummary(dataX,dataY,get_buckets(bins))
	col_labels, row_labels = get_labels(bins) #1...x    x...1
	table_vals = make_table(row_labels,col_labels,tempDict)


	ax = fig.add_subplot(subR,subC,subN)
	plt.title(title) #add plot title
	ax.set_frame_on(False)
	ax.set_xlabel('%')
	ax.set_ylabel('Abs')
	ax.get_xaxis().set_visible(False)
	ax.get_yaxis().set_visible(False)
	
	sum_table=ax.table(cellText=table_vals,rowLabels=row_labels,colLabels=col_labels,loc='center')
	sum_table.set_fontsize(10)
	return fig

def quantile_Summary(data,tempTitle,bins):
	fig = plt.figure(figsize=(12,9), dpi=100)
	fig.suptitle(tempTitle + ' (n = '+str(len(data.index))+')')

	peakList = getPeakList(data)
	countLocation = 0
	for entry in peakList:
		countLocation += 1
		try:
			calc_QuantileSummary(entry['p'],entry['q'],fig,7,7,countLocation,entry['peakName'],bins)
		except: pass

	plt.tight_layout()
	plt.subplots_adjust(top=0.92)
	return plt

def apply_stats(data,runTTest):
	peakList = getPeakList(data)
	tempList = list()

	colNames = ['Fatty Acid Type',           #1
				'Peak Name',                 #2
				'Pearson Coefficient',       #3
				'Pearson P Value',           #4
				'Spearman Coefficient',      #5
				'Spearman P Value',          #6
				'P Geometric Mean (%)',      #7
				'Q Geometric Mean (ug/ml)',  #8
				'P Mean (%)',                #9
				'P Stdev',                   #10
				'Q Mean (ug/ml)',            #11
				'Q Stdev',                   #12
				'P T-test',                  #13
				'P T-test P value',          #14
				'Q T-test',                  #15
				'Q T-test P value',          #16
				'Common Name']               #17


	for entry in peakList:
		try:
			pearson = pearsonr(entry['p'],entry['q'])
			spearman = spearmanr(entry['p'],entry['q'])
			if runTTest == 'y':
				ttestP = ttest_1samp(entry['p'],0)
				ttestQ = ttest_1samp(entry['q'],0)
			else:
				ttestP = ('-','-')
				ttestQ = ('-','-')
			tempList += [entry['FAtype'],            #1
						entry['peakName'],           #2
						pearson[0],                  #3
						pearson[1],                  #4
						spearman[0],                 #5
						spearman[1],                 #6
						gmean(entry['p']),           #7
						gmean(entry['q']),           #8
						np.mean(entry['p']),         #9
						np.std(entry['p'],ddof=1),   #10
						np.mean(entry['q']),         #11
						np.std(entry['q'],ddof=1),   #12
						ttestP[0],                   #13
						ttestP[1],                   #14
						ttestQ[0],                   #15
						ttestQ[1],                   #16
						entry['common']],            #17
		except: pass

	return pd.DataFrame(tempList, columns=colNames)




#------------MAIN------------
clear_terminal()
start_screen()

tempDataFile = raw_input('\nEnter the file to be processed (default:Data.csv): \n') or 'Data.csv'
tempName = raw_input('\nEnter Study Name for file export (default:test): \n') or 'test'
tempTitle = raw_input('\nEnter Description for title of report (default:test_title): \n') or 'test_title'
tempBins = int(raw_input('\nEnter number of bins (4=quartile, 5=quintile/default):') or 5) 

scatterPlot_run = raw_input('\nScatter Plot? y/n (default:y): ' ) or 'y'
histPlot_run = raw_input('\nHistogram Plot? y/n (default:y): ' ) or 'y'
quantPlot_run = raw_input('\nQuantile Plot? y/n (default:y): ' ) or 'y'
stats_run = raw_input('\nStatistics? y/n (default:y): ' ) or 'y'


MasterTime = time.time()
try:
	data = pd.read_csv(tempDataFile, index_col=0) #import csv as dataframe

	try: #=========scatter plot Percentage vs AbsoluteQuant=======
		if scatterPlot_run == 'y':
			startTime = time.time()
			scatter = scatterPlot(data,tempTitle)
			scatter.savefig('%s_Results_ScatterPlot.jpeg' % (tempName,))
			print '\nScatter plot successfully printed in %s seconds' % str(time.time() - startTime)
	except: print '\n...Scatter plot could not be printed...'

	try: #=========Histogram plot AbsoluteQuant=======================
		if histPlot_run == 'y':
			startTime = time.time()
			histo = histPlot(data,tempTitle)
			histo.savefig('%s_Results_HistPlot.jpeg' % (tempName,))
			print '\nHist plot successfully printed in %s seconds' % str(time.time() - startTime)
	except: print '\n...Hist plot could not be printed...'

	try: #=========Quantile (5-bin) summary=======================
		if quantPlot_run == 'y':
			startTime = time.time()
			dataQuant = apply_quantile(data,tempBins)
			dataQuant.to_csv('%s_Results_Quantile.csv' % (tempName,))
			summary = quantile_Summary(dataQuant,tempTitle,tempBins)
			summary.savefig('%s_Results_QuantileSummary.jpeg' % (tempName,))
			print '\nQuantile plot successfully printed in %s seconds' % str(time.time() - startTime)
	except: print '\n...Quantile summary could not be printed...'


	try: #=========Stats=======================
		if stats_run == 'y':
			ttest_run = raw_input('\n1 sample T-Test? y/n (default:y): ' ) or 'y'
			startTime = time.time()
			dataStat = apply_stats(data,ttest_run)
			dataStat.to_csv('%s_Results_Statistics.csv' % (tempName,),index=False)
			print '\nStatistics summary successfully printed in %s seconds' % str(time.time() - startTime)
	except: print '\n...Statistics summary could not be printed...'
	


	print '\nThe data has been successfully processed.'

except:
	print '''
\nSorry, the File could not be processed... 
Be sure the correct file name was entered
and the file has been saved in the correct location'''

print '\nTotal time: ' + str(time.time() - MasterTime) + ' seconds\n\n'