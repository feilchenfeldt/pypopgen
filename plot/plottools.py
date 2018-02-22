import io, re
import itertools
import matplotlib as mpl
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

mpl.rcParams['svg.fonttype'] = 'none'

#fix svg inkscape bug

def fixmiterlimit(svgdata, miterlimit = 10):
    # miterlimit variable sets the desired miterlimit
    mlfound = False
    svgout = ""
    for line in svgdata:
        if not mlfound:
            # searches the stroke-miterlimit within the current line and changes its value
            mlstring = re.subn(r'stroke-miterlimit:([0-9]+)', "stroke-miterlimit:" + str(miterlimit), line)
        #if mlstring[1]: # use number of changes made to the line to check whether anything was found
            #mlfound = True
            svgout += mlstring[0] + '\n'
        else:
            svgout += line + '\n'
    return svgout



def svg_save(fn,**kwa):
    imgdata = io.StringIO() # initiate StringIO to write figure data to
    # the same you would use to save your figure to svg, but instead of filename use StringIO object
    ax = plt.gca()
    plt.savefig(imgdata,
                 format='svg', bbox_inches='tight',**kwa)
    imgdata.seek(0)  # rewind the data
    svg_dta = imgdata.getvalue()  # this is svg data
    svgoutdata = fixmiterlimit(re.split(r'\n', svg_dta)) # pass as an array of lines
    svgfh = open(fn, 'w')
    svgfh.write(svgoutdata.encode('utf-8').strip())
    svgfh.close()
    return ax



def multiscatter(x_list, y_list, color_list, label_list, ax=None, **kwa):
    if ax is None:
        fig = plt.figure(figsize=(6,4))
        ax = plt.gca()
    data_l = []
    legend_elements = []
    for xs, ys, c, lbl in itertools.izip(x_list, y_list, color_list, label_list):
        d = pd.DataFrame({'d':ys,'c':c}, index=xs)
        data_l.append(d)
        legend_element = mpl.lines.Line2D([0], [0], marker='o', color='w', label=lbl,
                                          markerfacecolor=c, markersize=7)
        legend_elements.append(legend_element)
    #Thre driopna is crucial, otherwise color will not match values
    data =  pd.concat(data_l).sort_index().dropna()
    ax.scatter(data.index.values,data['d'].values, color=list(data['c'].values),**kwa)
    l = plt.legend(handles=legend_elements)
    return ax, l
