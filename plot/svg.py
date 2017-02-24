#fix svg inkscape bug
import io, re


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
