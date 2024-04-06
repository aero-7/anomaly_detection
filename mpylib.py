import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
#from matplotlib.pylab import rcParams
#from enum import Enum
from datetime import datetime, timedelta
from colorama import Fore


rc = mpl.rc
rc('font', size=24, weight='medium') #22
rc('figure', dpi=100, facecolor='whitesmoke', titlesize=27, titleweight='medium') #30--regular
rc('axes', titlecolor='k', titlesize=25, titleweight='regular', titlepad=15, labelpad=15)
#rc('axes', labelcolor='r', labelsize=10, labelweight='medium')
rc('axes', facecolor='#e6e6e6', labelcolor='r', grid=True) #'gainsboro' #'whitesmoke' #
rc('axes.spines', top=False, right=False, bottom=False, left=False)
rc('grid', color='white', linestyle='-', linewidth=1)
rc('legend', fontsize=17, facecolor='#ffffff50')
rc('xtick', color='dimgrey', labelsize=22, labelcolor='dimgrey')
rc('xtick.major', width=2)
rc('ytick', color='dimgrey', labelsize=22, labelcolor='dimgrey')
rc('ytick.major', width=2)

def formatPlotAx(ax:plt.Axes):
    ax.tick_params(length=5, rotation=0)

def getPlotAx()->plt.Axes:
    ax = plt.figure(figsize=(21,7)).add_subplot(111)
    formatPlotAx(ax=ax)
    return ax
    
"""--kp
def adjustPlotView(ax:plt.Axes=None, xMin:int|float|str=None, xMax:int|float|str=None, yMin:int|float|str=None, yMax:int|float|str=None):
    if ax==None:
        print(Fore.RED+"\nError! "+Fore.GREEN+"You have not provided the axes object of the plot you want to adjust.\n") #Error Message
    else:
        if type(xMin)==str: 
            try:
                xMin = datetime.strptime(xMin, '%Y-%m-%d %H:%M:%S')
            except:
                xMin=None
                print(Fore.RED+"\nError! "+Fore.GREEN+f"The xMin argument provided is invalid. It should be a date-time string of the form 'YYYY-MM-DD hh:mm:ss' (for a date-time horizontal axis) or a number of type 'int' or 'float'.\n") #Error Message
        elif type(xMin)==int or type(xMin)==float:
            xMin=float(xMin)
        #elif xMin==None:
            #xMin = ax.get_xlim()[0]
        elif xMin != None:
            xMin=None
            print(Fore.RED+"\nError! "+Fore.GREEN+f"The xMin argument provided is invalid. It should either be a date-time string of the form 'YYYY-MM-DD hh:mm:ss' (for a date-time horizontal axis) or a number of type 'int' or 'float'.\n") #Error Message
        
        if type(xMax)==str: 
            try:
                xMax = datetime.strptime(xMax, '%Y-%m-%d %H:%M:%S')
            except:
                xMax=None
                print(Fore.RED+"\nError! "+Fore.GREEN+f"The xMax argument provided is invalid. It should be a date-time string of the form 'YYYY-MM-DD hh:mm:ss' (for a date-time horizontal axis) or a number of type 'int' or 'float'.\n") #Error Message
        elif type(xMax)==int or type(xMax)==float:
            xMax=float(xMax)
        #elif xMax==None:
            #xMax = ax.get_xlim()[1]
        elif xMax != None:
            xMax=None
            print(Fore.RED+"\nError! "+Fore.GREEN+f"The xMax argument provided is invalid. It should be a date-time string of the form 'YYYY-MM-DD hh:mm:ss' (for a date-time horizontal axis) or a numbmer of type 'int' or 'float'.\n") #Error Message
        
        if (type(xMin) != type(xMax)) and not(xMin==None or xMax==None):
            xMin=None; xMax=None
            print(Fore.RED+"\nError! "+Fore.GREEN+f"The xMin and xMax arguments provided are invalid. They should all be either date-time strings of the form 'YYYY-MM-DD hh:mm:ss' (for a date-time horizontal axis) or numbers of type 'int' or 'float'.\n") #Error Message

        if type(yMin)==int or type(yMin)==float:
            yMin=float(yMin)
        #elif yMin==None:
            #yMin = ax.get_ylim()[0]
        elif yMin != None:
            yMin=None
            print(Fore.RED+"\nError! "+Fore.GREEN+f"Plot could not be adjusted vertically due to invalid yMin argument provided. It should be a number.\n") #Error Message

        if type(yMax)==int or type(yMax)==float:
            yMax=float(yMax)
        #elif yMax==None:
            #yMax = ax.get_ylim()[1]
        elif yMax != None:
            yMax=None
            print(Fore.RED+"\nError! "+Fore.GREEN+f"Plot could not be adjusted vertically due to invalid yMax argument provided. It should be a number.\n") #Error Message

        if (type(yMin) != type(yMax)) and not(yMin==None or yMax==None):
            yMin=None; yMax=None
            print(Fore.RED+"\nError! "+Fore.GREEN+f"The yMin and yMax arguments provided are invalid. They should be numbers.\n") #Error Message



        ax_xMn, ax_xMx = ax.get_xlim(); ax_dx = 0.1*(ax_xMx-ax_xMn)
        ax_yMn, ax_yMx = ax.get_ylim(); ax_dy = 0.1*(ax_yMx-ax_yMn)
        if xMin != None and xMax != None:
            if xMin<xMax: 
                ax.set_xlim(xmin=xMin,xmax=xMax)
            else:
                print(Fore.RED+"\nError! "+Fore.GREEN+"Plot could not be adjusted horizontally due to invalid x-limits. Provide a min value that is less than the max value.\n")    #Error Message
        
        elif xMin==None and xMax != None and xMax<=ax_xMn: #new
            ax.set_xlim(xmin=ax_xMn-ax_dx, xmax=ax_xMn)
        elif xMax==None and xMin != None and xMin>=ax_xMx:
            ax.set_xlim(xmin=ax_xMx, xmax=xMin+ax_dx)
        
        elif xMin != None:
            ax.set_xlim(xmin=xMin)
        elif xMax != None:
            ax.set_xlim(xmax=xMax)

        if yMin != None and yMax != None:
            if yMin<yMax: 
                ax.set_ylim(ymin=yMin,ymax=yMax)
            else:
                print(Fore.RED+"\nError! "+Fore.GREEN+"Plot could not be adjusted vertically due to invalid y-limits. Provide a min value that is less than the max value.\n")    #Error Message

        elif yMin==None and yMax != None and yMax<=ax_yMn: #new
            ax.set_ylim(ymin=ax_yMn-ax_dy, ymax=ax_yMn)
        elif yMax==None and yMin != None and yMin>=ax_yMx:
            ax.set_ylim(ymin=ax_yMx, ymax=yMin+ax_dy)

        elif yMin != None:
            ax.set_ylim(ymin=yMin)
        elif yMax != None:
            ax.set_ylim(ymax=yMax)
"""    


def adjustPlotView(ax:plt.Axes=None, xMin:int|float|str=None, xMax:int|float|str=None, yMin:int|float|str=None, yMax:int|float|str=None):
    if ax==None:
        print(Fore.RED+"\nError! "+Fore.GREEN+"You have not provided the axes object of the plot you want to adjust.\n") #Error Message
    else:
        xMn=None; xMx=None; yMn=None; yMx=None
        if type(xMin)==str: 
            try:
                xMn = datetime.strptime(xMin, '%Y-%m-%d %H:%M:%S')
            except:
                print(Fore.RED+"\nError! "+Fore.GREEN+f"The xMin argument provided is invalid. It should be a date-time string of the form 'YYYY-MM-DD hh:mm:ss' (for a date-time horizontal axis) or a number of type 'int' or 'float'.\n") #Error Message
        elif type(xMin)==int or type(xMin)==float:
            xMn=float(xMin)
        elif xMin != None:
            print(Fore.RED+"\nError! "+Fore.GREEN+f"The xMin argument provided is invalid. It should either be a date-time string of the form 'YYYY-MM-DD hh:mm:ss' (for a date-time horizontal axis) or a number of type 'int' or 'float'.\n") #Error Message
        
        if type(xMax)==str: 
            try:
                xMx = datetime.strptime(xMax, '%Y-%m-%d %H:%M:%S')
            except:
                print(Fore.RED+"\nError! "+Fore.GREEN+f"The xMax argument provided is invalid. It should be a date-time string of the form 'YYYY-MM-DD hh:mm:ss' (for a date-time horizontal axis) or a number of type 'int' or 'float'.\n") #Error Message
        elif type(xMax)==int or type(xMax)==float:
            xMx=float(xMax)
        elif xMax != None:
            print(Fore.RED+"\nError! "+Fore.GREEN+f"The xMax argument provided is invalid. It should be a date-time string of the form 'YYYY-MM-DD hh:mm:ss' (for a date-time horizontal axis) or a numbmer of type 'int' or 'float'.\n") #Error Message
        
        if (xMin != None) and (xMax != None) and (type(xMn) != type(xMx)):
            print(Fore.RED+"\nError! "+Fore.GREEN+f"The xMin and xMax arguments provided are invalid. They should all be either date-time strings of the form 'YYYY-MM-DD hh:mm:ss' (for a date-time horizontal axis) or numbers of type 'int' or 'float'.\n") #Error Message

        if type(yMin)==int or type(yMin)==float:
            yMn=float(yMin)
        elif yMin != None:
            print(Fore.RED+"\nError! "+Fore.GREEN+f"Plot could not be adjusted vertically due to invalid yMin argument provided. It should be a number.\n") #Error Message

        if type(yMax)==int or type(yMax)==float:
            yMx=float(yMax)
        elif yMax != None:
            print(Fore.RED+"\nError! "+Fore.GREEN+f"Plot could not be adjusted vertically due to invalid yMax argument provided. It should be a number.\n") #Error Message

        if (yMin != None) and (yMax != None) and (type(yMn) != type(yMx)):
            print(Fore.RED+"\nError! "+Fore.GREEN+f"The yMin and yMax arguments provided are invalid. They should be numbers of type 'int' or 'float'.\n") #Error Message

        ax_xMn, ax_xMx = float(ax.get_xlim()[0]), float(ax.get_xlim()[1]);  ax_dx = float(0.1*(ax_xMx-ax_xMn)) 
        ax_yMn, ax_yMx = float(ax.get_ylim()[0]), float(ax.get_ylim()[1]);  ax_dy = float(0.1*(ax_yMx-ax_yMn)) 
        if xMn != None and xMx != None:
            if xMn<xMx: 
                ax.set_xlim(xmin=xMn, xmax=xMx)
            else:
                print(Fore.RED+"\nError! "+Fore.GREEN+"Plot could not be adjusted horizontally due to invalid x-limits. Provide valid values for the parameters with xMin value less than xMax value.\n")    #Error Message
        elif xMin==None and (xMx != None) and ((type(xMx)==float and xMx<=ax_xMn) or (type(xMx)!=float and mpl.dates.date2num(xMx)<=ax_xMn)):# and mpl.dates.date2num(xMx)<=ax_xMn: #new
            if type(xMx)==float: 
                ax.set_xlim(xmin=xMx-ax_dx, xmax=xMx)
            else:   
                ax.set_xlim(xmin=xMx-mpl.dates.num2timedelta(ax_dx), xmax=xMx)
        elif xMax==None and (xMn != None) and ((type(xMn)==float and xMn>=ax_xMx) or (type(xMn)!=float and mpl.dates.date2num(xMn)>=ax_xMx)):
            if type(xMn)==float:  
                ax.set_xlim(xmin=xMn, xmax=xMn+ax_dx)
            else:
                ax.set_xlim(xmin=xMn, xmax=xMn+mpl.dates.num2timedelta(ax_dx))
        elif xMn != None:
            ax.set_xlim(xmin=xMn)
        elif xMx != None:
            ax.set_xlim(xmax=xMx)

        if yMin != None and yMax != None:
            if yMin<yMax: 
                ax.set_ylim(ymin=yMin,ymax=yMax)
            else:
                print(Fore.RED+"\nError! "+Fore.GREEN+"Plot could not be adjusted vertically due to invalid y-limits. Provide a min value that is less than the max value.\n")    #Error Message
        elif yMin==None and (yMx != None) and yMx<=ax_yMn: #new
            ax.set_ylim(ymin=ax_yMn-ax_dy, ymax=ax_yMn)
        elif yMax==None and (yMn != None) and yMn>=ax_yMx:
            ax.set_ylim(ymin=ax_yMx, ymax=ax_yMx+ax_dy)
        elif yMin != None:
            ax.set_ylim(ymin=yMin)
        elif yMax != None:
            ax.set_ylim(ymax=yMax)
    

    
def closestDateTimeTo(elapsed_time:int|float, in_dataframe:pd.DataFrame): #timeElapsed_in_minutes
    df = in_dataframe
    dtTmRslt = None
    if not(df.empty):
        if type(elapsed_time)==float or type(elapsed_time)==int:
            lowerLmts = df[df.Elapsed_Time <= elapsed_time].Elapsed_Time 
            upperLmts = df[df.Elapsed_Time >= elapsed_time].Elapsed_Time 
            timeVl = None 

            if len(lowerLmts)==0:
                timeVl = df.Elapsed_Time[0]
            elif len(upperLmts)==0:
                timeVl = df.Elapsed_Time[-1]
            else:
                lowerLmt = max(lowerLmts); upperLmt = min(upperLmts)

                if (elapsed_time-lowerLmt) < (upperLmt-elapsed_time):
                    timeVl = lowerLmt
                else:
                    timeVl = upperLmt
            dtTmRslt = df.index[df.Elapsed_Time==timeVl][0]
        else:
            print(Fore.RED+"\nError! "+Fore.GREEN+f"The times elapsed value should be a either an int or float type.\n") #Error Message
    else:
        print(Fore.RED+"\nError! "+Fore.GREEN+f"The dataset (in_dataframe argument) provided is empty. A pandas dataframe with a datetime ('YYYY-MM-DD hh:mm:ss') index and Elapsed_Time field (of type int or float) is expected.\n") #Error Message
    return dtTmRslt 

    
def closestElapsedTimeTo(dateTime_string:str, in_dataframe:pd.DataFrame):
    df = in_dataframe
    elpsdTmRslt = None
    if not(df.empty):
        if type(dateTime_string)==str:
            try:
                mDateTime = datetime.strptime(dateTime_string, '%Y-%m-%d %H:%M:%S')
                dtTm_1 = df.index[0]
                dtTm_2 = df.index[-1]
                if mDateTime < dtTm_1:
                    elpsdTmRslt = df.Elapsed_Time[0] 
                elif mDateTime > dtTm_2:
                    elpsdTmRslt = df.Elapsed_Time[-1] 
                else:
                    lowerLmt = (df[df.index <= mDateTime]).index[-1]
                    upperLmt = (df[df.index >= mDateTime]).index[0]
                    if (mDateTime-lowerLmt) < (upperLmt-mDateTime):
                        elpsdTmRslt = (df[df.index==lowerLmt]).Elapsed_Time[0]    
                    else:
                        elpsdTmRslt = (df[df.index==upperLmt]).Elapsed_Time[0]    
            except:
                print(Fore.RED+"\nError! "+Fore.GREEN+f"The date-time string provided is invalid. It should be in the format 'YYYY-MM-DD hh:mm:ss'.\n") #Error Message
        else:
            print(Fore.RED+"\nError! "+Fore.GREEN+f"The date-time value should be a string in the format 'YYYY-MM-DD hh:mm:ss'.\n") #Error Message
    else:
        print(Fore.RED+"\nError! "+Fore.GREEN+f"The dataset (in_dataframe argument) provided is empty. A pandas dataframe with a datetime ('YYYY-MM-DD hh:mm:ss') index and Elapsed_Time field (of type int or float) is expected.\n") #Error Message
    return elpsdTmRslt


def hilitePlots(on_ax:plt.Axes, for_data=pd.DataFrame(), from_elapsed_time=None, to_elapsed_time=None, is_dateTime_axis=False, show_time_windows=False, win_size=None, time_delta=None, min_y=None, max_y=None, hlt_clrs:list[str]=['lightsteelblue','aliceblue'], hlt_alpha:list[float]=[0.75,0.75], edge_clr='none'):        
    n_data = for_data 
    
    if from_elapsed_time != None: 
        n_data = n_data[n_data.Elapsed_Time >= from_elapsed_time]  
    if to_elapsed_time != None: 
        n_data = n_data[n_data.Elapsed_Time <= to_elapsed_time]

    if not(n_data.empty):   
        strtTm=None; endTm=None; 
        mrgn = time_delta/2    

        if is_dateTime_axis:
            mrgn = timedelta(minutes=float(mrgn)) 
            strtTm = n_data.index[0]
            endTm = n_data.index[-1]
        else:
            strtTm = n_data.Elapsed_Time[0]
            endTm = n_data.Elapsed_Time[-1]
        
        maxVl = max(for_data.Value.values)  
        minVl = min(for_data.Value.values) 
        
        if (min_y != None) and (min_y > minVl) and (min_y <= min(n_data.Value)):
            minVl = min_y
        if (max_y != None) and (max_y >= max(n_data.Value)) and (max_y < maxVl):
            maxVl = max_y
        
        htXtnd = 0.0075*(maxVl-minVl)

        if show_time_windows and not(win_size==None) and not(time_delta==None):
            clrStrs = hlt_clrs
            clr_alpha = hlt_alpha   
            clrStrs_sz = len(clrStrs)
            win_l_bound = strtTm     
            win_l_bound_limit = endTm   
            win_ndx = 0
            
            strt_pnt_x=None; wdth=None
            while win_l_bound <= win_l_bound_limit:                
                fclr_ndx = win_ndx % clrStrs_sz  
                strt_pnt_x = win_l_bound-mrgn
                wdth = win_size

                if is_dateTime_axis:
                    wdth = timedelta(minutes=float(win_size))

                on_ax.add_patch(
                    Rectangle(
                        (strt_pnt_x, minVl - htXtnd), wdth, (maxVl-minVl) + 2*htXtnd,
                        edgecolor='none',
                        lw=2,
                        facecolor=clrStrs[fclr_ndx],  
                        alpha=clr_alpha[fclr_ndx], 
                        fill=True,
                        zorder = -1
                    )
                )
                win_l_bound = win_l_bound + wdth 
                win_ndx += 1
        else:
            on_ax.add_patch(
                Rectangle(
                    (strtTm - mrgn, minVl - htXtnd), (endTm - strtTm + 2*mrgn), (maxVl - minVl + 2*htXtnd),
                    edgecolor=edge_clr,  
                    ls='--',
                    lw=2,
                    facecolor=hlt_clrs[0],  
                    alpha=hlt_alpha[0],
                    fill=True,
                    zorder = 1
                )
            )
    else:
        print(Fore.RED+"\nError! "+Fore.GREEN+f"You have no data in the date-time limits provided. 'to_elapsed_time' should be later than 'from_elapsed_time', in comformity to the provided dataset, in order to highlight a time period in the dataset.\n") #Error Message
