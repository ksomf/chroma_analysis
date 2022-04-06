import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as pl
import numpy as np

from hl import *

##-- TABLEUO COLOURS, GOOD SET OF COLOURS THAT ARE ALSO COLOUR BLIND FRIENDLY --##
tableuo_colours = [ (  31./255., 119./255., 180./255. )  ##-- MID BLUE     --##
                  , (  44./255., 160./255.,  44./255. )  ##-- GREEN        --##
                  , ( 214./255.,  39./255.,  40./255. )  ##-- RED          --##
                  , ( 148./255., 103./255., 189./255. )  ##-- PURPLE       --##
                  , ( 255./255., 127./255.,  14./255. )  ##-- ORANGE       --##
                  , ( 127./255., 127./255., 127./255. )  ##-- DARK GREY    --##
                  , ( 152./255., 223./255., 138./255. )  ##-- LIGHT GREEN  --##
                  , ( 158./255., 218./255., 229./255. )  ##-- GREY BLUE    --##
                  , (  23./255., 190./255., 207./255. )  ##-- LIGHT BLUE   --##
                  , ( 140./255.,  86./255.,  75./255. )  ##-- DARK BROWN   --##
                  , ( 196./255., 156./255., 148./255. )  ##-- LIGHT BROWN  --##
                  , ( 255./255., 152./255., 150./255. )  ##-- PINK ORANGE  --##
                  , ( 197./255., 176./255., 213./255. )  ##-- LIGHT PURPLE --##
                  , ( 174./255., 199./255., 232./255. )  ##-- GREY BLUE    --##
                  , ( 247./255., 182./255., 210./255. )  ##-- LIGHT PINK   --##
                  , ( 255./255., 187./255., 120./255. )  ##-- DARK ORANGE  --##
                  , ( 227./255., 119./255., 194./255. )  ##-- PINK         --##
                  , ( 188./255., 189./255.,  34./255. )  ##-- YELLOW       --##
                  , ( 219./255., 219./255., 141./255. )  ##-- TEAL         --##
                  , ( 199./255., 199./255., 199./255. )  ##-- GREY         --##
                  ]

for c in tableuo_colours:
  print "{:.3f} {:.3f} {:.3f}".format( c[0]*255, c[1]*255, c[2]*255  )

exit(1)

#bg=(95.0/100.0,95.0/100.0,95.0/100.0)
#bg=(249.0/255.0,249.0/255.0,249.0/255.0)
bg=(255.0/255.0,255.0/255.0,255.0/255.0)
  
red = ( 214./255.,  39./255.,  40./255. )  ##-- RED --##

def set_foregroundcolor(ax, color):
  for tl in ax.get_xticklines() + ax.get_yticklines():
    tl.set_color(color)
  for spine in ax.spines:
    ax.spines[spine].set_edgecolor(color)
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_color(color)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_color(color)
  ax.axes.xaxis.label.set_color(color)
  ax.axes.yaxis.label.set_color(color)
  ax.axes.xaxis.get_offset_text().set_color(color)
  ax.axes.yaxis.get_offset_text().set_color(color)
  ax.axes.title.set_color(color)
  lh = ax.get_legend()
  if lh != None:
    lh.get_title().set_color(color)
    lh.legendPatch.set_edgecolor('none')
    labels = lh.get_texts()
    for lab in labels:
      lab.set_color(color)
  for tl in ax.get_xticklabels():
    tl.set_color(color)
  for tl in ax.get_yticklabels():
    tl.set_color(color)

def rgba_to_rgb( background, source, alpha ):
  return ( (1-alpha)*background[0] + alpha*source[0]
         , (1-alpha)*background[1] + alpha*source[1]
         , (1-alpha)*background[2] + alpha*source[2] )


def plot_fill( ax, px, py, ps, c, lw, alpha=1.0, **kwargs ):
  if 'label' in kwargs:
    ax.plot( px, py, c=c, lw=lw, alpha=alpha, label=kwargs['label'] )
  else:
    ax.plot( px, py, c=c, lw=lw, alpha=alpha )
  if 'confidence_lower' in kwargs:
    lower, upper =  kwargs['confidence_lower'], kwargs['confidence_upper']
  else:
    lower, upper = py - ps, py + ps
  #ax.fill_between( px, lower, upper, facecolor=c, alpha=alpha*0.33, lw=0 )
  ax.fill_between( px, py - ps, py + ps, facecolor=rgba_to_rgb(bg,c,alpha*0.33), lw=0 )


def Plot( data, ax=False, **kwargs ):
  print dict(kwargs)
  lw = kwargs.get('lw',2)
  fs = kwargs.get('fs',14)
  ls = kwargs.get('ls',fs)
  colours = tableuo_colours*10
  if not ax:
    figsize = kwargs.get( 'figsize', (14,10.5) )
    plot = pl.figure( figsize=figsize )
    ax = pl.subplot( 111 )

    #ax.spines['top'].set_visible( False )
    #ax.spines['right'].set_visible( False )
    #ax.spines['bottom'].set_visible( False )
    #ax.spines['left'].set_visible( False )
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    pl.xticks( fontsize=fs )
    pl.yticks( fontsize=fs )
  pl.xscale( kwargs.get('xscale','linear') )
  pl.yscale( kwargs.get('yscale','linear') )
  #pl.rc( 'grid', linestyle=':', c='black' )
  #pl.rc( 'text', usetex=True )
  #pl.rc( 'font', family='serif' )

  
  for d in data:
    for k in ['x','y','s']:
      if k in d:
        d[k] = np.asarray(d[k])

    if 'y' in d:
      if np.any(np.iscomplex(d['y'])):
        if kwargs.get('data_type','real') == 'real':
          d['y'] = np.copy(d['y'].real)
          if 's' in d:
            d['s'] = np.copy(d['s'].real)
        else:
          d['y'] = np.copy(d['y'].imag)
          if 's' in d:
            d['s'] = np.copy(d['s'].imag)
      
  xs = np.asarray([ x for d in data for x in d.get('x',np.asarray([])).tolist() ])
  if len(xs.tolist()):
    xlims = MinMax( xs, scale=1.1 )
    xlims = kwargs.get('xlims',xlims)
    lowerx, upperx = xlims
    ys = []
    for d in data:
      if 'y' in d:
        mask = np.logical_and( lowerx <= d['x'], d['x'] <= upperx )
        if 's' in d:
          ys = ys + (d['y'][mask]-d['s'][mask]).tolist() + (d['y'][mask]+d['s'][mask]).tolist()
        else:
          ys = ys + d['y'][mask].tolist()
    ylims = MinMax(ys)
    if kwargs.get('includezero',False):
      ylims = [min(ylims[0],0),max(ylims[1],0)]
    ylims = kwargs.get('ylims',ylims)
    if 'ylimmod' in kwargs:
      ylims = list(ylims)
      ylims[0] = ylims[0] - kwargs['ylimmod']
      ylims[1] = ylims[1] + kwargs['ylimmod']
      ylims = tuple(ylims)
    kwargs['xlims'] = xlims
    kwargs['ylims'] = ylims
  else:
    xlims = kwargs['xlims']
    #ylims = kwargs['ylims']
  lowerx, upperx = xlims

  pl.xlabel( kwargs.get( 'xlabel', '' ), fontsize=ls )
  pl.ylabel( kwargs.get( 'ylabel', '' ), fontsize=ls )
  pl.title ( kwargs.get( 'title' , '' ), fontsize=ls )

  legend=0
  n = float(len(data))
  for i,d in enumerate(data):
    if 'c' in d:
      c = tableuo_colours[d['c']]
    else:
      c = colours.pop(0)
    args = {}
    if 'l' in d:
      args = { 'label' : d['l'] }
      legend = 1
    if 'y' in d and 'x' not in d:
      d['x'] = np.arange( d['y'].shape[0] )
    if 'displace' in kwargs and 'x' in d:
      amount=float(kwargs['displace'])
      if d['x'].shape[0] > 1:
        d['x'] = d['x'] + amount*(i/n)*(d['x'][1]-d['x'][0])
    if 'absdisplace' in kwargs and 'x' in d:
      d['x'] = d['x'] + i*float(kwargs['absdisplace'])
    if 'cl' in d:
      cl = d['cl']
      ax.text(cl[0],cl[1],cl[2],fontsize=fs,color=c,verticalalignment='center')
    if 'yl' in d:
      yl = d['yl']
      ax.text(yl[0],yl[1],yl[2],fonstize=fs,color=c,vecticalalignment='center',horizontalalignment='right')
    if 'y' in d:
      if 's' in d:
        if 'red_after' in kwargs:
          after = kwargs['red_after']
          ax.errorbar( d['x'][:after+1], d['y'][:after+1], d['s'][:after+1], lw=lw, ecolor=c  , mfc=(1.0,1.0,1.0), c=c  , mec=c  , fmt='o', mew=2.0, **args )
          ax.errorbar( d['x'][after+1:], d['y'][after+1:], d['s'][after+1:], lw=lw, ecolor=red, mfc=(1.0,1.0,1.0), c=red, mec=red, fmt='o', mew=2.0, **args )
        else:
          ax.errorbar( d['x'], d['y'], d['s'], lw=lw, ecolor=c, mfc=(1.0,1.0,1.0), c=c, mec=c, fmt='o', mew=2.0, **args )
      else:
        #print d.keys()
        ax.plot( d['x'], d['y'], d.get('fmt','o'), c=c, lw=lw, **args )
    #if 'b' in d:
      #print 'Bootstrap shape:', d['b'].shape
    if 'f' in d:
      fit = d['f']
      lower,upper = map(float, fit.get('range',xlims))
      if 'y' in d:
        args.pop('label','')
      if fit['type'] == 'poly':
        px   = np.linspace( lower, upper )
        pxl  = np.linspace( xlims[0], lower )
        pxu  = np.linspace( upper, xlims[1] )
        pyb  = boot_object( np.dot(fit['boot'], np.asarray([ px  ** coef for coef in fit['coefs'] ])) )
        pybl = boot_object( np.dot(fit['boot'], np.asarray([ pxl ** coef for coef in fit['coefs'] ])) )
        pybu = boot_object( np.dot(fit['boot'], np.asarray([ pxu ** coef for coef in fit['coefs'] ])) )
        plot_fill( ax, px , pyb['mean'] , pyb['std'] , c, lw            , **args )
        plot_fill( ax, pxl, pybl['mean'], pybl['std'], c, lw, alpha=0.33         )
        plot_fill( ax, pxu, pybu['mean'], pybu['std'], c, lw, alpha=0.33         )
      elif fit['type'] == 'eff_cosh':
        offset = fit['offset']
        px   = np.linspace( lower, upper )
        pxl  = np.linspace( xlims[0], lower )
        pxu  = np.linspace( upper, xlims[1] )
        pyb  = boot_object( -fit['boot'][:,1,None] * np.tanh( (px[None,:]-offset)  * fit['boot'][:,1,None] ) )
        pybl = boot_object( -fit['boot'][:,1,None] * np.tanh( (pxl[None,:]-offset) * fit['boot'][:,1,None] ) )
        pybu = boot_object( -fit['boot'][:,1,None] * np.tanh( (pxu[None,:]-offset) * fit['boot'][:,1,None] ) )
        plot_fill( ax, px,  pyb['mean'] , pyb['std'] , c, lw            , **args )
        plot_fill( ax, pxl, pybl['mean'], pybl['std'], c, lw, alpha=0.33         )
        plot_fill( ax, pxu, pybu['mean'], pybu['std'], c, lw, alpha=0.33         )
      elif fit['type'] == 'cosh':
        offset = fit['offset']
        px   = np.linspace( lower, upper )
        pxl  = np.linspace( xlims[0], lower )
        pxu  = np.linspace( upper, xlims[1] )
        pyb  = boot_object( fit['boot'][:,0,None] * np.cosh( (px[None,:]-offset)  * fit['boot'][:,1,None] ) )
        pybl = boot_object( fit['boot'][:,0,None] * np.cosh( (pxl[None,:]-offset) * fit['boot'][:,1,None] ) )
        pybu = boot_object( fit['boot'][:,0,None] * np.cosh( (pxu[None,:]-offset) * fit['boot'][:,1,None] ) )
        plot_fill( ax, px,  pyb['mean'] , pyb['std'] , c, lw            , **args )
        plot_fill( ax, pxl, pybl['mean'], pybl['std'], c, lw, alpha=0.33         )
        plot_fill( ax, pxu, pybu['mean'], pybu['std'], c, lw, alpha=0.33         )
      elif fit['type'] == 'function':
        function = eval(fit['function'])
        outer_list = outerlist(np.swapaxes(fit['boot'],0,1))
        pxm = np.linspace(lower,upper)
        pxl = np.linspace(xlims[0],lower)
        pxu = np.linspace(upper,xlims[1])
        pym = boot_object( np.swapaxes(np.asarray([ function(x,*outer_list) for x in pxm ]),0,1) )
        pyl = boot_object( np.swapaxes(np.asarray([ function(x,*outer_list) for x in pxl ]),0,1) )
        pyu = boot_object( np.swapaxes(np.asarray([ function(x,*outer_list) for x in pxu ]),0,1) )

        pyml      = np.percentile( pym['boot'], 10    , axis=0 )
        pymu      = np.percentile( pym['boot'], 100-10, axis=0 )
        pyll      = np.percentile( pyl['boot'], 10    , axis=0 )
        pylu      = np.percentile( pyl['boot'], 100-10, axis=0 )
        pyul      = np.percentile( pyu['boot'], 10    , axis=0 )
        pyuu      = np.percentile( pyu['boot'], 100-10, axis=0 )
        plot_fill( ax, pxm, pym['mean'], pym['std'], c, lw, confidence_upper=pymu, confidence_lower=pyml, **args )
        plot_fill( ax, pxl, pyl['mean'], pyl['std'], c, lw, confidence_upper=pylu, confidence_lower=pyll, alpha=0.33 )
        plot_fill( ax, pxu, pyu['mean'], pyu['std'], c, lw, confidence_upper=pyuu, confidence_lower=pyul, alpha=0.33 )
      else:
        print 'unknown fit type: ' + fit['type']
        exit(1)

  if legend:
    leg = ax.legend( fontsize=24, loc=kwargs.get('loc','best'), scatterpoints=1, numpoints=1, handlelength=1 )
    #leg = pl.legend( fontsize=24, loc=kwargs.get('loc','best') )
    #leg.get_frame().set_alpha(0.0)
  if 'xstride' in kwargs:
    pl.xticks( range(int(xlims[0]),int(xlims[1])+2,kwargs['xstride']), fontsize=fs )

  if 'dark' in kwargs:
    set_foregroundcolor(ax,'w')
    ax.grid(c='w')
  else:
    ax.grid()
  #pl.tight_layout()
  #print xlims
  #print ylims
  pl.xlim( xlims )
  pl.ylim( ylims )
  if 'save' in kwargs:
    pl.show()
    pl.savefig( kwargs['save'], transparent=kwargs.get('transparent',True), facecolor=bg, bbox_inches='tight' )
    pl.close()
  return pl, ax

def MultiPlot( data, **kwargs ):
  save = kwargs.pop('save','test.png')
  lowertime, uppertime = literal_eval( str(kwargs.get('time','[7,19]')) )
  lowertau , uppertau  = literal_eval( str(kwargs.get('tau' ,'[0,13]')) )
  #print lowertime, uppertime, lowertau, uppertau

  figure1, axes1 = pl.subplots( uppertime-lowertime+1, sharex=True, sharey=True, figsize=kwargs.pop('figsize',(14,18)) )
  for i, it in enumerate(range(lowertime,uppertime+1)):
    itplot_data = [dict(d) for d in data]
    #for data_index in range(len(itplot_data)):
    #  itplot_d = itplot_data[data_index]
    #  itplot_d['x'] = itplot_d['x2'][:it+1]
    #  itplot_d['y'] = itplot_d['y'][it,:it+1]
    #  itplot_d['s'] = itplot_d['s'][it,:it+1]
    #  if 'f' in itplot_d:
    #    points   = np.asarray([ (it,j) for j in range(it+1) ])
    #    x_points = map( Snd, PointsInside( itplot_d['f']['trapezoid'], points ) )
    #    if len(x_points):
    #      itplot_d['f']['range'] = (min(x_points),max(x_points))
    #    else:
    #      itplot_d['f']['range'] = (0,0)
    #if i == len(axes1)-1:
    #  kwargs['save']   = save
    #  kwargs['xlabel'] = 'tau'
    #figure1.subplots_adjust(hspace=0)
    #Plot( itplot_data, axes1[i], **kwargs )
    for data_index in range(len(itplot_data)):
      itplot_d = itplot_data[data_index]
      itplot_d['x'] = itplot_d['x2'][:]
      itplot_d['y'] = itplot_d['y'][it,:]
      itplot_d['s'] = itplot_d['s'][it,:]
      if 'f' in itplot_d:
        points   = np.asarray([ (it,j) for j in range(len(itplot_d['x'])) ])
        x_points = map( Snd, PointsInside( itplot_d['f']['trapezoid'], points ) )
        if len(x_points):
          itplot_d['f']['range'] = (min(x_points),max(x_points))
        else:
          itplot_d['f']['range'] = (0,0)
    if i == len(axes1)-1:
      kwargs['save']   = save
      kwargs['xlabel'] = 'tau'
    kwargs['red_after'] = it
    kwargs['xlims'] = [lowertau-1,max(uppertime-1,uppertau+1)]
    figure1.subplots_adjust(hspace=0)
    Plot( itplot_data, axes1[i], **kwargs )
  
  kwargs.pop('xlabel','')
  kwargs.pop('red_after','')
  figure2, axes2 = pl.subplots( uppertau-lowertau+1, sharex=True, sharey=True, figsize=kwargs.pop('figsize',(14,18)) )
  save2 = afterslash( kwargs.pop('save'), 'backwards_' )
  for i, itau in enumerate(range(lowertau,uppertau+1)):
    itplot_data = [dict(d) for d in data]
    for data_index in range(len(itplot_data)):
      itplot_d = itplot_data[data_index]
      itplot_d['x'] = itplot_d['x']
      itplot_d['y'] = itplot_d['y'][:,itau]
      itplot_d['s'] = itplot_d['s'][:,itau]
      if 'f' in itplot_d:
        points   = np.asarray([ (j,itau) for j in range(itplot_d['x'].shape[0]) ])
        x_points = map( Fst, PointsInside( itplot_d['f']['trapezoid'], points ) )
        if len(x_points):
          itplot_d['f']['range'] = (min(x_points),max(x_points))
        else:
          itplot_d['f']['range'] = (0,0)
    if i == len(axes2)-1:
      kwargs['save'] = save2
      kwargs['xlabel'] = 'time'
      kwargs['xlims'] = [lowertime-1,uppertime+1]
    figure2.subplots_adjust(hspace=0)
    Plot( itplot_data, axes2[i], **kwargs )

plot_extract_keys = {
  'time' : 'x'
, 'tau'  : 'x2'
, 'mean' : 'y'
, 'std'  : 's'
, 'fit' : 'f' 
, 'boot' : 'b'
, 'x' : 'x'
, 'y' : 'y'
, 's' : 's'
, 'p' : 'p'
, 'f' : 'f'
, 'l' : 'l'
, 'b' : 'b'
, 'o' : 'o'
, 'fmt' : 'fmt'
, 'cl' : 'cl'
, 'c'  : 'c'
}

def ChromaPlot( data, plotter=Plot, **kwargs ):
  plot_data = [{ v:d[k] for k,v in plot_extract_keys.items() if k in d } for d in data ]
  plotter( plot_data, **kwargs )

def _ChromaPlot( data, plotter=Plot, **kwargs ):
  if plotter == Plot:
    plotstr = 'plot'
  else:
    plotstr = 'mplot'
  yml_file = change_extension( kwargs['save'], 'yml' )
  open( yml_file, 'w+' ).write( json.dumps( { 'data':data, 'extra':dict(kwargs), 'plotter':plotstr }, cls=NumpyJsonEncoder ) )
