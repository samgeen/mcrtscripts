'''
Shade area of figure after 4 Myr
Sam Geen, May 2017
'''

def run(ax):
    xl = list(ax.get_xlim())
    yl = list(ax.get_ylim())
    xl[0] = 4.0
    ax.fill([xl[0],xl[1],xl[1],xl[0]],
            [yl[0],yl[0],yl[1],yl[1]],
            color=r"#888888",hatch='\\',alpha=0.2)
