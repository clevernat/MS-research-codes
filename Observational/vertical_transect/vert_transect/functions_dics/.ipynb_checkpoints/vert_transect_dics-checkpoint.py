#### File Paths to all data ####
def main_dics(folder, nc_file):

    f = {'kazr': [],
             'cloud': [],
             'mpl': [],
             'sur': [],
             'is': [],
             'mwr': [],
             'mwrlos': [],
             'vwp': [],
             'vpt': [],
             'iwp': [],
             'cld': [],
             'aeri': [],
             'micro': [],
             }
    fp ={'kazr': "KAZR/"+folder+"/"+nc_file,
                 'cloud': "KAZR/anxarsclkazrbnd1kolliasM1/anxarsclkazrbnd1kolliasM1.c1.20",
                 'mpl': "MPL/anx30smplcmask1zwangM1.c1.20",
                 'sur': "surface_cor/anxsurM1.20",
                 'is': "interpsonde/anxinterpolatedsondeM1.c1.20",
                 'mwr':"MWRRET/anxmwrret1liljclouM1.c2.20",
                 'mwrlos': "MWR/anxmwrlosM1.b1.20",
                 'vwp': "WindProfiler/anx1290rwpwindconM1.a1.20",
                 'vpt': "KASACR/VPT/anxkasacrcfrvptM1.a1.20",
                 'iwp': "iwc/anxiwcmicroM1.20",
                 'cld': "cldphase/anxthermocldphaseM1.c1.20",
                 'aeri': "aeri/anxtropoeM1.c1.20",
                 'micro': "microbaseplus/anxmicrobasekaplusM1.c1.20",
                }
    return f, fp

def vars_dics():
    v = {'kazr': {'h': [],
                 'ref': [],
                 'vel': [],
                 'sw': [],
                 'noise': [],
                 'time': [],
                },
            'vpt': {'h': [],
                'ref': [],
                'vel': [],
                'sw': [],
                'noise': [],
                'time': [],
                },
            'cloud': {'fcth': [],
                  'cb': [],
                  'cd': [],
                  'time': [],
                 },
            'mwr': {'lwp': [],
                'pwv': [],
                'time': [],
               },
            'is': {'h': [],
               'pwv': [],
               'q': [],
               't': [],
               'dp': [],
               'pres': [],
               'wspd': [],
               'u': [],
               'v': [],
               'ept': [],
               'time': [],
              },
            'sur': {'t': [],
                'wspd': [],
                'pres': [],
                'prec': [],
                'dirc': [],
                'time': [],
               },
            'iwp': {'iwp': [],
               },
            'cld': {'h': [],
                'phase': [],
                'time': [],
               },
            'aeri': {'h': [],
                     't': [],
                     'RH': [],
                     'gamma': [],
                     'rmsr': [],
                     'time': [],
                     },
            'micro': {'h': [],
                      'lwc': [],
                      'lwc_u': [],
                      'iwc': [],
                      'iwc_u': [],
                      'time': [],
                        }
            }

    fc = {'kazr': [],
          'cloud': [],
          'mpl': [],
          'sur': [],
          'is': [],
          'mwr': [],
          'mwrlos': [],
          'vwp': [],
          'vpt': [],
          'iwp': [],
          'cld': [],
          'aeri': [],
          'micro': [],
          }
    return v, fc

def panel_name_dic():
    pn = {"reflectivity": -1,
          "velocity": -1,
          "spectrum width": -1,
          "cloud phase": -1,
          "water path": -1,
          "surface": -1,
          "AERI": -1,
          "MPL": -1,
          "iwc": -1,
          }
    return pn
