'''
By M.F.L. Meersman.
In order to use this model, please first contact m.f.l.meersman@student.tudelft.nl to ask for my permission.
If you just want to have a look at the processed data, please check out file 'dataframe.csv'.
All the produced graphs can be viewed and downloaded as '.PNG' images. It is however recommended to download the '.html'
files since these allow you to interact with the figures such that you can find and look up specific molecules.

In order to process the data yourself, please uncomment '# main()' of line 336 in this file, then run this file.
Going to 'https://rcbs.org' leads you to the international protein databank. From here, you can copy the
respective 'pdb-code-names' and paste them into one of three 'xxxxx_code_list.csv' files to process those respective
molecules.

This model simulaties the time required for the disintegration of possible lifeforms on Enceladus.
Going to 'dataframe.csv'-file will lead you to an overview of all the simulated molecule structures.
You can open this file best in Excel or any spreadsheet viewer. The meaning of all columns is:
- Identifier: This indicates whether the molecule is DNA, fatty-acid or proten/amino-acid
- code: This is the universal naming of the molecule that is simulated. Typing in this code in Google can give you additional info about the molecule.
- Required Height [m]: This is the distance you need to be elevated above the surface in order to capture the molecule before it brakes down.
- id: binary indication for Identifier
- Area [Å2]: This is the geometric area of the molecule structure, i.e. the area that will capture radiation flux.
- Time required [s]: This is the time from exiting the vent untill the molecule breaks down.
- URL: the link that leads to all the information of the molecule.
- Molecular complexity [-]: expressed as a function of number of atoms and variations of atoms.
'''

import csv
import os
import urllib.request
import webbrowser
from os.path import isfile, join

import freesasa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.graph_objs as go
from Bio.PDB import PDBParser
from scipy.interpolate import make_interp_spline


def make_code_list(filename):
    file = open(filename)
    csvreader = csv.reader(file)
    rows = []
    for row in csvreader:
            rows.append(row)
    return np.array(rows).flatten()

def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    print("Downloading molecular data files...   ", outfnm)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        print("Couldn't find this molecule :(")
        return None

def url_list(pdbcode_list, viewurl= "https://www.rcsb.org/structure/"):
    lst = []
    for code in pdbcode_list:
        url = viewurl + code
        lst.append(url)
    return lst

def parser(name):
    parser = PDBParser()
    structure = parser.get_structure(name)
    result, sasa_classes = freesasa.calcBioPDB(structure)
    print('Calculating geometric area of molecule:', name)
    return result, sasa_classes

def smooth(x, y):
    X_Y_Spline = make_interp_spline(x, y)
    X_ = np.linspace(x.min(), x.max(), 20)
    Y_ = X_Y_Spline(X_)
    return(X_, Y_)

def do_click(trace, points, state):
    if points.point_inds:
        ind = points.point_inds[0]
        url = df.link.iloc[ind]
        webbrowser.open_new_tab(url)

def complexity(directory):
    lst = []
    for file in os.listdir(directory):
        b = os.path.getsize(directory + '/' +file)
        lst.append(b)
    return lst

def sasa(filename):
    structure = freesasa.Structure(filename)
    result = freesasa.calc(structure)
    area_classes = freesasa.classifyResults(result, structure)

    print("Total : %.2f A2" % result.totalArea())
    for key in area_classes:
        print (key, ": %.2f A2" % area_classes[key])
    return result.totalArea()

def integral(xarr, yarr):
    area = 0
    for i in range(5):
        avg = (yarr[i]+yarr[i+1])*0.5
        block = avg*(xarr[i+1]-xarr[i])
        area += block
    return area

def H(list):
    heights = []
    avg_cov = 6.239
    rad_per_sec_cm = 144675998.37885764
    for area in list:
        height = (avg_cov * 1000) / ((area/(10**16)) * rad_per_sec_cm * 10**6)
        heights.append(height)
    return heights



def main():
    # Enceladus
    pro = np.array([10**12, 10**5, 10**2, 10**1, 10**-1, 10**-3]) # [#/cm^2/day]
    el = np.array([10**12, 10**5, 10**2, 10**-1, 10**-3, 10**-6]) # [#/cm^2/day]
    x1 = np.array([0, 15, 75, 100, 300, 1000]) # [MeV]
    x2 = np.array([0, 10, 50, 100, 300, 1000]) # [MeV]

    daily_pro = integral(x1, pro)
    daily_el = integral(x2, el)
    daily_rad = daily_pro + daily_el
    rad_per_sec_cm = daily_rad/(24*60**2)

    plt.plot(x1, pro) #enceladus proton
    plt.plot(x2, el) #enceladus electron
    plt.title("Enceladus integral photon and electron")
    plt.ylabel("Integral Flux, #/cm2/day")
    plt.xlabel("Particle Energy [MeV]")
    plt.yscale("log")
    plt.grid()
    plt.figure()
    #plt.show()


    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x1, y=pro,
                             mode='lines+markers',
                             name='proton radiation', fill='tonexty'))
    fig.add_trace(go.Scatter(x=x2, y=el,
                             mode='lines+markers',
                             name='electron radiation', fill='tonexty'))
    fig.update_yaxes(type="log")
    fig.update_layout(title='Enceladus integral photon and electron by SATRAD model',
                      xaxis_title='Particle Energy [MeV]',
                      yaxis_title='Particle Energy [MeV]')
    #fig.show()





    # Molecule Areas
    avg_cov = 4.094 #[eV], average covalent bond strenghth for organic molecules with temperature correction to 200K


    ## DNA ##
    print("\nDNA")
    dna_code_lst = make_code_list('dna_code_list.csv')
    for code in dna_code_lst:
        download_pdb(code, 'datadirdna', downloadurl="https://files.rcsb.org/download/")
    DNA_area_lst = []
    for file in os.listdir('datadirdna'):
        f = os.path.join('datadirdna', file)
        # checking if it is a file
        if os.path.isfile(f):
            DNA_area = sasa(f)
            DNA_area_lst.append(DNA_area)
    DNA_area = np.average(DNA_area_lst)
    DNA_area_min = min(DNA_area_lst)



    ## Lipids ##
    print("\nLipids")
    lipid_code_lst = make_code_list('lipid_code_list.csv')
    for code in lipid_code_lst:
        download_pdb(code, 'datadirlipid', downloadurl="https://files.rcsb.org/download/")
    lipid_area_lst = []
    for file in os.listdir('datadirlipid'):
        f = os.path.join('datadirlipid', file)
        # checking if it is a file
        if os.path.isfile(f):
            lipid_area = sasa(f)
            lipid_area_lst.append(lipid_area)
    lipid_area = np.average(lipid_area_lst)
    lipid_area_min = min(lipid_area_lst)


    ## Amino-acids ##
    print("\nAmino-Acids")
    amino_code_lst = make_code_list('amino_code_list.csv')
    for code in amino_code_lst:
        download_pdb(code, 'datadiramino', downloadurl="https://files.rcsb.org/download/")
    amino_area_lst = []
    for file in os.listdir('datadiramino'):
        f = os.path.join('datadiramino', file)
        # checking if it is a file
        if os.path.isfile(f):
            amino_area = sasa(f)
            amino_area_lst.append(amino_area)
    amino_area = np.average(amino_area_lst)
    amino_area_min = min(amino_area_lst)


    ## Requirements on time and height
    print("Radiation energy (Electron and proton) per second per cm2 [MeV]:", rad_per_sec_cm)

    print('\nDNA')
    DNA_energy_received = (DNA_area/(10**16)) * rad_per_sec_cm * 10**6
    print('Energy received per second [eV]: ', DNA_energy_received)
    t_dna = avg_cov/DNA_energy_received
    print('Time req. for dissociation [s]: ', t_dna)
    print('Height requirement [m]', t_dna * 1000)
    DNA_height_min = (avg_cov * 1000) / ((DNA_area_min/(10**16)) * rad_per_sec_cm * 10**6)
    print('Height requirement var [m]', t_dna*1000 + DNA_height_min, t_dna*1000 - DNA_height_min)

    print('\nLIPIDS')
    lipid_energy_received = (lipid_area/(10**16))*rad_per_sec_cm * 10**6
    print('Energy received per second [eV]: ', lipid_energy_received)
    t_lipid = avg_cov/lipid_energy_received
    print('Time req. for dissociation [s]: ', t_lipid)
    print('Height requirement [m]', t_lipid*1000)
    lipid_height_min = (avg_cov * 1000) / ((lipid_area_min/(10**16)) * rad_per_sec_cm * 10**6)
    print('Height requirement var [m]', t_lipid*1000 + lipid_height_min, t_lipid*1000 - lipid_height_min)

    print('\nAMINO-ACIDS')
    amino_energy_received = (amino_area/(10**16))*rad_per_sec_cm * 10**6
    print('Energy received per second [eV]: ', amino_energy_received)
    t_amino = avg_cov/amino_energy_received
    print('Time req. for dissociation [s]: ', t_amino)
    print('Height requirement [m]', t_amino * 1000)
    amino_height_min = (avg_cov * 1000) / ((amino_area_min/(10**16)) * rad_per_sec_cm * 10**6)
    print('Height requirement var [m]', t_amino*1000 + amino_height_min, t_amino*1000 - amino_height_min)


    dna_names = ['DNA'] * len(DNA_area_lst)
    dna_id = [0]*len(dna_names)
    lipid_names = ['lipid'] * len(lipid_area_lst)
    lipid_id = [0.1]*len(lipid_names)
    amino_names = ['amino'] * len(amino_area_lst)
    amino_id = [0.2] *len(amino_names)
    dna_codes = [f[:-4] for f in os.listdir('datadirdna') if isfile(join('datadirdna', f))]
    dna_urls = url_list(dna_codes, viewurl="https://www.rcsb.org/structure/")
    lipid_codes = [f[:-4] for f in os.listdir('datadirlipid') if isfile(join('datadirlipid', f))]
    lipid_urls = url_list(lipid_codes, viewurl="https://www.rcsb.org/structure/")
    amino_codes = [f[:-4] for f in os.listdir('datadiramino') if isfile(join('datadiramino', f))]
    amino_urls = url_list(amino_codes, viewurl="https://www.rcsb.org/structure/")
    urls = dna_urls + lipid_urls + amino_urls
    codes = dna_codes + lipid_codes + amino_codes
    dna_comp = complexity('datadirdna')
    lipid_comp = complexity('datadirlipid')
    amino_comp = complexity('datadiramino')


    d = {'Identifier': dna_names + lipid_names + amino_names, 'code': codes,
         'Required Height before dissociation [m]': H(DNA_area_lst) + H(lipid_area_lst) + H(amino_area_lst),
         "id": dna_id + lipid_id + amino_id, "Area [Å2]": DNA_area_lst + lipid_area_lst + amino_area_lst,
         "Time to measure before dissociation [t]": [i / 1000 for i in H(DNA_area_lst) + H(lipid_area_lst) + H(amino_area_lst)],  "URL": urls,
         "Molecular Complexity [-]": dna_comp+lipid_comp+amino_comp}
    df = pd.DataFrame(data=d)


    df.to_csv('dataframe.csv')


def generateresults():
    df = pd.read_csv('dataframe.csv')

    fig = px.scatter(df, x="Molecular Complexity [-]", y="Required Height before dissociation [m]",
                     size="Area [Å2]", color="Identifier",
                     hover_name="code", size_max=35, hover_data=["URL", "Time to measure before dissociation [t]"],
                     title="Required height above Enceladus to capture molecules before dissociating in a 1 [km/s] flow")
    fig.update_yaxes(type="log")
    # fig.update_xaxes(type="log")
    fig.show()
    fig.write_html("Height_Complexity.html")

    fig = px.histogram(df, x="Time to measure before dissociation [t]", color="Identifier", marginal="rug",
                       hover_data=df.columns, title="Dissociation time on Enceladus in a 1 [km/s] flow")
    fig.show()
    fig.write_html("Time_Histogram.html")

    ''' This is an interactive application that allows people 
    to click different points in plots and gets them redirected 
    to all the information related to the molecule. For now in Alpha state

    app = dash.Dash(__name__)

    app.layout = html.Div(
        [
            html.A(
                children="Please click this link",
                id="link",
                href="https://dash.plot.ly",
                target="_blank",
            ),
            dcc.Graph(
                id="figure",
                # specify custom_data to be the "urls" column
                figure=px.scatter(df, x="molecule_complexity", y="Required Height before dissociation [m]",
                     size="Area [angstron2]",  color="Identifier",
                     hover_name="Identifier", size_max=35, custom_data=("urls",)),
            ),
        ]
    )

    @app.callback(
        Output('link', 'href'),
        [Input('figure', 'hoverData')])
    def display_hover_data(hoverData):
        if hoverData:
            target = hoverData['points'][0]['customdata']
            return target
        else:
            raise PreventUpdate

    app.run_server(debug=False)
    '''


# Run once in order to generate all the data
if __name__ == "__main__":
    # main() # Comment this line out if you ran main() the simulations.
    generateresults()
