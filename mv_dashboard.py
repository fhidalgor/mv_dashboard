#!/usr/bin/env python
# coding: utf-8

# # Dashboard

# ## Import modules

# In[ ]:


import dash
import dash_bio as dashbio
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
#from dash_bio_utils import pdb_parser as parser
from dash.dependencies import Input, Output
import json
import tempfile
import re
import os
from shutil import copy2

import parmed as pmd


import copy
import pandas as pd

try:
    import mutagenesis_visualization as mut
    jupyterlab = False
    hras_RBD=mut.hras_RBD()
    pdb = 'data/5p21.pdb'
except ModuleNotFoundError:  # This step is only for when I run the notebooks locally
    import sys
    sys.path.append('../../../')
    import mutagenesis_visualization as mut
    __name__ = '__main__'
    jupyterlab = True # for local use in jupyter lab
    hras_RBD=mut.hras_RBD()
    pdb = '../../data/5p21.pdb'


# # Auxiliary functions

# In[ ]:


def _parse_pdb(pdb):
    '''return the pdb in jason format'''
    # Parse pdb file
    modata = parser(pdb)
    #modata = parser.create_data(pdb)

    # Put in jason format
    fmodel = _files_data_style(modata)
    with open(fmodel) as fm:
        model_data = json.load(fm)

    return model_data


def _parse_styles_data(
    self,
    model_data,
    df,
    gof,
    lof,
    mode,
    position_correction,
    chain,
):
    '''
    From a dataframe with enrichment scores, this function will return a jason file
    with the color of each atom.
    
    Returns
    -------
    styles_data : jason file
    '''

    # Create empty dict
    styles_data = {}

    # Calculate df with colors
    df_color = _add_color(
        self.dataframe.copy(), gof, lof, mode, position_correction=0
    )

    # Iterate over parsed pdb
    for item in model_data['atoms']:
        #if item['chain'] != 'A':  # only color atoms from selected chain
        #break
        try:
            style_atom = {
                'color': _assign_color_jason(df_color, item['residue_index']),
                'visualization_type': 'cartoon'
            }
            styles_data[str(item['serial'])] = style_atom
        except IndexError:  # in case we runt out of index
            pass
    return styles_data


def _assign_color_jason(df, residue):
    '''
    Give a color to the atom based on the enrichment score of that residue.
    As an input takes the dataframe that _add_color returns.
    '''
    return df.loc[df['Position_Corrected'] == residue, 'Color'].iloc[0]


def _add_color(df, gof, lof, mode, position_correction):
    '''You input the dataframe. Removes stop codons. 
    Returns the positions that are going to be colored blue,red and white'''

    # Correct position
    df['Position_Corrected'] = df['Position'] + position_correction

    # Add dummy color column
    red = '#FD3216'
    blue = '#6A76FC'
    green = '#16FF32'

    # Select grouping
    if mode == 'mean':
        df_grouped = df.groupby(['Position'], as_index=False).mean()
    else:
        df_grouped = df.loc[df['Aminoacid'] == mode]

    # Color of mutations
    df_grouped['Color'] = green
    df_grouped.loc[df_grouped['Score'] < lof, 'Color'] = blue
    df_grouped.loc[df_grouped['Score'] > gof, 'Color'] = red

    return df_grouped


def _files_data_style(content):
    '''
    Function to create the modelData and style files for molecule visualization
    '''
    fdat = tempfile.NamedTemporaryFile(suffix=".js", delete=False, mode='w+')
    fdat.write(content)
    dataFile = fdat.name
    fdat.close()
    return dataFile


# ## PDB parser

# In[ ]:


def parser(pdb_path):
    """
    Parse the protein data bank (PDB) file to generate
    input modelData
    @param pdb_path
    Name of the biomolecular structure file in PDB format
    """

    # Create local copy of temp file
    copy2(pdb_path, './tmp.pdb')

    # Use parmed to read the bond information from temp file
    top = pmd.load_file('tmp.pdb')

    # Remove the created temp file
    #os.remove('tmp.pdb') #was giving error when uploading to huroku

    # Read PDB file to create atom/bond information
    with open(pdb_path, 'r') as infile:
        # store only non-empty lines
        lines = [l.strip() for l in infile if l.strip()]

    # Initialize all variables
    var_nchains = []
    serial = []
    atm_name = []
    res_name = []
    chain = []
    res_id = []
    positions = []
    occupancy = []
    temp_factor = []
    atom_type = []
    ct = 0

    datb = {
        'atoms': [],
        'bonds': []
    }

    # Variables that store the character positions of different
    # parameters from the molecule PDB file
    serialpos = [6, 11]
    atm_namepos = [12, 16]
    r_namepos = [17, 20]
    chainpos = [21, 22]
    r_idpos = [22, 26]
    xpos = [30, 38]
    ypos = [38, 46]
    zpos = [46, 54]
    occupos = [54, 60]
    bfacpos = [60, 66]
    atm_typepos = [77, 79]

    for l in lines:
        line = l.split()
        if "ATOM" in line[0] or "HETATM" in line[0]:
            serial.append(int(l[serialpos[0]:serialpos[1]]))
            atm_name.append(l[atm_namepos[0]:atm_namepos[1]].strip())
            val_r_name = l[r_namepos[0]:r_namepos[1]].strip()
            res_name.append(val_r_name)
            chain_val = l[chainpos[0]:chainpos[1]].strip()
            chain.append(chain_val)
            if chain_val not in var_nchains:
                var_nchains.append(chain_val)
            val_r_id = int(l[r_idpos[0]:r_idpos[1]])
            res_id.append(val_r_id)
            x = float(l[xpos[0]:xpos[1]])
            y = float(l[ypos[0]:ypos[1]])
            z = float(l[zpos[0]:zpos[1]])
            positions.append([x, y, z])
            occupancy.append(l[occupos[0]:occupos[1]].strip())
            temp_factor.append(l[bfacpos[0]:bfacpos[1]].strip())
            atom_type.append(l[atm_typepos[0]:atm_typepos[1]].strip())
            ct += 1

    # Create list of atoms
    tmp_res = res_id[0]
    resct = 1
    for i in range(len(chain)):  # pylint: disable=consider-using-enumerate
        if tmp_res != res_id[i]:
            tmp_res = res_id[i]
            resct += 1
        datb['atoms'].append({
            "name": atm_name[i],
            "chain": chain[i],
            "positions": positions[i],
            "residue_index": resct,
            "element": atom_type[i],
            "residue_name": res_name[i] + str(res_id[i]),
            "serial": i,
        })

    # Create list of bonds using the parmed module
    for i in range(len(top.bonds)):
        bondpair = top.bonds[i].__dict__
        atom1 = re.findall(r"\[(\d+)\]", str(bondpair['atom1']))
        atom2 = re.findall(r"\[(\d+)\]", str(bondpair['atom2']))
        datb['bonds'].append({
            'atom2_index': int(atom1[0]),
            'atom1_index': int(atom2[0])
        })

    return json.dumps(datb)


# # Main dashboard function

# In[ ]:


self = hras_RBD
chain='A'
position_correction=0
# update kwargs
temp_kwargs = copy.deepcopy(mut.code_kwargs.kwargs())

# Load data from pdb file
model_data = _parse_pdb(pdb)

# Open app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
server = app.server

# App layout
app.layout = html.Div([
    dbc.Row(
        dbc.Col(
            html.H3("Site-saturation mutagenesis of H-Ras"),  # title
            width={'size': 12},
        ),
    ),
    dbc.Row(
        dbc.Col(
            dcc.Graph(  # heatmap
                id="heatmap",
                figure={},
                config={'displayModeBar': False},
            ),
            width=12,
        )
    ),
    dbc.Row([
        dbc.Col(
            html.Div("Filter the data by:"),
            width={'size': 4, 'offset': 3}
        ),
        dbc.Col(
            dcc.Dropdown(
                id='dropdown',
                placeholder='choose amino acid',
                clearable=False,
                value='mean',
                options=[
                    {'label': 'Mean', 'value': 'mean'},
                    {'label': 'Alanine', 'value': 'A'},
                    {'label': 'Arginine', 'value': 'R'},
                    {'label': 'Asparagine', 'value': 'N'},
                    {'label': 'Aspartic acid ', 'value': 'D'},
                    {'label': 'Cysteine', 'value': 'C'},
                    {'label': 'Glutamine', 'value': 'Q'},
                    {'label': 'Glutamic acid ', 'value': 'E'},
                    {'label': 'Glycine', 'value': 'G'},
                    {'label': 'Histidine', 'value': 'H'},
                    {'label': 'Isoleucine', 'value': 'I'},
                    {'label': 'Leucine', 'value': 'L'},
                    {'label': 'Lysine', 'value': 'K'},
                    {'label': 'Methionine', 'value': 'M'},
                    {'label': 'Phenylalanine', 'value': 'F'},
                    {'label': 'Proline', 'value': 'P'},
                    {'label': 'Serine', 'value': 'S'},
                    {'label': 'Threonine', 'value': 'T'},
                    {'label': 'Tryptophan', 'value': 'W'},
                    {'label': 'Tyrosine', 'value': 'Y'},
                    {'label': 'Valine', 'value': 'V'},
                    #{'label': 'Stop codon', 'value': '*'},
                ]
            ),
            width=4,
        ),
    ]),
    dbc.Row([
        dbc.Col(
            id='moleculeviewer',
            children={},
            width={"size": 5, 'order':1, "offset": -1},
        ),
        dbc.Col(
            [dbc.Row(
            dcc.Graph(
                id="mean",
                figure={},
                config={'displayModeBar': False},
            ),
        ),
             dbc.Row([dbc.Col(
            dcc.Graph(
                id="scatter_3d",
                figure={},
                config={'displayModeBar': False},
            ),width={"size": 7, 'order':1},
             ),dbc.Col(
            dcc.Graph(
                id="histogram",
                figure={},
                config={'displayModeBar': False},
            ),width={"size": 5, 'order':2},
             )],no_gutters=True),],width={"size": 7, 'order':2},)
    ],
            no_gutters=True),
])

@app.callback([Output('moleculeviewer', 'children')],
              [Input('dropdown', 'value')])
def update_molecule3d(mode):
    '''
    Call the dashbio.molecule3dviewer and updated the coloring based on user input.
    '''
    # Calculate styles based on enrichment scores for the 3d viewer
    styles_data = _parse_styles_data(
        self,
        model_data,
        self.dataframe.copy(),
        0.1,
        -0.3,
        mode,
        position_correction,
        chain,
    )

    return [html.Div(dashbio.Molecule3dViewer( # 3d molecule
                        modelData=model_data,
                        styles=styles_data,
                        selectionType='Chain',
                    ))]

@app.callback([
    Output('heatmap', 'figure'),
    Output('mean', 'figure'),
    Output('scatter_3d', 'figure'),
    Output('histogram', 'figure'),
], [
    Input('dropdown', 'value'),
])
def update_graphs(mode='mean'):
    '''
    Aux function to update the plotly figures based on the user input.
    '''

    heatmap = self.heatmap_plotly(
        return_plot_object=True,
        show=False,
        title='Heatmap',
        figsize=(8, 2.5),
    )
    mean = self.mean_plotly(
        mode=mode,
        return_plot_object=True,
        show=False,
        figsize=(6, 2.5),
    )

    scatter_3d = self.scatter_3D_plotly(
        mode=mode,
        pdb_path=pdb,
        return_plot_object=True,
        show=False,
        figsize=(3, 3),
        title='C-α carbons',
    )

    histogram = self.histogram_plotly(
        mode=mode,
        return_plot_object=True,
        show=False,
        figsize=(3, 3),
    )
    return heatmap, mean, scatter_3d, histogram

# Run server
if jupyterlab:
    app.run_server(port=8083)        
else:
    if __name__ == '__main__':
        app.run_server(debug=True)


# # Function to call dashboard

# In[ ]:





# # Run Dashboard

# In[ ]:





# In[ ]:




