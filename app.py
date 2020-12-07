#!/usr/bin/env python
# coding: utf-8

# # Dashboard

# ## Import modules

# In[4]:


import six.moves.urllib.request as urlreq

import dash
import dash_bio as dashbio
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash_bio_utils import pdb_parser as parser
from dash.dependencies import Input, Output
import json
import tempfile

import copy
import pandas as pd

try:
    import mutagenesis_visualization as mut
except ModuleNotFoundError:  # This step is only for when I run the notebooks locally
    import sys
    sys.path.append('../../../')
    import mutagenesis_visualization as mut


# In[5]:


hras_RBD = mut.hras_RBD()


# # Main dashboard function

# In[30]:


def dashboard_3d_protein(self, pdb, chain='A', position_correction=0, **kwargs):
    '''
    Docstring
    '''

    # update kwargs
    temp_kwargs = copy.deepcopy(mut.code_kwargs.kwargs())
    temp_kwargs.update(kwargs)

    # Load data from pdb file
    model_data = _parse_pdb(pdb)

    # Run dashhoard
    _run_dashboard_3d_protein(
        self, model_data, temp_kwargs, position_correction, pdb, chain
    )


def _run_dashboard_3d_protein(
    self, model_data, temp_kwargs, position_correction, pdb, chain
):
    '''
    Create dashboard for main function.
    '''

    # Open app
    app = dash.Dash('__main__', external_stylesheets=[dbc.themes.BOOTSTRAP])
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
            model_data, self.dataframe.copy(), temp_kwargs['gof'],
            temp_kwargs['lof'], mode, position_correction, chain
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
    app.run_server(port=8082, )


# # Auxiliary functions

# In[7]:


def _parse_pdb(pdb):
    '''return the pdb in jason format'''
    # Parse pdb file
    modata = parser.create_data(pdb)

    # Put in jason format
    fmodel = files_data_style(modata)
    with open(fmodel) as fm:
        model_data = json.load(fm)

    return model_data


def _parse_styles_data(
    model_data, df, gof, lof, mode, position_correction, chain
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
        hras_RBD.dataframe.copy(), 0.1, -0.3, mode, position_correction=0
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


def files_data_style(content):
    '''
    Function to create the modelData and style files for molecule visualization
    '''
    fdat = tempfile.NamedTemporaryFile(suffix=".js", delete=False, mode='w+')
    fdat.write(content)
    dataFile = fdat.name
    fdat.close()
    return dataFile


# # Run Dashboard

# In[ ]:


pdb = 'data/5p21.pdb'
dashboard_3d_protein(hras_RBD, pdb)

