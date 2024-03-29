import dash
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
from dash import Input, Output, dcc, html
import scanpy as sc
import numpy as np
import pandas as pd
import sys

adata_path = sys.argv[1]

adata = sc.read_h5ad(adata_path)

app = dash.Dash(external_stylesheets=[dbc.themes.SLATE])

gene_list = adata.raw.var_names.to_list()
cat_list = adata.obs.columns.to_list()
obsm_list = [X for X in adata.obsm]
total_list = gene_list + cat_list
df = sc.get.obs_df(adata, keys=total_list, use_raw=True)

controls = dbc.Card(
    [
        html.Div(
            [
                dbc.Label("Categorical Data"),
                dcc.Dropdown(
                    id="cat_data",
                    options=[
                        {"label": col, "value": col} for col in cat_list
                    ],
                    value="cell_type",
                ),
            ]
        ),
        html.Br(),
        html.Div(
            [
                dbc.Label("Gene"),
                dcc.Dropdown(
                    id="gene_var",
                    options=[
                        {"label": col, "value": col} for col in gene_list
                    ],
                    value="TPSB2",
                ),
            ]
        ),
        html.Br(),
        html.Div(
            [
                dbc.Label('Gene Value Scaling:'),
                dcc.RadioItems(
                    options=[{"label":"Default","value":"Default"},
                             {"label":"Manual","value":"Manual"}],
                    id='value_scale',
                    value="Default",
)
            ]
        ),
        html.Br(),
        html.Div(
            [
                dbc.Label('Gene Value Manual'),
                dcc.RangeSlider(
                    id='value_input',
                    min=0,
                    max=10,
                    step=1,
                    value=[0, 5],
)
            ]
        )
    ],
    body=True,
)

app.layout = dbc.Container(
    [
        html.H1("scGeneViewer", style={'textAlign': 'center'}),
        html.Div(children=html.P([f"Data Visualizer for Wu C, Boey D, et al."]),
                 style={
        'textAlign': 'center',
        }),
        dbc.Row([dbc.Col(dcc.Graph(id="cluster-graph"), md=4), dbc.Col(controls, md=2)],
                align='center',
                className="h-95",
                justify='center',
                ),
        html.Div(children='Usage: Select either [Categorical Data] or [Gene], and clear the other field.', style={
                'textAlign': 'center',
                }),
        html.Br(),
        html.Div([html.Button("Download annotated h5ad", id="btn_download_data"),
                  dcc.Download(id="download_data")],
                 style={'textAlign': 'center'}
                 ),
        html.Div(children='All figures and data are provided through a CC-BY 4.0 license, please cite the original article as referenced:', style={
                'textAlign': 'center',
                }),
        html.Div(children=html.P([f"Wu C.*, Boey D.*, Bril O., Grootens J., Vijayabaskar M.S., Sorini C., Ekoff M., Wilson N.K., Ungerstedt J.S., Nilsson G. & J.S. Dahlin.",
                                 html.Br(),f"Single-cell transcriptomics reveals the identity and regulators of human mast cell progenitors.",
                                 html.Br(), f"Blood Advances, 6: 4439-4449 (2022)."]),
                 style={
            'textAlign': 'center',
            'font-size' : '6'
        }),
        html.Div(children=['Developed by Daryl Boey, Team Dahlin, Karolinska Institutet, using Dash framework 2.1. Current version of dashboard available at ',
                           dcc.Link('Github', href='https://github.com/boeydaryl/sc_geneviewer')], style={
        'textAlign': 'center',
        })
        ],
        fluid=True,
        style={"height":"100vh"}
)

@app.callback(
    Output("cluster-graph", "figure"),
    Input("cat_data", "value"),
    Input("gene_var", "value"),
    Input("value_scale", "value"),
    Input("value_input", "value")
)

def update_graph(cat_data, gene_var, value_scale, value_input):
    rep = adata.obsm['X_umap']
    x_coords = rep[:, 0]
    y_coords = rep[:, 1]

    if gene_var in gene_list:
        vals = df[gene_var]
        if value_scale == "Default":
            val_input = [0, np.max(vals)]
        elif value_scale == "Manual":
            val_input = value_input
        fig = px.scatter(x=x_coords,
                         y=y_coords,
                         color=vals,
                         color_continuous_scale=['lightgrey', 'tomato', 'firebrick', 'darkred'],
                         range_color=val_input,
                         title=gene_var,
                         opacity=0.5)

    if cat_data in cat_list:
        if df[cat_data].dtype == 'category':
            color_dict = sc.plotting._tools.scatterplots._get_palette(adata, cat_data)
            vals = df[cat_data]
            fig = px.scatter(x=x_coords,
                             y=y_coords,
                             color=vals,
                             title=cat_data,
                             color_discrete_map=color_dict
                             )
        else:
            vals = df[cat_data]
            fig = px.scatter(x=x_coords,
                             y=y_coords,
                             color=vals,
                             title=cat_data,
                             )
    fig.update_layout(
        xaxis_title=str("UMAP" + " 1"),
        yaxis_title=str("UMAP" + " 2"),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        title_x=0.5,
        font=dict(
            family="Arial, monospace",
            size=14,),
        plot_bgcolor='rgba(0,0,0,0)',
        autosize=True
    )

    fig.update_traces(marker=dict(size=5))

    return fig

@app.callback(
    Output("download_data", "data"),
    Input("btn_download_data", "n_clicks"),
    prevent_initial_call=True,
)

def func(n_clicks):
    if n_clicks is not None:
        return dcc.send_file(
            "./data/Buffy210208_processed_labelled_9Nov23.h5ad"
        )

if __name__ == "__main__":
    app.run_server(debug=True)
