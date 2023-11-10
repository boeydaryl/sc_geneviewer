import dash
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
from dash import Input, Output, dcc, html
import scanpy as sc
import numpy as np
import pandas as pd

adata_path = "./genespector/genespector/data/Buffy210208_processed_labelled_9Nov23.h5ad"

adata = sc.read_h5ad(adata_path)

app = dash.Dash(external_stylesheets=[dbc.themes.SLATE])

gene_list = adata.raw.var_names.to_list()
cat_list = adata.obs.columns.to_list()
total_list = gene_list + cat_list
df = sc.get.obs_df(adata, keys=total_list, use_raw=True)
UMAP = adata.obsm['X_umap']
x_coords = UMAP[:, 0]
y_coords = UMAP[:, 1]

controls = dbc.Card(
    [
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
        html.Div(
            [
                dbc.Label('Gene Value Scale'),
                dcc.RangeSlider(
                    id='range_slider',
                    min=0,
                    max=10,
                    step=1,
                    value=[0,5])
            ]
        )
    ],
    body=True,
)

app.layout = dbc.Container(
    [
        html.H1("Wu, Boey, et al, 2022", style={'textAlign': 'center'}),
        html.Div(children='Data Visualizer for doi: 10.1182/bloodadvances.2022006969', style={
        'textAlign': 'center',
        }),
        dbc.Row(dbc.Col(dcc.Graph(id="cluster-graph"), md=5), align='center', className="h-85", justify='center'),
        dbc.Row(dbc.Col(controls, md=2), align='start', className='h-15', justify='center'),
        html.Div(children='Developed by Daryl Boey, Team Dahlin, Karolinska Institutet', style={
        'textAlign': 'center',
        })
        ],
        fluid=True,
        style={"height":"100vh"}
)

@app.callback(
    Output("cluster-graph", "figure"),
    Input("gene_var", "value"),
    Input("cat_data", "value"),
    Input("range_slider", "value")
)

def update_UMAP(gene_var, cat_data, range_slider):
    if gene_var in gene_list:
        vals = df[gene_var]
        fig = px.scatter(x=x_coords,
                         y=y_coords,
                         color=vals,
                         color_continuous_scale=['lightgrey', 'tomato', 'firebrick', 'darkred'],
                         range_color=range_slider,
                         title=gene_var,
                         opacity=0.5)
    if cat_data in cat_list:
        vals = df[cat_data]
        fig = px.scatter(x=x_coords,
                         y=y_coords,
                         color=vals,
                         title=cat_data,)
    fig.update_layout(
        xaxis_title="UMAP 1",
        yaxis_title="UMAP2",
        font=dict(
            family="Courier New, monospace",
            size=16,),
    )

    fig.update_traces(marker=dict(size=5))

    return fig

if __name__ == "__main__":
    app.run_server(debug=True)