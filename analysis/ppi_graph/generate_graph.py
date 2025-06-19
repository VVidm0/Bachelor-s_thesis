import sys
import os
import webbrowser
import networkx as nx
import plotly.graph_objects as go
import plotly.io as pio
import json
import time
import dash
from dash import Dash, dcc, html, Input, Output
from threading import Timer

def extract_protein_info(json_file, fasta_file):
    with open(json_file, 'r') as i:
        data = json.load(i)

    names = {}
    lengths = {}
    tax = {}

    for record in data["results"]:
        uniprot_id = record["to"]["primaryAccession"]
        try:
            full_name = record["to"]["proteinDescription"]["recommendedName"]["fullName"]["value"]
        except KeyError:
            full_name = "N/A"
        try:
            seq_len = record["to"]["sequence"]["length"]
        except KeyError:
            seq_len = None

        names[uniprot_id] = full_name
        lengths[uniprot_id] = seq_len
        
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                line = line.strip()
                parts = line.split('|')
                if len(parts) >= 2:
                    uniprot_id = parts[1]
                else:
                    continue
                if "OS=" in line:
                    organism = line.split("OS=")[1].split("OX=")[0].strip()
                else:
                    organism = "Nieznany"
                tax[uniprot_id] = organism

    return names, lengths, tax
    
class PPIList:
  def __init__(self):
    self.list = {}

  def add_index(self, index):
    if index not in self.list:
      self.list[index] = []

  def add_info(self, index, info):
      self.add_index(index)
      if info not in self.list[index]:
        self.list[index].append(info)

  def display(self):
    for index, info_list in self.list.items():
      print(f"{index}: {info_list}")

  def get_indexes(self):
    return list(self.list.keys())

  def get_info_list(self):
    return [info for info_list in self.list.values() for info in info_list]

  def get_info(self, index):
    return self.list.get(index, [])

  def del_index(self,index):
    if index in self.list:
      del self.list[index]

def ppi_diff_sum(id_file1,id_file2):
  ppi1 = PPIList()
  ppi2 = PPIList()
  out = PPIList()

  with open(id_file1, 'r') as f1, open(id_file2, 'r') as f2:
    for line1 in f1:
      id1, id2 = line1.strip().split("\t")
      ppi1.add_info(id1, id2)
    for line2 in f2:
      id1, id2 = line2.strip().split("\t")
      ppi2.add_info(id1, id2)

    for host_protein in set(ppi1.get_indexes()) & set(ppi2.get_indexes()):
      pathogens1 = set(ppi1.get_info(host_protein))
      pathogens2 = set(ppi2.get_info(host_protein))

      if pathogens1 != pathogens2:
        union_pathogens = pathogens1.union(pathogens2)
        for pathogen in union_pathogens:
          out.add_info(host_protein, pathogen)

  return out

def ppi_diff(id_file1,id_file2):
  ppi1 = PPIList()
  ppi2 = PPIList()
  out = PPIList()
  with open(id_file1, 'r') as f1, open(id_file2, 'r') as f2:
    for line1 in f1:
      id1, id2 = line1.strip().split("\t")
      ppi1.add_info(id1, id2)
    for line2 in f2:
      id1, id2 = line2.strip().split("\t")
      ppi2.add_info(id1, id2)

    for host_protein, pathogen_proteins1 in ppi1.list.items():
        if host_protein in ppi2.list:
            set1 = set(pathogen_proteins1)
            set2 = set(ppi2.get_info(host_protein))
            diff_pathogens = set1.symmetric_difference(set2)  # czyli (set1 - set2) ∪ (set2 - set1)
            for pathogen in diff_pathogens:
                out.add_info(host_protein, pathogen)

    return out

#-------------------------------------------------------------------------------------------------------------
def ppi_to_graph(ppi):
    G = nx.Graph()
    for host_protein in ppi.get_indexes():
      info_list = ppi.get_info(host_protein)
      if info_list:
        for pathogen_protein in info_list:
          G.add_edge(host_protein, pathogen_protein)

    # Usuwanie węzłów z zerową liczba połączeń
    G.remove_nodes_from([n for n in G.nodes if G.degree[n] == 0])

    return G

def edge_trace(edge_list, color, name):
    return go.Scatter(
        x=edge_list[0], y=edge_list[1],
        mode='lines',
        line=dict(width=1, color=color),
        hoverinfo='none',
        name=name
    )

def create_dash_app(G, ppi_source, ppi_list1, ppi_list2, json_file1, json_file2, fasta_file1, fasta_file2, protein_type):
    pos = nx.kamada_kawai_layout(G, scale=60)

    if protein_type == "h":
        host_proteins = ppi_source.get_indexes()
        pathogen_proteins = ppi_list1.get_info_list() + ppi_list2.get_info_list()
    elif protein_type == "p":
        host_proteins = ppi_list1.get_indexes() + ppi_list2.get_indexes()
        pathogen_proteins = ppi_source.get_indexes()

    names1, lengths1, tax1 = extract_protein_info(json_file1, fasta_file1)
    names2, lengths2, tax2 = extract_protein_info(json_file2, fasta_file2)

    #Operator ** to rozpakowanie słownika. Dodaje wszystkie pary klucz-wartość
    names = {**names1, **names2}
    lengths = {**lengths1, **lengths2}
    tax = {**tax1, **tax2}

    nodes_to_remove = []
    # Sprawdzanie węzłów gospodarzy, którzy mają połączenia z innymi gospodarzami
    for node in G.nodes():
        if node in host_proteins:
            if any(neighbor in host_proteins for neighbor in G.neighbors(node)):
                nodes_to_remove.append(node)
    # Usuwanie węzłów gospodarzy, które mają połączenia z innymi gospodarzami
    if nodes_to_remove:
        G.remove_nodes_from(nodes_to_remove)
        nodes_to_remove.clear()


    node_infos = []

    for node in G.nodes():
        neighbors = list(G.neighbors(node))
        total_connections = len(neighbors)
        pathogen_connections = sum(1 for neighbor in neighbors if neighbor in pathogen_proteins)
        host_connections = total_connections - pathogen_connections
        seq_len = lengths.get(node, 0)
        full_name = names.get(node, 'Brak nazwy')
        tax_info = tax.get(node, 'Nieznany')
        node_infos.append({
            'id': node,
            'x': pos[node][0],
            'y': pos[node][1],
            'total': total_connections,
            'pathogen': pathogen_connections,
            'host': host_connections,
            'len': seq_len,
            'organism': tax_info,
            'hover': (
                f"{node}<br>"
                f"Total interactions: {total_connections}<br>"
                f"Interactions with pathogen proteins: {pathogen_connections}<br>"
                f"Interactions with host proteins: {host_connections}<br>"
                f"Full protein name: {full_name}<br>"
                f"Sequence length: {seq_len}<br>"
                f"Organism: {tax_info}"
            ),
            'group': 'host' if node in host_proteins else 'pathogen'
        })

    degrees = [info['total'] for info in node_infos]
    lengths_list = [info['len'] for info in node_infos if isinstance(info['len'], (int, float))]

    app = Dash()

    app.layout = html.Div([
        html.Div([
            html.Label('Minimum node degree:'),
            dcc.Slider(
                id='degree-slider',
                min=min(degrees),
                max=max(degrees),
                step=1,
                value=min(degrees),
                marks={deg: str(deg) for deg in range(min(degrees), max(degrees)+1, 5)},
            )
        ], style={'padding': '20px'}),

        html.Div([
            html.Label('Minimum sequence length:'),
            dcc.Slider(
                id='length-slider',
                min=min(lengths_list),
                max=max(lengths_list),
                step=10,
                value=min(lengths_list),
                marks={l: str(l) for l in range(min(lengths_list), max(lengths_list)+1, 200)},
            )
        ], style={'padding': '20px'}),
        
        html.Div([
            html.Label('Filter by organism:'),
            dcc.Dropdown(
                id='organism-dropdown',
                options=[{'label': org, 'value': org} for org in sorted(set(info['organism'] for info in node_infos))],
                value=[],
                multi=True,
                placeholder="Choose organism"
            )
        ], style={'padding': '20px'}),

        dcc.Graph(id='ppi-graph')
    ])

    @app.callback(
        Output('ppi-graph', 'figure'),
        Input('degree-slider', 'value'),
        Input('length-slider', 'value'),
        Input('organism-dropdown', 'value')
    )

    def update_graph(degree_threshold, len_threshold, organism_values):
        x_host, y_host, text_host, hover_host = [], [], [], []
        x_path, y_path, text_path, hover_path = [], [], [], []

        select_nodes = []

        for info in node_infos:
            if (info['total'] >= degree_threshold 
            and info['len'] >= len_threshold
            and (not organism_values or info['organism'] in organism_values)):
                select_nodes.append(info['id'])

                if info['group'] == 'host':
                    x_host.append(info['x'])
                    y_host.append(info['y'])
                    text_host.append(info['id'])
                    hover_host.append(info['hover'])
                else:
                    x_path.append(info['x'])
                    y_path.append(info['y'])
                    text_path.append(info['id'])
                    hover_path.append(info['hover'])

        inter_x, inter_y, phi_x, phi_y, hpidb_x, hpidb_y = [], [], [], [], [], []
        for node1, node2 in G.edges():
            if node1 in select_nodes or node2 in select_nodes:
                if node1 in host_proteins and node2 in pathogen_proteins:
                    host, pathogen = node1, node2
                elif node2 in host_proteins and node1 in pathogen_proteins:
                    host, pathogen = node2, node1
                else:
                    continue
            else:
                continue

            from_phi = pathogen in ppi_list1.get_info(host)
            from_hpidb = pathogen in ppi_list2.get_info(host)

            x_pair = [pos[host][0], pos[pathogen][0], None]
            y_pair = [pos[host][1], pos[pathogen][1], None]

            if from_phi and from_hpidb:
                inter_x += x_pair
                inter_y += y_pair
            elif from_phi:
                phi_x += x_pair
                phi_y += y_pair
            elif from_hpidb:
                hpidb_x += x_pair
                hpidb_y += y_pair

        edge_traces = [
            edge_trace((inter_x, inter_y), '#7D2CF4', f'Common interaction of datasets ({len(inter_x) // 3})'),
            edge_trace((phi_x, phi_y), '#F55D13', f'PHISTO ({len(phi_x) // 3})'),
            edge_trace((hpidb_x, hpidb_y), '#755856', f'HPIDB ({len(hpidb_x) // 3})')
        ]

        fig = go.Figure(
            data=edge_traces
            +
             [
                go.Scatter(x=x_host, y=y_host, mode='markers', name='Host',
                           text=text_host,
                           hovertext=hover_host,
                           hoverinfo='text',
                           marker=dict(color='#40F596', size=10)),
                go.Scatter(x=x_path, y=y_path, mode='markers', name='Pathogen',
                           text=text_path,
                           hovertext=hover_path,
                           hoverinfo='text',
                           marker=dict(color='#F52F17', size=10))
            ],
            layout=go.Layout(
                title='PPIs network graph',
                showlegend=True,
                hovermode='closest',
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                width=1600,
                height=1200,
                paper_bgcolor='rgba(0,0,0,0)',
                plot_bgcolor='rgba(0,0,0,0)'
            )
        )
        return fig

    return app
    
def open_browser():
    if not os.environ.get("WERKZEUG_RUN_MAIN"):
        webbrowser.open_new('http://127.0.0.1:1222/')
        
phisto_ids_final_file = "phisto_human_bacteria_id_uniprot_final.txt"
hpidb3_ids_final_file = "hpidb3_human_bacteria_id_uniprot_final.txt"

phisto_json_final_file = "phisto_human_bacteria_id_uniprot_final.json"
hpidb3_json_final_file = "hpidb3_human_bacteria_id_uniprot_final.json"

phisto_fasta_final_file = "phisto_human_bacteria_id_uniprot_final.fasta"
hpidb3_fasta_final_file = "hpidb3_human_bacteria_id_uniprot_final.fasta"

merged_ids_file = "merged_clustered_ids_file.txt"

phi = PPIList()
with open(phisto_ids_final_file, 'r') as f:
  for line in f:
    id1, id2 = line.strip().split("\t")
    phi.add_info(id1, id2)

hpidb = PPIList()
with open(hpidb3_ids_final_file, 'r') as f:
  for line in f:
    id1, id2 = line.strip().split("\t")
    hpidb.add_info(id1, id2)
    
graph_type = sys.argv[1]
if graph_type == "sum":
    inter = ppi_diff_sum(phisto_ids_final_file, hpidb3_ids_final_file)
    app = create_dash_app(ppi_to_graph(inter), inter, phi, hpidb, phisto_json_final_file, hpidb3_json_final_file, phisto_fasta_final_file, hpidb3_fasta_final_file, "h")
elif graph_type == "diff":
    diff = ppi_diff(phisto_ids_final_file, hpidb3_ids_final_file)
    app = create_dash_app(ppi_to_graph(diff), diff, phi, hpidb, phisto_json_final_file, hpidb3_json_final_file, phisto_fasta_final_file, hpidb3_fasta_final_file, "h")
elif graph_type == "all":
    merged = PPIList()
    with open(merged_ids_file, 'r') as f:
        for line in f:
            id1, id2 = line.strip().split("\t")
            merged.add_info(id1, id2)
    app = create_dash_app(ppi_to_graph(merged), merged, phi, hpidb, phisto_json_final_file, hpidb3_json_final_file, phisto_fasta_final_file, hpidb3_fasta_final_file, "h")
        
if __name__ == "__main__":
    Timer(1, open_browser).start()
    app.run(debug=True, port=1222)
