from flask import Flask, request, jsonify, make_response
from flask_cors import CORS
from ProteinToMRNA import ProteinToMRNA
from ProteinToSiRNA import ProteinToSiRNA
from ProteinToCRISPR import ProteinToCRISPR
import requests
from Bio import Entrez
from datetime import datetime

app = Flask(__name__)
CORS(app)

# Initialize converters (assuming these are real implementations)
mrna_converter = ProteinToMRNA(email="jshaigler07@gmail.com")
sirna_converter = ProteinToSiRNA(email="jshaigler07@gmail.com")
crispr_converter = ProteinToCRISPR(email="jshaigler07@gmail.com")

@app.route('/get_mrna', methods=['POST'])
def get_mrna():
    """
    API endpoint to receive a protein name and return the optimized mRNA sequence.
    """
    data = request.get_json()
    protein_name = data.get('protein_name')
    if not protein_name:
        return jsonify({"error": "Protein name is required"}), 400

    result = mrna_converter.get_optimized_mrna(protein_name, use_structure_aware=True, use_pseudouridine=True)

    if result["success"] and "Homo sapiens" in result["description"]:
        return jsonify(result), 200
    else:
        return jsonify({"error": "Protein not found for Homo sapiens or other error occurred."}), 500

@app.route('/get_sirna', methods=['POST'])
def get_sirna():
    """
    API endpoint to receive a protein name and return siRNA candidates.
    """
    data = request.get_json()
    protein_name = data.get('protein_name')
    if not protein_name:
        return jsonify({"error": "Protein name is required"}), 400

    result = sirna_converter.get_sirna_from_protein(protein_name)

    if result["success"] and "Homo sapiens" in result["description"]:
        return jsonify(result), 200
    else:
        return jsonify({"error": "Protein not found for Homo sapiens or other error occurred."}), 500

@app.route('/design_crispr', methods=['POST'])
def design_crispr():
    """
    API endpoint to design a CRISPR complex for a given gene and edit function.
    """
    data = request.get_json()
    gene_name = data.get('gene_name')
    edit_function = data.get('edit_function')
    if not gene_name or not edit_function:
        return jsonify({"error": "Gene name and edit function are required"}), 400

    try:
        result = crispr_converter.get_crispr_design(gene_name, edit_function)
        if result["success"]:
            return jsonify(result), 200
        else:
            return jsonify({"error": result["error"]}), 500
    except Exception as e:
        print(f"Error in /design_crispr: {e}")
        return jsonify({"error": "An internal server error occurred. Please try again later."}), 500

@app.route('/literature_check', methods=['POST'])
def literature_check():
    """
    API endpoint to retrieve related publications and patents for a protein or gene.
    """
    data = request.get_json()
    query = data.get('protein_name') or data.get('gene_name')
    if not query:
        return jsonify({"error": "Protein or gene name is required"}), 400

    try:
        europe_pmc_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
        params = {
            "query": f"{query} AND (CRISPR OR siRNA OR mRNA)",
            "format": "json",
            "pageSize": 5,
        }
        response = requests.get(europe_pmc_url, params=params)
        results = response.json().get("resultList", {}).get("result", [])
        papers = []
        for item in results:
            papers.append({
                "title": item.get("title"),
                "journal": item.get("journalTitle", item.get("source", "")),
                "year": item.get("pubYear"),
                "type": item.get("pubType"),
                "link": f"https://europepmc.org/article/{item.get('source')}/{item.get('id')}"
            })
        return jsonify({"query": query, "results": papers}), 200

    except Exception as e:
        print(f"Error in /literature_check: {e}")
        return jsonify({"error": "An internal error occurred while searching Europe PMC."}), 500

# -------------------------------------------------------------------
# Advanced Features Inspired by Google Moonshot Labs
# -------------------------------------------------------------------

@app.route('/visualize_structure', methods=['POST'])
def visualize_structure():
    """
    Enhanced 3D Molecular Visualization.
    Given a molecule ID (e.g., PDB ID), fetch metadata from the RCSB PDB API and return a viewer URL.
    """
    data = request.get_json()
    molecule_id = data.get('molecule_id')
    if not molecule_id:
        return jsonify({"error": "Molecule ID is required"}), 400

    try:
        pdb_info = fetch_pdb_metadata(molecule_id)
        viewer_url = f"https://nglviewer.org/app/?pdb={molecule_id}"
        return jsonify({
            "molecule_id": molecule_id,
            "structure_info": pdb_info,
            "viewer_url": viewer_url,
            "timestamp": datetime.now().isoformat()
        }), 200

    except Exception as e:
        print(f"Error in structure visualization: {e}")
        return jsonify({"error": str(e)}), 500

def fetch_pdb_metadata(molecule_id: str):
    """
    Fetch metadata for a PDB structure using the RCSB PDB API.
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{molecule_id}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        resolution = data.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
        method = data.get("exptl", [{}])[0].get("method", "Unknown")
        return {"id": molecule_id, "resolution": resolution, "method": method}
    else:
        return {"id": molecule_id, "error": "PDB entry not found."}

def build_preflight_response():
    response = make_response()
    response.headers.add('Access-Control-Allow-Origin', '*')
    response.headers.add('Access-Control-Allow-Headers', 'Content-Type')
    response.headers.add('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
    return response

@app.route('/simulate_gene_editing', methods=['POST', 'OPTIONS'])
def simulate_gene_editing():
    if request.method == 'OPTIONS':
        return build_preflight_response()
    
    data = request.get_json()
    gene = data.get('gene')
    design_type = data.get('design_type')
    if not gene or not design_type:
        return jsonify({'error': 'Gene and design type are required.'}), 400
    
    # Dummy simulation result (adjust with your actual logic)
    response_data = {
        'off_target_risk': 'Low',
        'predicted_outcome': 'Successful editing',
        'simulated_efficiency': '90%'
    }
    return jsonify(response_data), 200

@app.route('/delivery_strategy', methods=['POST', 'OPTIONS'])
def delivery_strategy():
    if request.method == 'OPTIONS':
        return build_preflight_response()
    
    data = request.get_json()
    tissue_type = data.get('tissue_type')
    if not tissue_type:
        return jsonify({'error': 'Tissue type is required.'}), 400
    
    # Dummy delivery strategy (adjust with your actual logic)
    response_data = {
        'recommended_strategies': ['Lipid nanoparticles', 'AAV vectors', 'Electroporation']
    }
    return jsonify(response_data), 200

if __name__ == '__main__':
    app.run(debug=True)
