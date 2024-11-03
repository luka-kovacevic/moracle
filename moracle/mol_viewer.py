import zipfile
import os
import collections
import logging

viewer_description = """
PDB viewer using 3Dmol.js  
If using please cite:  
3Dmol.js: Molecular visualization with WebGL, Nicholas Rego, David Koes , Bioinformatics, Volume 31, Issue 8, April 2015, Pages 1322–1324, https://doi.org/10.1093/bioinformatics/btu829

See also:
https://huggingface.co/blog/spaces_3dmoljs
https://huggingface.co/spaces/simonduerr/3dmol.js/blob/main/app.py
"""

def parse_ligand_filename(filename: str) -> dict:
    """
    Parses an sdf filename to extract information.
    """
    if not filename.endswith(".sdf"):
        return {}

    basename = os.path.basename(filename).replace(".sdf", "")
    tokens = basename.split("_")
    rank = tokens[0]
    rank = int(rank.replace("rank", ""))
    if len(tokens) == 1:
        return {"filename": basename, "rank": rank, "confidence": None}

    con_str = tokens[1]
    conf_val = float(con_str.replace("confidence", ""))

    return {"filename": basename, "rank": rank, "confidence": conf_val}


def process_zip_file(zip_path: str):
    pdb_file = []
    sdf_files = []
    with zipfile.ZipFile(open(zip_path, "rb")) as my_zip_file:
        for filename in my_zip_file.namelist():
            # print(f"Processing file {filename}")
            if filename.endswith("/"):
                continue

            if filename.endswith(".pdb"):
                content = my_zip_file.read(filename).decode("utf-8")
                pdb_file.append({"path": filename, "content": content})

            if filename.endswith(".sdf"):
                info = parse_ligand_filename(filename)
                info["content"] = my_zip_file.read(filename).decode("utf-8")
                info["path"] = filename
                sdf_files.append(info)

    sdf_files = sorted(sdf_files, key=lambda x: x.get("rank", 1_000))

    return pdb_file, sdf_files


def gen_3dmol_vis(pdb_text: str, sdf_text: str):
    x = (
            """<!DOCTYPE html>
            <html>
            <head>    
        <meta http-equiv="content-type" content="text/html; charset=UTF-8" />
        <style>
        body{
            font-family:sans-serif
        }
        .mol-container {
        width: 100%;
        height: 600px;
        position: relative;
        }
        .mol-container select{
            background-image:None;
        }
        </style>
         <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.3/jquery.min.js" integrity="sha512-STof4xm1wgkfm7heWqFJVn58Hm3EtS31XFaagaa8VMReCXAkQnJZ+jEy8PCC/iT18dFy95WcExNHFTqLyp72eQ==" crossorigin="anonymous" referrerpolicy="no-referrer"></script>
        <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
        </head>
        <body>  
    
        <div id="container" class="mol-container"></div>
    
                <script>
                   let pdb = `""" + pdb_text + """`  
            
                   let sdf = `""" + sdf_text + """`

             $(document).ready(function () {
                let element = $("#container");
                let config = { backgroundColor: "white" };
                let viewer = $3Dmol.createViewer(element, config);
                
                viewer.addModel(pdb, "pdb", { format: "pdb" });
                
                pdb_model = viewer.getModel(0);
                // Cartoon view for protein
                // First argument is a selector
                pdb_model.setStyle({}, {cartoon: {color: "red", opacity: 0.8}});

                
                viewer.addModel(sdf, "sdf", {format: "sdf"});
                sdf_model = viewer.getModel(1);
                // Stick view for SDF
                sdf_model.setStyle({}, {stick: {color: "blue", opacity: 0.8}});

                
                
                viewer.zoomTo();
                viewer.render();
                // viewer.zoom(0.8, 2000);
              })
        </script>
        </body></html>"""
    )

    return f"""<iframe style="width: 100%; height: 600px" name="result" sandbox="allow-modals allow-forms 
    allow-scripts allow-same-origin allow-popups 
    allow-top-navigation-by-user-activation allow-downloads" allowfullscreen="" 
    allowpaymentrequest="" frameborder="0" srcdoc='{x}'></iframe>"""


def run_wrapper(output_file, *args) -> str:

    pdb_files, sdf_files = process_zip_file(output_file)
    # print(f"PDB file: {pdb_files}")
    pdb_file = pdb_files[0] if pdb_files else None
    for sdf_file in sdf_files:
        confidence = sdf_file.get("confidence", None)
        # rank1 has no confidence
        if confidence is None:
            continue
        label = f"Rank {sdf_file['rank']}. Confidence {confidence:.2f}"
        pdb_text = pdb_file['content'] if pdb_file else None
        sdf_text = sdf_file['content']
        output_viz = "Output visualisation unavailable"
        if pdb_text:
            logging.debug(f"Creating 3D visualisation")
            output_viz = gen_3dmol_vis(pdb_text, sdf_text)

    return output_viz