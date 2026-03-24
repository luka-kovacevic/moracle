import zipfile
import os
import html as _html
import logging


def parse_ligand_filename(filename: str) -> dict:
    if not filename.endswith(".sdf"):
        return {}
    basename = os.path.basename(filename).replace(".sdf", "")
    tokens = basename.split("_")
    rank = int(tokens[0].replace("rank", ""))
    if len(tokens) == 1:
        return {"filename": basename, "rank": rank, "confidence": None}
    conf_val = float(tokens[1].replace("confidence", ""))
    return {"filename": basename, "rank": rank, "confidence": conf_val}


def process_zip_file(zip_path: str):
    pdb_file = []
    sdf_files = []
    with zipfile.ZipFile(open(zip_path, "rb")) as zf:
        for filename in zf.namelist():
            if filename.endswith("/"):
                continue
            if filename.endswith(".pdb"):
                pdb_file.append({"path": filename, "content": zf.read(filename).decode("utf-8")})
            if filename.endswith(".sdf"):
                info = parse_ligand_filename(filename)
                info["content"] = zf.read(filename).decode("utf-8")
                info["path"] = filename
                sdf_files.append(info)
    sdf_files = sorted(sdf_files, key=lambda x: x.get("rank", 1_000))
    return pdb_file, sdf_files


_UNAVAILABLE_HTML = """
<div style="display:flex; align-items:center; justify-content:center;
            height:560px; background:#f8fafc; border-radius:10px;
            font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;">
  <div style="text-align:center; color:#94a3b8;">
    <div style="font-size:36px; margin-bottom:12px;">🔬</div>
    <div style="font-size:14px; font-weight:600; color:#64748b; margin-bottom:4px;">
      No DiffDock data available
    </div>
    <div style="font-size:12px;">
      This molecule has no pre-computed docking result for the selected protein.
    </div>
  </div>
</div>
"""


def gen_3dmol_vis(pdb_text: str, sdf_text: str) -> str:
    # Embed PDB/SDF as JSON strings to avoid quote/backtick escaping issues
    import json
    pdb_json = json.dumps(pdb_text)
    sdf_json = json.dumps(sdf_text)

    inner_html = f"""<!DOCTYPE html>
<html>
<head>
  <meta http-equiv="content-type" content="text/html; charset=UTF-8"/>
  <style>
    html, body {{ margin:0; padding:0; background:#f8fafc; }}
    #container {{ width:100%; height:600px; position:relative; }}
  </style>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.6.3/jquery.min.js"
    integrity="sha512-STof4xm1wgkfm7heWqFJVn58Hm3EtS31XFaagaa8VMReCXAkQnJZ+jEy8PCC/iT18dFy95WcExNHFTqLyp72eQ=="
    crossorigin="anonymous" referrerpolicy="no-referrer"></script>
  <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
</head>
<body>
  <div id="container"></div>
  <script>
    var pdb = {pdb_json};
    var sdf = {sdf_json};
    $(document).ready(function() {{
      var viewer = $3Dmol.createViewer($("#container"), {{ backgroundColor: "#f8fafc" }});
      viewer.addModel(pdb, "pdb");
      viewer.getModel(0).setStyle({{}}, {{ cartoon: {{ color: "spectrum", opacity: 0.85 }} }});
      viewer.addModel(sdf, "sdf");
      viewer.getModel(1).setStyle({{}}, {{ stick: {{ colorscheme: "cyanCarbon", radius: 0.15 }},
                                           sphere: {{ colorscheme: "cyanCarbon", radius: 0.3 }} }});
      viewer.zoomTo();
      viewer.render();
    }});
  </script>
</body>
</html>"""

    # HTML-encode for use as a srcdoc attribute value (double-quote delimited)
    encoded = _html.escape(inner_html, quote=True)
    return (
        f'<iframe style="width:100%;height:600px;border:none;" name="result" '
        f'sandbox="allow-modals allow-forms allow-scripts allow-same-origin '
        f'allow-popups allow-top-navigation-by-user-activation allow-downloads" '
        f'allowfullscreen="" frameborder="0" srcdoc="{encoded}"></iframe>'
    )


def run_wrapper(output_file, *args) -> str:
    if not output_file:
        return _UNAVAILABLE_HTML

    try:
        pdb_files, sdf_files = process_zip_file(output_file)
    except Exception:
        return _UNAVAILABLE_HTML

    pdb_file = pdb_files[0] if pdb_files else None

    # Use the highest-ranked SDF that has a confidence score
    for sdf_file in sdf_files:
        if sdf_file.get("confidence") is None:
            continue
        pdb_text = pdb_file["content"] if pdb_file else None
        sdf_text = sdf_file["content"]
        if pdb_text:
            logging.debug("Creating 3D visualisation")
            return gen_3dmol_vis(pdb_text, sdf_text)

    return _UNAVAILABLE_HTML
