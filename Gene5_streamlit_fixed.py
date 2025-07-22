
[14:31:09] 🖥 Provisioning machine...
[14:31:09] 🎛 Preparing system...
[14:31:09] ⛓ Spinning up manager process...
[14:31:07] 🚀 Starting up repository: 'project-fish', branch: 'main', main module: 'Gene5.py'
[14:31:07] 🐙 Cloning repository...
[14:31:08] 🐙 Cloning into '/mount/src/project-fish'...

[14:31:08] 🐙 Cloned repository!
[14:31:08] 🐙 Pulling code changes from Github...
[14:31:09] 📦 Processing dependencies...

──────────────────────────────────────── uv ───────────────────────────────────────────

Using uv pip install.
Using Python 3.13.5 environment at /home/adminuser/venv
Resolved 66 packages in 2.19s
Prepared 66 packages in 4.17s
Installed 66 packages in 155ms
 + altair==5.5.0
 + asttokens==3.0.0
 + attrs==25.3.0
 + biopython==1.85
 + blinker==1.9.0
 + cachetools==6.1.0
 + certifi==2025.7.14
 + charset-normalizer==3.4.2
 + click==8.2.1
 + comm==0.2.2
 + contourpy==1.3.2
 + cycler==0.12.1
 + decorator==5.2.1
 + executing==2.2.0
 + fonttools==4.59.0
 + fpdf==1.7.2
 + gitdb==4.0.12
 + gitpython==3.1.44
 + idna==3.10
 + ipython==9.4.0
 + ipython-pygments-lexers==1.1.1
 + ipywidgets==8.1.7
 + jedi==0.19.2
 + jinja2==3.1.6
 + jsonschema==4.25.0
 + jsonschema-specifications==2025.4.1
 + jupyterlab-widgets==3.0.15
 + kiwisolver==1.4.8
 + markupsafe==3.0.2
 + matplotlib==3.10.3
 + matplotlib-inline==0.1.7
 + narwhals==1.48.0
 + numpy==2.3.1
 + packaging==25.0
 + pandas==2.3.1
 + parso==0.8.4
 + pexpect==4.9.0
 + pillow==11.3.0
 + prompt-toolkit==3.0.51
 + protobuf==6.31.1
 + ptyprocess==0.7.0
 + pure-eval==0.2.3
 + pyarrow==21.0.0
 + pydeck==0.9.1
 + pygments==2.19.2
 + pyparsing==3.2.3
 + python-dateutil==2.9.0.post0
 + pytz==2025.2
 + referencing==0.36.2
 + requests==2.32.4
 + rpds-py==0.26.0
 + seaborn==0.13.2
 + six==1.17.0
 + smmap==5.0.2
 + stack-data==0.6.3
 + streamlit==1.47.0
 + tenacity==9.1.2
 + toml==0.10.2
 + tornado==6.5.1
 + traitlets==5.14.3
 + typing-extensions==4.14.1
 + tzdata==2025.2
 + urllib3==2.5.0
 + watchdog==6.0.0
 + wcwidth==0.2.13
 + widgetsnbextension==4.0.14
Checking if Streamlit is installed
Found Streamlit version 1.47.0 in the environment
Installing rich for an improved exception logging
Using uv pip install.
Using Python 3.13.5 environment at /home/adminuser/venv
Resolved [2025-07-22 14:31:18.318023] 4 packages in 164ms
Prepared [2025-07-22 14:31:18.373097] 3 packages in 55ms
Installed 3 packages in 8ms
 + markdown-it-py==3.0.0
 + mdurl==0.1.2
 + rich==14.0.0

────────────────────────────────────────────────────────────────────────────────────────

[14:31:19] 🐍 Python dependencies were installed from /mount/src/project-fish/requirements.txt using uv.
Check if streamlit is installed
Streamlit is already installed
[14:31:21] 📦 Processed dependencies!



[14:31:25] 🐙 Pulling code changes from Github...
[14:31:25] 📦 Processing dependencies...
[14:31:25] 📦 Processed dependencies!
[14:31:27] 🔄 Updated app!
/mount/src/project-fish/Gene5.py:208: FutureWarning: 

Passing `palette` without assigning `hue` is deprecated and will be removed in v0.14.0. Assign the `x` variable to `hue` and set `legend=False` for the same effect.

  sns.barplot(x='Accession', y='Length', data=df, palette='viridis')
VBox(children=(Checkbox(value=False, description='Display Table'), Button(description='Download CSV', style=ButtonStyle()), Button(description='Show Plot', style=ButtonStyle()), Output()))
VBox(children=(Checkbox(value=False, description='Display Table'), Button(description='Download CSV', style=ButtonStyle()), Button(description='Download GFF3', style=ButtonStyle()), Button(description='Show Plot', style=ButtonStyle()), Output()))
CSV report generated at: p4a_gene_report.csv
Download your CSV report from: p4a_gene_report.csv
     Accession Gene  Start  End  Strand Scaffold/Chromosome
0   MN316643.1  p4a      0  357      -1            MN316643
1   MW915530.1  p4a      0  451      -1            MW915530
2   MW660860.1  p4a      0  440      -1            MW660860
3   MW660859.1  p4a      0  413      -1            MW660859
4   MW660858.1  p4a      0  414      -1            MW660858
5   MW660857.1  p4a      0  414      -1            MW660857
6   MW660856.1  p4a      0  453      -1            MW660856
7   MK990715.1  P4a      0  478      -1            MK990715
8   MK990714.1  P4a      0  478      -1            MK990714
9   MK990713.1  P4a      0  478      -1            MK990713
10  MK990712.1  P4a      0  478      -1            MK990712
11  MK990711.1  P4a      0  478      -1            MK990711
12  MK990709.1  P4a      0  478      -1            MK990709
/mount/src/project-fish/p4a_gene_coordinates.csv
/mount/src/project-fish/gene_lengths_plot.png
────────────────────── Traceback (most recent call last) ───────────────────────
  /home/adminuser/venv/lib/python3.13/site-packages/streamlit/runtime/scriptru  
  nner/exec_code.py:128 in exec_func_with_error_handling                        
                                                                                
  /home/adminuser/venv/lib/python3.13/site-packages/streamlit/runtime/scriptru  
  nner/script_runner.py:669 in code_to_exec                                     
                                                                                
  /mount/src/project-fish/Gene5.py:368 in <module>                              
                                                                                
    365 import matplotlib.pyplot as plt                                         
    366 import streamlit as st                                                  
    367 ...                                                                     
  ❱ 368 st.download_button(...)  # Now it's safe to use                         
    369                                                                         
    370 import os                                                               
    371                                                                         
                                                                                
  /home/adminuser/venv/lib/python3.13/site-packages/streamlit/runtime/metrics_  
  util.py:443 in wrapped_func                                                   
────────────────────────────────────────────────────────────────────────────────
TypeError: ButtonMixin.download_button() missing 1 required positional argument:
'data'
