import streamlit as st
from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt
from fpdf import FPDF
import os
from io import StringIO

# Configure Entrez email
Entrez.email = "your_email@example.com"  # Replace with your actual email

st.set_page_config(page_title="Fish Protein Explorer", layout="centered")
st.title("🐟 Fish Protein Sequence Explorer")

st.markdown("Use the dropdowns to select parameters, then click Search.")

# Dropdowns
species = st.selectbox("Select Fish Species", [
    "Common Carp (Cyprinus carpio)", 
    "Grass Carp (Ctenopharyngodon idella)",
    "Silver Carp (Hypophthalmichthys molitrix)", 
    "Catla carp (Cyprinus catla)",
    "Bighead Carp (Hypophthalmichthys nobilis)"
])

organ = st.selectbox("Select Organ", ["Gill", "Kidney", "Liver", "Spleen", "Heart", "Intestine"])
gene = st.selectbox("Select Gene", ["p4a", "vp8", "ORF16", "CyHV1", "CEV_F1", "CEV_R1", "CEV_F2"])
disease = st.selectbox("Select Disease", ["Carp Edema Virus", "koi herpesvirus", "cyprinid herpesvirus-2"])

protein_search_term = st.selectbox("Select Protein", [
    "Carp edema virus isolate Agnes_661xx2018 P4a core protein (p4a) gene, partial cds",
    "Carp edema virus isolate IR-UT-ZGNARM 40020 4a (p4a) gene, partial cds",
    "Carp edema virus isolate IR-UT-ZGNARM 9360 4a (p4a) gene, partial cds",
    "Carp edema virus isolate IR UT-ZGNARM 9429 4a (p4a) gene, partial cds",
    "Carp edema virus isolate IR-UT-ZGNARM 9171 4a (p4a) gene, partial cds",
    "Carp edema virus isolate IR-UT-ZGNARM 9325 4a (p4a) gene, partial cds",
    "Carp edema virus isolate IR-UT-ZGNARM 9299 4a (p4a) gene, partial cds",
    "Carp edema virus 4a protein (P4a) gene, partial cds",
    "Carp edema virus isolate F238 4a protein (P4a) gene, partial cds",
    "Carp edema virus isolate LE2 4a protein (P4a) gene, partial cds",
    "Carp edema virus isolate LE1 4a protein (P4a) gene, partial cds",
    "Carp edema virus isolate F129 4a protein (P4a) gene, partial cds",
    "Carp edema virus isolate F14 4a protein (P4a) gene, partial cds"
])

# Options
show_text = st.checkbox("Show GenBank Record")
show_plot = st.checkbox("Show Protein Length Plot")
download_pdf = st.checkbox("Download PDF Report")

# Search and display
if st.button("🔍 Search Protein"):
    with st.spinner("Searching NCBI and fetching data..."):
        handle = Entrez.esearch(db="protein", term=protein_search_term, retmax=1)
        result = Entrez.read(handle)
        handle.close()

        if not result['IdList']:
            st.error("❌ No matching protein record found.")
        else:
            protein_id = result['IdList'][0]
            handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text")
            protein_text = handle.read()
            handle.close()

            seq_record = SeqIO.read(StringIO(protein_text), "genbank")
            seq_len = len(seq_record.seq)

            if show_text:
                st.subheader("📄 GenBank Record (truncated):")
                st.text(protein_text[:1200] + "\n...[truncated]...")

            if show_plot:
                st.subheader("🧬 Protein Length Plot")
                fig, ax = plt.subplots(figsize=(5, 1.5))
                ax.barh(["Selected Protein"], [seq_len], color='skyblue')
                ax.set_xlabel("Amino Acid Length")
                st.pyplot(fig)

            if download_pdf:
                pdf = FPDF()
                pdf.add_page()
                pdf.set_font("Arial", "B", 14)
                pdf.cell(200, 10, "Selected Protein Report", ln=True)
                pdf.set_font("Arial", size=10)
                pdf.multi_cell(0, 5, protein_text[:2500] + "\n...[truncated]...")
                pdf_path = "/tmp/protein_report.pdf"
                pdf.output(pdf_path)
                with open(pdf_path, "rb") as f:
                    st.download_button("📄 Download PDF Report", f, file_name="protein_report.pdf")