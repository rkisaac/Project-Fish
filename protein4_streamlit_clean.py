import streamlit as st
from Bio import Entrez, SeqIO
from fpdf import FPDF
import matplotlib.pyplot as plt
from io import StringIO
import base64
import os

# Email for NCBI Entrez
Entrez.email = "your_email@example.com"

st.set_page_config(page_title="Protein Explorer", layout="centered")
st.title("🧬 Fish Protein Sequence Explorer")

# Dropdown options
fish_species = [
    "Common Carp (Cyprinus carpio)", "Grass Carp (Ctenopharyngodon idella)",
    "Silver Carp (Hypophthalmichthys molitrix)", "Catla carp (Cyprinus catla)",
    "Bighead Carp (Hypophthalmichthys nobilis)"
]

fish_organs = ["Gill", "Kidney", "Liver", "Spleen", "Heart", "Intestine"]
genes = ["p4a", "vp8", "ORF16", "CyHV1", "CEV_F1", "CEV_R1", "CEV_F2"]
diseases = ["Carp Edema Virus", "koi herpesvirus", "cyprinid herpesvirus-2"]

protein_list = [
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
]

# UI Inputs
species = st.selectbox("🐟 Fish Species", fish_species)
organ = st.selectbox("🧫 Organ", fish_organs)
gene = st.selectbox("🧬 Gene", genes)
disease = st.selectbox("🦠 Disease", diseases)
protein_query = st.selectbox("🔎 Select Protein Query", protein_list)

if st.button("🔍 Search Protein"):
    if species != "Common Carp (Cyprinus carpio)" or organ != "Gill" or gene != "p4a" or disease != "Carp Edema Virus":
        st.warning("Please select: Common Carp, Gill, p4a, and Carp Edema Virus to proceed.")
    else:
        st.info(f"🔍 Searching NCBI for: {protein_query}")
        handle = Entrez.esearch(db="protein", term=protein_query, retmax=1)
        result = Entrez.read(handle)
        handle.close()

        if not result["IdList"]:
            st.error("❌ No matching protein found.")
        else:
            protein_id = result["IdList"][0]
            st.success(f"✅ Found Protein ID: {protein_id}")
            handle = Entrez.efetch(db="protein", id=protein_id, rettype="gb", retmode="text")
            protein_text = handle.read()
            handle.close()

            # Save GenBank text
            genbank_path = "protein_record.txt"
            with open(genbank_path, "w") as f:
                f.write(protein_text)

            # Load and visualize sequence
            record = SeqIO.read(StringIO(protein_text), "genbank")
            seq_len = len(record.seq)

            fig, ax = plt.subplots(figsize=(6, 1.5))
            ax.barh(["Selected Protein"], [seq_len], color='skyblue')
            ax.set_xlabel("Amino Acid Length")
            ax.set_title("Protein Length")
            st.pyplot(fig)

            # Show GenBank text
            with st.expander("📄 View GenBank Record"):
                st.text(protein_text[:2000] + "\n...[truncated]...")

            # PDF Report
            pdf_path = "protein_report.pdf"
            pdf = FPDF()
            pdf.add_page()
            pdf.set_font("Arial", "B", 14)
            pdf.cell(200, 10, "Protein Report", ln=1)
            pdf.set_font("Arial", size=10)
            pdf.multi_cell(0, 5, protein_text[:2500] + "\n...[truncated]...")
            pdf.output(pdf_path)

            # Download PDF
            with open(pdf_path, "rb") as f:
                st.download_button("📥 Download PDF Report", f, file_name="protein_report.pdf", mime="application/pdf")