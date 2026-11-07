"""
SeqAlign Web : Streamlit Interactive Bioinformatics Workbench
============================================================

An interactive web application providing both in-browser pairwise sequence
alignment (Needleman-Wunsch & Smith-Waterman) and cloud-based Multiple
Sequence Alignment (MSA) powered by the EMBL-EBI Clustal Omega REST API.

Tabs:
1. Global Alignment (Needleman-Wunsch) - Pure Python with interactive DP Matrix.
2. Local Alignment (Smith-Waterman) - Pure Python with conserved motif highlighting.
3. Multiple Sequence Alignment - Integration with EMBL-EBI Clustal Omega web service.

Author: Azheem Syahmee
License: MIT
"""

import sys
import os
import time
from typing import Optional

# Ensure package root is in path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, "../.."))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

try:
    import streamlit as st
    import requests
    import pandas as pd
except ImportError as e:
    print(f"Error importing web dependencies: {e}")
    print("Please install via: pip install streamlit requests pandas")

from src.algorithms.needleman_wunsch import needleman_wunsch, format_matrix as format_nw_matrix
from src.algorithms.smith_waterman import smith_waterman, format_matrix as format_sw_matrix

# EMBL-EBI Clustal Omega REST API Endpoint
EBI_CLUSTALO_URL = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/"


def submit_clustalo_job(sequences: str, email: str, stype: str = "protein") -> Optional[str]:
    """
    Submit a multiple sequence alignment job to the EMBL-EBI Clustal Omega REST API.

    Args:
        sequences: Raw sequences in FASTA format.
        email: Valid user email address (required by EBI Fair Use policy).
        stype: Sequence type ('protein' or 'dna').

    Returns:
        Job ID string on success, or None on failure.
    """
    try:
        response = requests.post(
            EBI_CLUSTALO_URL + "run/",
            data={
                "email": email.strip(),
                "sequence": sequences.strip(),
                "stype": stype,
            },
            timeout=30,
        )
        if response.status_code == 200:
            return response.text.strip()
        else:
            st.error(f"EBI API Error (Status {response.status_code}): {response.text}")
            return None
    except requests.exceptions.RequestException as exc:
        st.error(f"Network request failed: {exc}")
        return None


def check_clustalo_status(job_id: str) -> Optional[str]:
    """Check status of submitted EBI job ('RUNNING', 'FINISHED', 'ERROR', 'NOT_FOUND')."""
    try:
        response = requests.get(f"{EBI_CLUSTALO_URL}status/{job_id}", timeout=15)
        return response.text.strip() if response.status_code == 200 else None
    except requests.exceptions.RequestException:
        return None


def get_clustalo_results(job_id: str, out_format: str = "aln-clustal") -> Optional[str]:
    """Retrieve completed multiple sequence alignment in specified format."""
    try:
        response = requests.get(f"{EBI_CLUSTALO_URL}result/{job_id}/{out_format}", timeout=30)
        return response.text if response.status_code == 200 else None
    except requests.exceptions.RequestException as exc:
        st.error(f"Failed to retrieve results: {exc}")
        return None


def main():
    st.set_page_config(
        page_title="SeqAlign - Bioinformatics Sequence Alignment",
        page_icon="🧬",
        layout="wide",
    )

    st.title("🧬 SeqAlign : Biological Sequence Alignment Workbench")
    st.markdown(
        "Interactive suite for pairwise alignment (implemented from scratch in pure Python) "
        "and multiple sequence alignment (interfacing with the EMBL-EBI Clustal Omega REST service)."
    )

    tab_nw, tab_sw, tab_msa = st.tabs([
        "🌐 Needleman-Wunsch (Global)",
        "🎯 Smith-Waterman (Local)",
        "🔬 Clustal Omega (Multiple Alignment)",
    ])

    # -------------------------------------------------------------
    # TAB 1: Needleman-Wunsch Global Alignment
    # -------------------------------------------------------------
    with tab_nw:
        st.header("Needleman-Wunsch Global Sequence Alignment")
        st.caption("Aligns two nucleotide or amino acid sequences end-to-end using dynamic programming.")

        col1, col2 = st.columns(2)
        with col1:
            nw_seq1 = st.text_input("Sequence 1", value="FKHMEDPLE", key="nw_s1")
        with col2:
            nw_seq2 = st.text_input("Sequence 2", value="FMDTPLNE", key="nw_s2")

        with st.expander("Scoring Scheme Parameters", expanded=False):
            pcol1, pcol2, pcol3 = st.columns(3)
            with pcol1:
                nw_match = st.number_input("Match Score", value=1, step=1, key="nw_m")
            with pcol2:
                nw_mismatch = st.number_input("Mismatch Penalty", value=-1, step=1, key="nw_mm")
            with pcol3:
                nw_gap = st.number_input("Gap Penalty", value=-2, step=1, key="nw_gap")

        if st.button("Run Global Alignment", type="primary", key="btn_nw"):
            if nw_seq1 and nw_seq2:
                res = needleman_wunsch(
                    nw_seq1, nw_seq2,
                    match_score=int(nw_match),
                    mismatch_penalty=int(nw_mismatch),
                    gap_penalty=int(nw_gap),
                )

                st.subheader("Alignment Result")
                metric_col1, metric_col2, metric_col3 = st.columns(3)
                metric_col1.metric("Alignment Score", res["score"])
                metric_col2.metric("Sequence Identity", f"{res['identity'] * 100:.1f}%")
                metric_col3.metric("Alignment Length", len(res["aligned_seq1"]))

                alignment_box = (
                    f"Seq1: {res['aligned_seq1']}\n"
                    f"      {res['similarity_str']}\n"
                    f"Seq2: {res['aligned_seq2']}"
                )
                st.code(alignment_box, language="text")

                with st.expander("Dynamic Programming Matrix Table", expanded=True):
                    # Convert matrix into pandas dataframe for rich table view
                    cols = ["-"] + list(nw_seq2.upper())
                    rows = ["-"] + list(nw_seq1.upper())
                    df = pd.DataFrame(res["matrix"], index=rows, columns=cols)
                    st.dataframe(df, use_container_width=True)
            else:
                st.warning("Please provide both sequence 1 and sequence 2.")

    # -------------------------------------------------------------
    # TAB 2: Smith-Waterman Local Alignment
    # -------------------------------------------------------------
    with tab_sw:
        st.header("Smith-Waterman Local Sequence Alignment")
        st.caption("Identifies the highest-scoring local conserved motif between two sequences with zero-floor DP.")

        col1, col2 = st.columns(2)
        with col1:
            sw_seq1 = st.text_input("Sequence 1 (Target)", value="ACACACTA", key="sw_s1")
        with col2:
            sw_seq2 = st.text_input("Sequence 2 (Query)", value="AGCACACA", key="sw_s2")

        with st.expander("Scoring Scheme Parameters", expanded=False):
            pcol1, pcol2, pcol3 = st.columns(3)
            with pcol1:
                sw_match = st.number_input("Match Score", value=3, step=1, key="sw_m")
            with pcol2:
                sw_mismatch = st.number_input("Mismatch Penalty", value=-1, step=1, key="sw_mm")
            with pcol3:
                sw_gap = st.number_input("Gap Penalty", value=-2, step=1, key="sw_gap")

        if st.button("Run Local Alignment", type="primary", key="btn_sw"):
            if sw_seq1 and sw_seq2:
                res = smith_waterman(
                    sw_seq1, sw_seq2,
                    match_score=int(sw_match),
                    mismatch_penalty=int(sw_mismatch),
                    gap_penalty=int(sw_gap),
                )

                if res["score"] == 0:
                    st.info("No positive local alignment found with the given scoring scheme (Score = 0).")
                else:
                    st.subheader("Conserved Local Motif")
                    mcol1, mcol2, mcol3 = st.columns(3)
                    mcol1.metric("Max Local Score", res["score"])
                    mcol2.metric("Motif Identity", f"{res['identity'] * 100:.1f}%")
                    mcol3.metric("Motif Length", len(res["aligned_seq1"]))

                    st.markdown(
                        f"**Motif Coordinates:** Seq1 `[{res['start1']}..{res['end1']}]` | "
                        f"Seq2 `[{res['start2']}..{res['end2']}]`"
                    )
                    alignment_box = (
                        f"Seq1 [{res['start1']}]: {res['aligned_seq1']} [{res['end1']}]\n"
                        f"         {res['similarity_str']}\n"
                        f"Seq2 [{res['start2']}]: {res['aligned_seq2']} [{res['end2']}]"
                    )
                    st.code(alignment_box, language="text")

                    with st.expander("Dynamic Programming Matrix Table", expanded=True):
                        cols = ["-"] + list(sw_seq2.upper())
                        rows = ["-"] + list(sw_seq1.upper())
                        df = pd.DataFrame(res["matrix"], index=rows, columns=cols)
                        st.dataframe(df, use_container_width=True)
            else:
                st.warning("Please provide both sequence 1 and sequence 2.")

    # -------------------------------------------------------------
    # TAB 3: Multiple Sequence Alignment (EMBL-EBI Clustal Omega)
    # -------------------------------------------------------------
    with tab_msa:
        st.header("Multiple Sequence Alignment (EMBL-EBI Clustal Omega)")
        st.caption(
            "Submits multi-sequence FASTA inputs to the European Bioinformatics Institute (EMBL-EBI) "
            "Clustal Omega REST API service."
        )

        sample_fasta = (
            ">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha\n"
            "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR\n"
            ">sp|P68871|HBB_HUMAN Hemoglobin subunit beta\n"
            "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH\n"
            ">sp|P02042|HBD_HUMAN Hemoglobin subunit delta\n"
            "MVHLTPEEKTAVNALWGKVNVDAVGGEALGRLLVVYPWTQRFFESFGDLSSPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFSQLSELHCDKLHVDPENFRLLGNVLVCVLARNFGKEFTPQMQAAYQKVVAGVANALAHKYH\n"
        )

        user_email = st.text_input(
            "Contact Email (Required by EMBL-EBI Fair Use Policy):",
            placeholder="researcher@university.edu",
            help="EMBL-EBI requires a valid email for job notification and server capacity tracking.",
        )
        seq_type = st.radio("Sequence Type:", ["protein", "dna"], horizontal=True)
        fasta_input = st.text_area("Sequences in FASTA Format (minimum 3 sequences):", value=sample_fasta, height=220)

        if st.button("Submit to Clustal Omega", type="primary", key="btn_msa"):
            if not user_email or "@" not in user_email:
                st.error("Please enter a valid email address as required by EMBL-EBI.")
            elif not fasta_input.strip().startswith(">"):
                st.error("Input must be in standard FASTA format starting with '>' header lines.")
            else:
                with st.spinner("Submitting job to EMBL-EBI Clustal Omega..."):
                    job_id = submit_clustalo_job(fasta_input, user_email, stype=seq_type)

                if job_id:
                    st.success(f"Job submitted successfully! Tracking Job ID: `{job_id}`")
                    status_placeholder = st.empty()
                    progress_bar = st.progress(10)

                    poll_count = 0
                    max_polls = 30
                    status = "RUNNING"

                    while status in ("RUNNING", "QUEUED") and poll_count < max_polls:
                        time.sleep(4)
                        poll_count += 1
                        status = check_clustalo_status(job_id)
                        progress_val = min(20 + poll_count * 5, 90)
                        progress_bar.progress(progress_val)
                        status_placeholder.info(f"Checking job status... State: **{status}** (Poll {poll_count}/{max_polls})")

                    if status == "FINISHED":
                        progress_bar.progress(100)
                        status_placeholder.success("Alignment calculation completed!")
                        results = get_clustalo_results(job_id, out_format="aln-clustal")
                        if results:
                            st.subheader("Clustal Omega Alignment Output")
                            st.code(results, language="text")
                            st.download_button(
                                "Download Clustal Alignment",
                                data=results,
                                file_name=f"clustalo_alignment_{job_id}.aln",
                                mime="text/plain",
                            )
                    elif status == "ERROR":
                        status_placeholder.error("Job encountered a server error on EMBL-EBI.")
                    else:
                        status_placeholder.warning("Job is still running. You can check results later using the Job ID.")


if __name__ == "__main__":
    main()
