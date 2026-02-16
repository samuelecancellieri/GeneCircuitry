"""
TRNspot Reporting Module
========================

This module generates comprehensive HTML and PDF reports from pipeline results,
including input data summary, configuration settings, operations performed,
and interactive plot galleries.
"""

from .generator import (
    ReportGenerator,
    generate_report,
    generate_stratified_report,
    generate_html_report,
    generate_pdf_report,
)
from .sections import (
    create_data_summary_section,
    create_settings_section,
    create_qc_section,
    create_preprocessing_section,
    create_clustering_section,
    create_celloracle_section,
    create_hotspot_section,
    create_operations_log_section,
    create_plot_gallery_section,
    create_grn_deep_analysis_section,
    create_stratified_clustering_section,
    create_stratified_celloracle_section,
    create_stratified_hotspot_section,
    create_stratification_summary_section,
)

__all__ = [
    # Main generator
    "ReportGenerator",
    "generate_report",
    "generate_stratified_report",
    "generate_html_report",
    "generate_pdf_report",
    # Section builders
    "create_data_summary_section",
    "create_settings_section",
    "create_qc_section",
    "create_preprocessing_section",
    "create_clustering_section",
    "create_celloracle_section",
    "create_hotspot_section",
    "create_operations_log_section",
    "create_plot_gallery_section",
    "create_grn_deep_analysis_section",
    # Stratified section builders
    "create_stratified_clustering_section",
    "create_stratified_celloracle_section",
    "create_stratified_hotspot_section",
    "create_stratification_summary_section",
]
